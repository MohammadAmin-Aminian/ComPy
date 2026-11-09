"""Shape of the DPG low-frequency response from four tidal constituents.

s03 compared M2 with EOT20 and found the deconvolved tide 8.6 dB too large at
the median DPG, with M2 and K1 ratios that differ. Here the ratio
R_c = (observed amplitude / |H_meta(f_c)| x gain) / (EOT20 ocean tide x rho g)
is measured for O1, K1, M2 and S2 at every DPG with a clear tide, and a power
law R ~ f^n is fitted across the four constituents.

* The metadata leak corner scaled by k (s03, both models): n = 0 (both
  responses are proportional to f far below the corner), R = 1/k.
* One more first-order high-pass in the true response than in the metadata
  (corner well above the tides): n = +1.
* A flatter true response (e.g. a parallel slow path): n < 0.

Output: ``13_audit/research/tide_constituents.csv``.
"""

import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))
import s03_dpg_tide as s03  # noqa: E402

ARCHIVE = s03.ARCHIVE
OUT = ARCHIVE / "13_audit/research"
CONS = ("O1", "K1", "M2", "S2")
FREQ = {c: s03.SPEED[c] / 360 / 3600 for c in CONS}


def one(task):
    sid, cha, t_mid, gain, lat, lon = task
    from obspy import UTCDateTime, read_inventory

    net, sta = sid.split(".")
    out = {"sid": sid}
    try:
        inv = read_inventory(
            str(ARCHIVE / "00_inventory/stationxml" / f"{net}.xml")
        ).select(station=sta)
        chans = [c for n in inv for s in n for c in s if c.code == cha]
        if t_mid is not None:  # pick the epoch covering the record when several exist
            t = UTCDateTime(t_mid)
            chans = [
                c
                for c in chans
                if c.start_date <= t and (c.end_date is None or c.end_date >= t)
            ] or chans
        ch = chans[0]
        stages = s03.paz(ch)
        tide = s03.tide_counts(net, sta, cha)
        if tide is None:
            return out
        r_n = abs(
            ch.response.get_evalresp_response_for_frequencies(
                np.array([0.1]), output="DEF"
            )[0]
        )
        for c in CONS:
            resp = r_n * abs(
                s03.analog(stages, FREQ[c])[0] / s03.analog(stages, 0.1)[0]
            )
            model = s03.model_amp(lat, lon, c, "ocean") * s03.RHO_G
            out[f"R_{c}"] = tide["amp_counts"][c] / resp * gain / model
            out[f"snr_{c}"] = tide["amp_counts"][c] / tide["noise_counts"]
        ok = [c for c in CONS if out.get(f"snr_{c}", 0) > 5 and out[f"R_{c}"] > 0]
        if len(ok) >= 3:
            x = np.log([FREQ[c] for c in ok])
            y = np.log([out[f"R_{c}"] for c in ok])
            n, a = np.polyfit(x, y, 1)
            out.update(
                n_slope=float(n), n_used=len(ok), rms_ln=float(np.std(y - (n * x + a)))
            )
    except Exception as exc:  # noqa: BLE001
        out["error"] = str(exc)
    return out


def main():
    rows = {
        r["sid"]: r
        for r in csv.DictReader(open(ARCHIVE / "13_audit/dpg_tide/dpg_tide.csv"))
    }
    phase = {
        r["sid"]: r
        for r in csv.DictReader(open(ARCHIVE / "13_audit/phase/phase_metadata.csv"))
    }
    st = {
        f"{r['network']}.{r['station']}": r
        for r in csv.DictReader(open(ARCHIVE / "00_inventory/stations.csv"))
    }
    stc = {
        r["sid"]: r
        for r in csv.DictReader(open(ARCHIVE / "13_audit/anomaly/stations.csv"))
    }
    tasks = []
    for sid, r in rows.items():
        if (
            r.get("tide_status") != "ok"
            or r["pressure_type"] != "DPG"
            or float(r["snr"]) < 20
        ):
            continue
        tasks.append((sid, phase[sid]["p_channel"], None, float(stc[sid]["gain"]),
                      float(st[sid]["latitude"]), float(st[sid]["longitude"])))  # fmt: skip
    with ProcessPoolExecutor(int(os.environ.get("WORKERS", 6))) as ex:
        out = list(ex.map(one, tasks))
    for o in out:
        o.update(
            {
                k: rows[o["sid"]][k]
                for k in ("network", "leak_corners_mhz", "k_lowest", "R_pred_lowest")
            }
        )
    keys = sorted({k for o in out for k in o})
    with open(OUT / "tide_constituents.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(out)
    n = np.array([o["n_slope"] for o in out if "n_slope" in o])
    print(
        f"{len(out)} DPGs, slope fitted at {len(n)}: median n = {np.median(n):.2f} (IQR {np.percentile(n, 25):.2f} to {np.percentile(n, 75):.2f})"
    )


if __name__ == "__main__":
    main()
