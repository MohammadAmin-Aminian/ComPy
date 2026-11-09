"""Two independent measurements of the DPG leak corner: compliance phase and the M2 tide.

(a) Phase. The residual Z/P phase in the compliance band is fitted with one factor k
    that scales every metadata pole below ``LEAK_HZ`` of the pressure channel (the
    leak and coupling corners): H_true(f) = H_meta(f; poles * k). Residual =
    phase(H_meta) - phase(H_true). k < 1 means a longer true time constant.
(b) Tide (Doran et al. 2019). A DPG still records the tide, attenuated by its leak.
    The M2 amplitude of the decimated counts, deconvolved with the metadata response
    at the M2 frequency and multiplied by the station's calibration gain, is divided
    by the EOT20 bottom-pressure tide (ocean tide x rho g). If the metadata were
    right, R_M2 = 1; with the corners scaled by k, R_M2 = |H_true/H_meta| at M2,
    which (a) predicts without using the tide.

Output: ``<OUT>/dpg_tide/dpg_tide.csv``. Needs ObsPy and xarray (EOT20 from the
pressure-atlas cache, ``EOT20_DIR``).
"""

import csv
import glob
import json
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from compy import response_check as rc  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
OUT = Path(os.environ.get("COMPY_AUDIT_OUT", ARCHIVE / "13_audit")) / "dpg_tide"
EOT20 = Path(
    os.environ.get(
        "EOT20_DIR",
        Path.home()
        / "Desktop/BRUIT-FM/ComPy-pressure-atlas/pressure_atlas/cache/eot20",
    )
)
LEAK_HZ = 0.05
SPEED = {
    "M2": 28.9841042,
    "S2": 30.0,
    "N2": 28.4397295,
    "K1": 15.0410686,
    "O1": 13.9430356,
}  # deg/h
F_M2 = rc.F_M2
RHO_G = 1035.0 * 9.81  # Pa per m, as in the pressure atlas tide check
warnings.filterwarnings("ignore")


# The response model and the harmonic fit live in compy.response_check (tested).
paz = rc.paz_stages


def analog(stages, f, k=1.0, mode="all"):
    return rc.analog_response(stages, f, k, mode=mode, leak_hz=LEAK_HZ)


def fit_k(stages, f, residual, mode="all"):
    fit = rc.fit_leak_factor(stages, f, residual, mode=mode, leak_hz=LEAK_HZ)
    return fit["k"], fit["rms_deg"], fit["k_low"], fit["k_high"]


def model_amp(lat, lon, cons, kind):
    import xarray as xr

    ds = xr.open_dataset(
        glob.glob(str(EOT20 / kind / f"{cons}_*"))[0], engine="netcdf4"
    )
    a = ds["amplitude"]
    lon = lon % 360
    x = a.sel(lat=slice(lat - 0.5, lat + 0.5), lon=slice(lon - 0.5, lon + 0.5))
    v = x.values.astype(float)
    la, lo = np.meshgrid(x.lat.values, x.lon.values, indexing="ij")
    ok = np.isfinite(v)
    if not ok.any():
        return np.nan
    d = np.hypot(la - lat, (lo - lon) * np.cos(np.radians(lat)))
    w = np.exp(-((d[ok] / 0.125) ** 2))
    return float((v[ok] * w).sum() / w.sum()) / 100.0  # cm -> m


def nodal(mid_epoch):
    T = (mid_epoch - 946728000.0) / (86400 * 36525)
    N = np.radians((125.04 - 1934.136 * T) % 360)
    return {"M2": 1 - 0.037 * np.cos(N), "S2": 1.0, "N2": 1 - 0.037 * np.cos(N),
            "K1": 1.006 + 0.115 * np.cos(N), "O1": 1.009 + 0.187 * np.cos(N)}  # fmt: skip


def tide_counts(net, sta, cha):
    """Hourly means of the counts, then ``compy.response_check.tidal_harmonics``."""
    import obspy

    ts, vs = [], []
    for path in sorted(
        glob.glob(str(ARCHIVE / f"01_preprocessed/{net}/{sta}/counts/*.mseed"))
    ):
        try:
            st = obspy.read(path).select(channel=cha)
        except Exception:
            continue
        for tr in st:
            n = int(round(3600 * tr.stats.sampling_rate))
            if tr.stats.npts < n:
                continue
            m = tr.stats.npts // n
            vs.append(tr.data[: m * n].reshape(m, n).astype(float).mean(axis=1))
            ts.append(tr.stats.starttime.timestamp + (np.arange(m) + 0.5) * 3600)
    if not vs:
        return None
    t, v = np.concatenate(ts), np.concatenate(vs)
    if (t.max() - t.min()) / 86400 < 30:
        return None
    out = rc.tidal_harmonics(t, v)
    nod = nodal(float(np.mean(t)))
    return {
        "days": out["days"],
        "amp_counts": {c: a / nod[c] for c, a in out["amplitude"].items()},
        "noise_counts": out["noise"],
    }


def one(task):
    net, sta, row = task
    from obspy import UTCDateTime, read_inventory

    out = {"sid": f"{net}.{sta}"}
    try:
        inv = read_inventory(
            str(ARCHIVE / "00_inventory/stationxml" / f"{net}.xml")
        ).select(station=sta)
        t = UTCDateTime(float(row["t_mid"]))
        chans = [c for n in inv for s in n for c in s if c.code == row["p_channel"]]
        inside = [
            c
            for c in chans
            if c.start_date <= t and (c.end_date is None or c.end_date >= t)
        ]
        ch = (inside or chans)[0]
        stages = paz(ch)
        f = np.array(json.loads(row["f"]))
        res = np.array(json.loads(row["residual"]))
        k, rms, klo, khi = fit_k(stages, f, res)
        pred = lambda kk, mode="all": float(
            np.abs(analog(stages, F_M2, kk, mode) / analog(stages, F_M2))[0]
        )  # noqa: E731
        k1, rms1, k1lo, k1hi = fit_k(stages, f, res, "lowest")
        out.update(k_lowest=k1, k_lowest_lo=k1lo, k_lowest_hi=k1hi, rms_k_lowest=rms1,
                   R_pred_lowest=pred(k1, "lowest"), R_pred_lowest_lo=pred(k1lo, "lowest"),
                   R_pred_lowest_hi=pred(k1hi, "lowest"))  # fmt: skip
        out.update(k=k, k_lo=klo, k_hi=khi, rms_k=rms, R_pred=pred(k), R_pred_lo=pred(klo), R_pred_hi=pred(khi),
                   leak_corners_mhz=";".join(f"{1e3 * abs(p) / (2 * np.pi if rad else 1):.3g}" for _, ps, rad in stages
                                             for p in ps if 0 < abs(p) / (2 * np.pi if rad else 1) < LEAK_HZ))  # fmt: skip
        tide = tide_counts(net, sta, row["p_channel"])
        if tide is None:
            return out | {"tide_status": "no_counts"}
        # evalresp underflows to 0 at 22 µHz for some responses: take the full response at
        # 0.1 Hz and scale it with the analog part (digital stages are flat at these periods).
        f_n = 0.1
        r_n = ch.response.get_evalresp_response_for_frequencies(
            np.array([f_n]), output="DEF"
        )[0]
        resp = abs(r_n) * abs(analog(stages, F_M2)[0] / analog(stages, f_n)[0])
        gain = float(row["gain"])
        p_obs = tide["amp_counts"]["M2"] / abs(resp) * gain
        p_noise = tide["noise_counts"] / abs(resp) * gain
        p_mod = model_amp(float(row["lat"]), float(row["lon"]), "M2", "ocean") * RHO_G
        out.update(tide_status="ok", tide_days=tide["days"], p_obs_M2=p_obs, p_noise=p_noise, p_mod_M2=p_mod,
                   snr=tide["amp_counts"]["M2"] / tide["noise_counts"], R_M2=p_obs / p_mod,
                   R_K1=tide["amp_counts"]["K1"] / abs(resp) * gain / (model_amp(float(row["lat"]), float(row["lon"]), "K1", "ocean") * RHO_G),
                   resp_M2_counts_per_pa=abs(resp))  # fmt: skip
    except Exception as exc:  # noqa: BLE001
        out["tide_status"] = f"error: {exc}"
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    phase = {
        r["sid"]: r
        for r in csv.DictReader(open(ARCHIVE / "13_audit/phase/phase_metadata.csv"))
    }
    curves = json.loads((ARCHIVE / "13_audit/phase/phase_curves.json").read_text())
    stations = {
        f"{r['network']}.{r['station']}": r
        for r in csv.DictReader(open(ARCHIVE / "00_inventory/stations.csv"))
    }
    plan = {
        r["sid"]: r
        for r in csv.DictReader(
            open(ARCHIVE / "12_inversion_all_stations/tables/station_plan.csv")
        )
    }
    from compy.results import load_station

    tasks = []
    only = set(sys.argv[1:])
    for sid, p in phase.items():
        if p["pressure_type"] not in ("DPG", "Hydrophone") or (
            only and sid not in only
        ):
            continue
        if abs(float(p["res_lo"])) > 45 or abs(float(p["res_hi"])) > 45:
            continue
        net, sta = sid.split(".")
        folder = ARCHIVE / "06_compliance_wide" / net / "results" / sta
        if plan[sid]["epoch_end"]:
            folder = folder / "epoch"
        s = load_station(folder, sta, plan[sid]["wide_gain_used"] or None)
        row = {"p_channel": p["p_channel"], "t_mid": float(np.median(s.window_starts)), "gain": s.gain,
               "lat": stations[sid]["latitude"], "lon": stations[sid]["longitude"],
               "f": json.dumps(curves[sid]["f"]), "residual": json.dumps(curves[sid]["residual"])}  # fmt: skip
        tasks.append((net, sta, row))
    with ProcessPoolExecutor(int(os.environ.get("WORKERS", 4))) as ex:
        rows = list(ex.map(one, tasks))
    for r in rows:
        r.update(
            {
                k: phase[r["sid"]][k]
                for k in (
                    "network",
                    "pressure_type",
                    "sensor_family",
                    "depth_m",
                    "res_lo",
                    "res_hi",
                )
            }
        )
        r["gain_label"] = plan[r["sid"]]["wide_gain_used"]
    keys = sorted({k for r in rows for k in r})
    with open(OUT / "dpg_tide.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} stations")


if __name__ == "__main__":
    main()
