"""Do earthquakes change the instruments? Superposed-epoch analysis on every station.

For every station and every catalogue earthquake (USGS, M >= 5.85) during its
record, a shaking proxy is computed (magnitude minus 1.66 log10 of the epicentral
distance in km, i.e. a surface-wave magnitude-distance scaling; larger = stronger).
Around each of the strongest events per station, two changes are measured:

* tilt: median hourly tilt angle (``03_tilt``) in the 48 h after (from +3 h)
  minus the 48 h before (to -1 h), and the circular change of azimuth;
* \\|Z/P\\|: median log10 \\|S_pz\\|/S_pp of coherent windows (gamma > 0.8,
  ``13_audit/stage_c/.../window_table``) in the 4 days after minus before.

The null distribution of each statistic comes from 300 random times per station
at least 3 days from any M >= 6.5 earthquake. Output:
``13_audit/research/earthquake_epochs.csv`` (station-event pairs with robust
z-scores against the station's own null) and a summary by shaking class.
"""

import csv
import glob
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from compy import anomaly  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
OUT = ARCHIVE / "13_audit/research"
H = 3600.0


def ts(s):
    return datetime.fromisoformat(s.replace("Z", "+00:00")).timestamp()


def dist_km(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = (
        np.sin((la2 - la1) / 2) ** 2
        + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    )
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def change(t, x, t0, before, after, circular=False):
    b = x[(t >= t0 - before) & (t < t0 - 1 * H)]
    a = x[(t >= t0 + 3 * H) & (t < t0 + after)]
    if (
        len(b) < 0.5 * before / H / (4 if before > 3 * 86400 else 1) * 0.25
        or len(a) < 5
        or len(b) < 5
    ):
        return np.nan
    if circular:
        mb = np.degrees(np.angle(np.mean(np.exp(1j * np.radians(b)))))
        ma = np.degrees(np.angle(np.mean(np.exp(1j * np.radians(a)))))
        return (ma - mb + 180) % 360 - 180
    return float(np.median(a) - np.median(b))


def one(task):
    sid, lat, lon, quakes = task
    net, sta = sid.split(".")
    series = {}
    path = ARCHIVE / "03_tilt" / net / f"{sid}_tilt_hourly.csv"
    if path.exists():
        t, ang, az = [], [], []
        for r in csv.DictReader(open(path)):
            try:
                t.append(ts(r["centre_utc"]))
                ang.append(float(r["angle_deg"]))
                az.append(float(r["azimuth_deg"]))
            except ValueError:
                pass
        series["tilt_angle"] = (np.array(t), np.array(ang), 2 * 86400, False)
        series["tilt_azimuth"] = (np.array(t), np.array(az), 2 * 86400, True)
    tabs = glob.glob(
        str(
            ARCHIVE
            / f"13_audit/stage_c/wide_inclusive/{net}/results/{sta}/{sta}_window_table_*.npz"
        )
    )
    if tabs:
        w = np.load(tabs[0])
        ok = w["coherence"] > 0.8
        series["zp_level"] = (w["time"][ok], w["log10_tf"][ok], 4 * 86400, False)
    if not series:
        return []
    t_all = np.concatenate([s[0] for s in series.values()])
    t0, t1 = t_all.min(), t_all.max()
    big = np.array([q[0] for q in quakes if q[1] >= 6.5])
    rng = np.random.default_rng(abs(hash(sid)) % 2**32)
    cand = rng.uniform(t0 + 5 * 86400, t1 - 5 * 86400, 2000)
    if len(big):
        cand = cand[np.min(np.abs(cand[:, None] - big[None, :]), axis=1) > 3 * 86400]
    cand = cand[:300]
    null = {
        k: np.array([change(s[0], s[1], c, s[2], s[2], s[3]) for c in cand])
        for k, s in series.items()
    }
    rows = []
    for qt, mag, qlat, qlon in quakes:
        if not t0 + 2 * 86400 < qt < t1 - 2 * 86400:
            continue
        d = dist_km(lat, lon, qlat, qlon)
        shake = mag - 1.66 * np.log10(max(d, 10.0))
        if shake < 1.0:  # e.g. M7 at 6000 km gives 0.7; M6 at 100 km gives 2.7
            continue
        row = {"sid": sid, "time": datetime.utcfromtimestamp(qt).isoformat()[:16], "mag": mag,
               "dist_km": round(float(d)), "shaking": round(float(shake), 2)}  # fmt: skip
        for k, s in series.items():
            v = change(s[0], s[1], qt, s[2], s[2], s[3])
            n = null[k][np.isfinite(null[k])]
            if np.isfinite(v) and len(n) > 30:
                row[k] = v
                row[f"{k}_z"] = (v - np.median(n)) / max(anomaly.mad_sigma(n), 1e-9)
                row[f"{k}_p"] = float(
                    np.mean(np.abs(n - np.median(n)) >= abs(v - np.median(n)))
                )
        rows.append(row)
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    st = {
        f"{r['network']}.{r['station']}": r
        for r in csv.DictReader(open(ARCHIVE / "00_inventory/stations.csv"))
    }
    plan = [
        r["sid"]
        for r in csv.DictReader(
            open(ARCHIVE / "12_inversion_all_stations/tables/station_plan.csv")
        )
    ]
    cat = [(ts(c["time"]), float(c["magnitude"]), float(c["latitude"]), float(c["longitude"]))
           for c in csv.DictReader(open(ARCHIVE / "00_inventory/catalog_usgs_M5.85.csv"))]  # fmt: skip
    tasks = [
        (sid, float(st[sid]["latitude"]), float(st[sid]["longitude"]), cat)
        for sid in plan
        if sid in st
    ]
    with ProcessPoolExecutor(int(os.environ.get("WORKERS", 8))) as ex:
        rows = [r for rs in ex.map(one, tasks) for r in rs]
    keys = sorted(
        {k for r in rows for k in r},
        key=lambda k: (k not in ("sid", "time", "mag", "dist_km", "shaking"), k),
    )
    with open(OUT / "earthquake_epochs.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} station-event pairs")


if __name__ == "__main__":
    main()
