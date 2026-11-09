"""Which sensor stepped? The tide as a pressure-only reference, month by month.

The tide is a known, stationary signal recorded by the pressure channel only. Its
M2 amplitude in the raw counts, fitted in 30-day windows (15-day step; M2, S2, N2,
K1, O1 and a spline drift, ``compy.response_check.tidal_harmonics``), follows the
pressure sensitivity and nothing else. At a |Z/P| step (s05):

* a pressure-sensitivity change of x % changes M2 by x % and |Z/P| by -x %;
* a seismometer change leaves M2 flat.

YV.RR50 (vertical step on 2013-07-30, s04) is the control: M2 must stay flat.

Output: ``13_audit/research/monthly_tide.csv`` (per window) and
``tide_at_steps.csv`` (per step: M2 change and its scatter).
"""

import csv
import glob
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))
from compy import response_check as rc  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
OUT = ARCHIVE / "13_audit/research"
WIN_D, STEP_D = 30.0, 15.0


def hourly(net, sta, cha):
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
        return None, None
    t, v = np.concatenate(ts), np.concatenate(vs)
    o = np.argsort(t)
    return t[o], v[o]


def one(task):
    sid, cha = task
    net, sta = sid.split(".")
    t, v = hourly(net, sta, cha)
    rows = []
    if t is None:
        return sid, rows
    start = t[0]
    while start + WIN_D * 86400 <= t[-1] + 3600:
        sel = (t >= start) & (t < start + WIN_D * 86400)
        if sel.sum() >= 0.7 * WIN_D * 24:
            try:
                out = rc.tidal_harmonics(t[sel], v[sel])
                rows.append({"sid": sid, "t_mid": float(start + WIN_D * 43200), "hours": int(sel.sum()),
                             "M2": out["amplitude"]["M2"], "S2": out["amplitude"]["S2"],
                             "K1": out["amplitude"]["K1"], "noise": out["noise"]})  # fmt: skip
            except ValueError:
                pass
        start += STEP_D * 86400
    return sid, rows


def iso(t):
    return datetime.fromtimestamp(t, tz=timezone.utc).strftime("%Y-%m-%d")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    stations = list(csv.DictReader(open(ARCHIVE / "13_audit/anomaly/stations.csv")))
    phase = {
        r["sid"]: r
        for r in csv.DictReader(open(ARCHIVE / "13_audit/phase/phase_metadata.csv"))
    }
    tide = {
        r["sid"]: r
        for r in csv.DictReader(open(ARCHIVE / "13_audit/dpg_tide/dpg_tide.csv"))
    }
    todo = []
    for r in stations:
        sid = r["sid"]
        snr = float(tide.get(sid, {}).get("snr") or 0)
        if (int(r["n_steps"] or 0) > 0 and snr >= 15) or sid.startswith("YV."):
            if sid in phase:
                todo.append((sid, phase[sid]["p_channel"]))
    with ProcessPoolExecutor(int(os.environ.get("WORKERS", 6))) as ex:
        results = dict(ex.map(one, todo))
    allrows = [x for rows in results.values() for x in rows]
    with open(OUT / "monthly_tide.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(allrows[0]))
        w.writeheader()
        w.writerows(allrows)
    # M2 change across each |Z/P| step, from windows entirely before and after it
    out = []
    for r in stations:
        rows = results.get(r["sid"])
        if not rows:
            continue
        t = np.array([x["t_mid"] for x in rows])
        m2 = np.array([x["M2"] for x in rows])
        scatter = (
            100 * np.std(np.diff(np.log(m2))) / np.sqrt(2) if len(m2) > 3 else np.nan
        )
        for s in json.loads(r["steps"] or "[]"):
            ts = (
                datetime.fromisoformat(s["time"])
                .replace(tzinfo=timezone.utc)
                .timestamp()
            )
            half = WIN_D * 43200
            before = m2[t + half <= ts]
            after = m2[t - half >= ts]
            if len(before) < 1 or len(after) < 1:
                continue
            change = 100 * (np.median(after) / np.median(before) - 1)
            out.append({"sid": r["sid"], "step_time": s["time"], "tf_shift_pct": s["tf_shift_pct"],
                        "m2_change_pct": change, "m2_window_scatter_pct": scatter,
                        "n_before": len(before), "n_after": len(after),
                        "stage_d_gain_ratio": s.get("stage_d_gain_ratio"), "level_sensor": s["sensor"],
                        "verdict": ("P" if abs(change) > max(3 * scatter, 2) and np.sign(change) == -np.sign(s["tf_shift_pct"])
                                    else "not P" if abs(change) < max(0.3 * abs(s["tf_shift_pct"]), scatter) else "unclear")})  # fmt: skip
    with open(OUT / "tide_at_steps.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"{len(todo)} stations, {len(allrows)} windows, {len(out)} steps")


if __name__ == "__main__":
    main()
