"""Which sensor changed? Z and P levels of each station relative to its network.

A step in |Z/P| can come from the pressure gauge or the seismometer. The
whole-archive PPSD files (``04_ppsd/<NET>/<STA>/<NET>.<STA>_cleaned.npz``) hold
per-window (3 h) PSDs of Z and P. In the microseism band (0.1-0.3 Hz) both follow
the regional sea state, so the station's level minus the median of the other
stations of the network (same 6-h bin) removes the ocean and leaves the
instrument. A sensitivity change of one sensor gives a step in that channel's
relative level only.

For each station: changepoints (``compy.anomaly.changepoints``) of the relative
Z and P levels in the microseism and the compliance band, and the step at given
dates (median of the 20 days after minus the 20 days before, with a robust
z-score against the same statistic at all other dates).

Output: ``<OUT>/level_steps/<NET>_relative_levels.npz``, ``level_steps.csv``
(changepoints), ``dated_steps.csv``.
"""

import csv
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from compy import anomaly  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
OUT = Path(os.environ.get("COMPY_AUDIT_OUT", ARCHIVE / "13_audit")) / "level_steps"
BANDS = {"ms": (0.1, 0.3), "ig": (0.007, 0.015)}
BIN_S = 6 * 3600
HALF_DAYS = 20
#: Steps reported by the ComPy 3 audit (PR #21) to test here.
DATED = {
    "YV": {"2013-09-14": ("RR28", "RR29", "RR40", "RR52"), "2013-06-06": ("RR50",),
           "2013-07-30": ("RR50",), "2013-08-05": ("RR36",), "2013-10-21": ("RR36",)},
}  # fmt: skip


def band_levels(path):
    d = np.load(path, allow_pickle=True)
    f = d["freqs"]
    out = {}
    for ch in ("Z", "P"):
        if f"{ch}_psd_db" not in d.files:
            return None
        t = d[f"{ch}_times"]
        db = d[f"{ch}_psd_db"].astype(float)
        for name, (lo, hi) in BANDS.items():
            sel = (f >= lo) & (f <= hi)
            # mean in power, back to dB
            out[(ch, name)] = (
                t,
                10 * np.log10(np.nanmean(10 ** (db[:, sel] / 10), axis=1)),
            )
    return out


def network(net, stations):
    levels = {}
    for sta in stations:
        path = ARCHIVE / "04_ppsd" / net / sta / f"{net}.{sta}_cleaned.npz"
        if path.exists():
            lv = band_levels(path)
            if lv:
                levels[sta] = lv
    if len(levels) < 4:
        return None
    t0 = min(v[("Z", "ms")][0].min() for v in levels.values())
    t1 = max(v[("Z", "ms")][0].max() for v in levels.values())
    edges = np.arange(t0, t1 + BIN_S, BIN_S)
    centres = edges[:-1] + BIN_S / 2
    grid = {}
    for key in [(c, b) for c in ("Z", "P") for b in BANDS]:
        m = np.full((len(levels), len(centres)), np.nan)
        for i, lv in enumerate(levels.values()):
            t, v = lv[key]
            idx = np.searchsorted(edges, t, side="right") - 1
            ok = (idx >= 0) & (idx < len(centres)) & np.isfinite(v)
            for j in np.unique(idx[ok]):
                m[i, j] = np.median(v[ok][idx[ok] == j])
        grid[key] = m
    rel = {}
    names = list(levels)
    for key, m in grid.items():
        r = np.full_like(m, np.nan)
        for i in range(len(names)):
            others = np.delete(m, i, axis=0)
            n_ok = np.sum(np.isfinite(others), axis=0)
            with np.errstate(all="ignore"):
                med = np.nanmedian(others, axis=0)
            r[i] = np.where(n_ok >= 3, m[i] - med, np.nan)
        rel[key] = r
    return names, centres, rel


def dated_step(t, x, when, half=HALF_DAYS):
    ok = np.isfinite(x)
    t, x = t[ok], x[ok]

    def stat(c):
        a = x[(t >= c - half * 86400) & (t < c)]
        b = x[(t >= c) & (t < c + half * 86400)]
        if len(a) < 20 or len(b) < 20:
            return np.nan
        return np.median(b) - np.median(a)

    s = stat(when)
    null = np.array([stat(c) for c in t[:: max(1, len(t) // 400)]])
    null = null[np.isfinite(null)]
    z = (
        (s - np.median(null)) / max(anomaly.mad_sigma(null), 1e-6)
        if len(null) > 20
        else np.nan
    )
    return s, z


def main():
    from datetime import datetime, timezone

    OUT.mkdir(parents=True, exist_ok=True)
    plan = list(
        csv.DictReader(
            open(ARCHIVE / "12_inversion_all_stations/tables/station_plan.csv")
        )
    )
    nets = sorted({r["network"] for r in plan})
    only = set(sys.argv[1:])
    cp_rows, dated_rows = [], []
    for net in nets:
        if only and net not in only:
            continue
        # every station of the network with a PPSD (not only inverted ones) as reference
        stations = sorted(
            p.name for p in (ARCHIVE / "04_ppsd" / net).iterdir() if p.is_dir()
        )
        res = network(net, stations)
        if res is None:
            continue
        names, centres, rel = res
        np.savez(OUT / f"{net}_relative_levels.npz", stations=np.array(names), times=centres,
                 **{f"{c}_{b}": v for (c, b), v in rel.items()})  # fmt: skip
        inverted = {r["station"] for r in plan if r["network"] == net}
        for i, sta in enumerate(names):
            if sta not in inverted:
                continue
            for (ch, band), m in rel.items():
                x = m[i]
                ok = np.isfinite(x)
                if ok.sum() < 120:
                    continue
                for c in anomaly.changepoints(
                    centres[ok], x[ok], min_size=40, threshold=6.0
                ):
                    if abs(c.shift) < 1.0:
                        continue
                    cp_rows.append({"network": net, "station": sta, "channel": ch, "band": band,
                                    "time": anomaly._iso(c.time), "shift_db": round(c.shift, 2),
                                    "sigma": round(c.significance, 1), "n_before": c.n_before,
                                    "n_after": c.n_after})  # fmt: skip
            for day, stas in DATED.get(net, {}).items():
                when = (
                    datetime.fromisoformat(day).replace(tzinfo=timezone.utc).timestamp()
                )
                row = {
                    "network": net,
                    "station": sta,
                    "date": day,
                    "listed": sta in stas,
                }
                for (ch, band), m in rel.items():
                    s, z = dated_step(centres, m[i], when)
                    row[f"{ch}_{band}_db"] = (
                        round(float(s), 2) if np.isfinite(s) else ""
                    )
                    row[f"{ch}_{band}_z"] = round(float(z), 1) if np.isfinite(z) else ""
                dated_rows.append(row)
        print(net, len(names), "stations", flush=True)
    for name, rows in (("level_steps.csv", cp_rows), ("dated_steps.csv", dated_rows)):
        if rows:
            with open(OUT / name, "w", newline="") as fh:
                w = csv.DictWriter(fh, list(rows[0]))
                w.writeheader()
                w.writerows(rows)


if __name__ == "__main__":
    main()
