"""Search for unusual coherence in every station's 3-h windows (04_ppsd, about 800 stations).

Per station, from the per-window cross-spectra of the cleaned data
(``04_ppsd/<NET>/<STA>/<NET>.<STA>_cleaned.npz``):

* band coherence gamma^2 of Z-P in the infragravity band (5 mHz to 0.8 f_c),
  in the "gap" above f_c (1.5 f_c to 60 mHz, where seafloor pressure from
  surface waves is attenuated and Z-P coherence is not expected), below the
  gauge corner (1-2.5 mHz) and in the microseism band (0.1-0.3 Hz);
* Z-H coherence (the larger of Z-1, Z-2) at 5-35 mHz (tilt);
* Z/P phase in the infragravity band;
* narrow lines in the median Z-P coherence spectrum.

Network level: each station's per-window series becomes a robust z-score
(its own median and MAD); in 6-h bins the number of stations with
|z| > 3 is compared with a null made by circularly shifting every station's
series (500 permutations). Bins above the null's 99.9th percentile of the
maximum are "synchronous episodes" and are matched with the USGS catalogue
(M >= 6.5, from 1 h before to 12 h after) and with the network's pressure
level (storms).

Output in ``13_audit/research/coherence``: ``stations.csv``, ``lines.csv``,
``episodes.csv``, ``series/<NET>.npz``.
"""

import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from compy import anomaly  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
OUT = ARCHIVE / "13_audit/research/coherence"
BIN_S = 6 * 3600
G = 9.81


def band(f, lo, hi):
    return (f >= lo) & (f <= hi)


def gamma2(szo, szz, soo, sel):
    with np.errstate(invalid="ignore", divide="ignore"):
        g = np.abs(szo[:, sel]) ** 2 / (szz[:, sel] * soo[:, sel])
    return np.nanmedian(g, axis=1)


def one(task):
    net, sta, depth = task
    path = ARCHIVE / "04_ppsd" / net / sta / f"{net}.{sta}_cleaned.npz"
    try:
        d = np.load(path, allow_pickle=True)
    except Exception:
        return None
    if "ZP_szo" not in d.files:
        return None
    f = d["freqs"]
    fc = np.sqrt(G / (2 * np.pi * depth)) if depth > 0 else np.nan
    t = d["ZP_times"]
    szo, szz, soo = (
        d["ZP_szo"].astype(complex),
        d["ZP_szz"].astype(float),
        d["ZP_soo"].astype(float),
    )
    out = {
        "sid": f"{net}.{sta}",
        "network": net,
        "depth_m": depth,
        "n_windows": len(t),
        "fc_hz": fc,
    }
    bands = {
        "ig": (0.005, 0.8 * fc),
        "gap": (1.5 * fc, 0.06),
        "low": (0.001, 0.0025),
        "ms": (0.1, 0.3),
    }
    series = {"t": t}
    for name, (lo, hi) in bands.items():
        sel = band(f, lo, hi)
        if not np.isfinite(lo) or lo >= hi or sel.sum() < 1:
            continue
        g = gamma2(szo, szz, soo, sel)
        series[f"zp_{name}"] = g
        out[f"zp_{name}_median"] = float(np.nanmedian(g))
        out[f"zp_{name}_p90"] = float(np.nanpercentile(g, 90))
        stacked = np.abs(szo[:, sel].sum(0)) ** 2 / (
            szz[:, sel].sum(0) * soo[:, sel].sum(0)
        )
        out[f"zp_{name}_stacked"] = float(np.nanmedian(stacked))
    sel = band(f, 0.005, min(0.035, 0.8 * fc if np.isfinite(fc) else 0.035))
    if sel.any():
        with np.errstate(invalid="ignore"):
            ph = -np.degrees(
                np.angle(szo[:, sel].sum(1))
            )  # Z/P phase = angle(conj(S_zp))
        series["zp_phase"] = ph
    zh = []
    for key in ("Z1", "Z2"):
        if f"{key}_szo" in d.files and len(d[f"{key}_times"]) == len(t):
            s = band(f, 0.005, 0.035)
            zh.append(gamma2(d[f"{key}_szo"].astype(complex), d[f"{key}_szz"].astype(float),
                             d[f"{key}_soo"].astype(float), s))  # fmt: skip
    if zh:
        series["zh"] = np.nanmax(np.vstack(zh), axis=0)
        out["zh_median"] = float(np.nanmedian(series["zh"]))
    if "P_psd_db" in d.files and len(d["P_times"]) == len(t):
        series["p_ig_db"] = np.nanmean(
            d["P_psd_db"][:, band(f, 0.005, 0.02)].astype(float), axis=1
        )
    # narrow lines in the median coherence spectrum
    with np.errstate(invalid="ignore", divide="ignore"):
        gs = np.nanmedian(np.abs(szo) ** 2 / (szz * soo), axis=0)
    lines = []
    for i in range(3, len(f) - 3):
        nb = np.r_[gs[i - 3 : i - 1], gs[i + 2 : i + 4]]
        if (
            np.isfinite(gs[i])
            and gs[i] - np.nanmedian(nb) > 0.2
            and gs[i] == np.nanmax(gs[i - 1 : i + 2])
        ):
            lines.append({"sid": out["sid"], "freq_hz": float(f[i]), "coherence": float(gs[i]),
                          "excess": float(gs[i] - np.nanmedian(nb))})  # fmt: skip
    out["median_coherence_spectrum"] = gs
    return out, series, lines


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "series").mkdir(exist_ok=True)
    st = {
        (r["network"], r["station"]): r
        for r in csv.DictReader(open(ARCHIVE / "00_inventory/stations.csv"))
    }
    tasks = []
    for path in sorted((ARCHIVE / "04_ppsd").glob("*/*/*_cleaned.npz")):
        net, sta = path.parent.parent.name, path.parent.name
        try:
            depth = float(st.get((net, sta), {}).get("water_depth_m") or "nan")
        except ValueError:
            depth = float("nan")
        tasks.append((net, sta, depth))
    with ProcessPoolExecutor(int(os.environ.get("WORKERS", 8))) as ex:
        res = [r for r in ex.map(one, tasks, chunksize=4) if r]
    rows = [r[0] for r in res]
    freqs = np.load(next((ARCHIVE / "04_ppsd").glob("*/*/*_cleaned.npz")))["freqs"]
    np.savez(OUT / "median_coherence_spectra.npz", freqs=freqs, sid=np.array([r["sid"] for r in rows]),
             gamma2=np.array([r.pop("median_coherence_spectrum") for r in rows]))  # fmt: skip
    keys = sorted({k for r in rows for k in r})
    with open(OUT / "stations.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    lines = [x for r in res for x in r[2]]
    with open(OUT / "lines.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, ["sid", "freq_hz", "coherence", "excess"])
        w.writeheader()
        w.writerows(lines)
    # ---- network-synchronous episodes
    cat = list(csv.DictReader(open(ARCHIVE / "00_inventory/catalog_usgs_M5.85.csv")))
    quakes = [(datetime.fromisoformat(c["time"].replace("Z", "+00:00")).timestamp(), float(c["magnitude"]),
               c["latitude"], c["longitude"]) for c in cat if float(c["magnitude"]) >= 6.5]  # fmt: skip
    qt = np.array([q[0] for q in quakes])
    episodes = []
    rng = np.random.default_rng(0)
    by_net = {}
    for r in res:
        by_net.setdefault(r[0]["network"], []).append(r)
    for net, items in by_net.items():
        if len(items) < 5:
            continue
        t0 = min(np.nanmin(s["t"]) for _, s, _ in items)
        t1 = max(np.nanmax(s["t"]) for _, s, _ in items)
        edges = np.arange(t0, t1 + BIN_S, BIN_S)
        nb = len(edges) - 1
        save = {"edges": edges, "stations": np.array([o["sid"] for o, _, _ in items])}
        for metric in ("zp_ig", "zp_gap", "zh", "zp_phase"):
            Z = np.full((len(items), nb), np.nan)
            for i, (_, s, _) in enumerate(items):
                if metric not in s:
                    continue
                x = s[metric]
                ok = np.isfinite(x)
                if ok.sum() < 100:
                    continue
                if metric == "zp_phase":
                    x = (x - np.nanmedian(x[ok]) + 180) % 360 - 180
                z = (x - np.nanmedian(x[ok])) / max(anomaly.mad_sigma(x[ok]), 1e-6)
                idx = np.clip(
                    np.searchsorted(edges, s["t"], side="right") - 1, 0, nb - 1
                )
                for j in np.unique(idx[ok]):
                    Z[i, j] = np.nanmedian(z[ok][idx[ok] == j])
            save[metric] = Z
            active = np.isfinite(Z).sum(0)
            for sign, label in ((-1, "drop"), (1, "rise")):
                hits = np.nansum(sign * Z > 3, axis=0)
                frac = np.where(active >= 5, hits / np.maximum(active, 1), 0)
                null = []
                for _ in range(200):
                    Zs = np.array([np.roll(row, rng.integers(nb)) for row in Z])
                    h = np.nansum(sign * Zs > 3, axis=0)
                    null.append(
                        np.max(np.where(active >= 5, h / np.maximum(active, 1), 0))
                    )
                thr = np.percentile(null, 99.9)
                for j in np.where(frac > thr)[0]:
                    tc = 0.5 * (edges[j] + edges[j + 1])
                    near = [
                        (q[1], round((tc - q[0]) / 3600, 1))
                        for q, dt in zip(quakes, tc - qt)
                        if -6 * 3600 < dt < 18 * 3600
                    ]
                    episodes.append({"network": net, "metric": metric, "direction": label,
                                     "time": datetime.fromtimestamp(tc, tz=timezone.utc).isoformat()[:16],
                                     "fraction": round(float(frac[j]), 3), "active": int(active[j]),
                                     "null_p999": round(float(thr), 3),
                                     "quakes_m6.5_hours_before": str(near)})  # fmt: skip
        np.savez(OUT / "series" / f"{net}.npz", **save)
        print(net, len(items), "stations", flush=True)
    with open(OUT / "episodes.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(episodes[0]) if episodes else ["network"])
        w.writeheader()
        w.writerows(episodes)
    print(f"{len(rows)} stations, {len(lines)} lines, {len(episodes)} episode bins")


if __name__ == "__main__":
    main()
