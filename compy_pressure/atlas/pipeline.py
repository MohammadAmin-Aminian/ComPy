"""The atlas from a channel table and the native-rate spectra: selection, screens, models and figures.

Input channel table (CSV), one row per channel. Required columns:
    seed_id, net, sta, cha, lat, lon, depth_m, fs, type, brand, product (path of the native-rate file)
(columns named key, in_atlas or drop_reason are dropped: the pipeline writes them)
Optional columns (a missing column disables its step):
    gauge          gauge identifier (one channel per gauge is kept: ?DH first, then the widest band)
    deployment     "seafloor" for gauges on the seafloor (others are removed)
    pa_units       True when the StationXML input unit is a pressure unit
    exclude_reason non-empty: removed upstream with this reason (metadata, instrument tests)
    gain_db        gain correction applied to the levels [dB] (default 0)
    f_lo_hz, f_hi_hz  usable band from the response or validated cut-offs [Hz]
    pz_level_db    P-Z calibration level at 0.2 Hz (``screens.pz_level``), deep gauges
    coast_km, shelf_break_km   distances (``geography.gauge_distances``)

Outputs in ``out``: tables/ (population, cascade, window counts, screens, model index,
season pairs), models/<axis>/<group>.{csv,json,_density.npz}, cubes/ (memory-mapped
histograms), figures/.
"""

from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from . import classify, cube, figures, native, screens
from .config import db_edges
from .models import GRID, Atlas, band_taper, slug

LAT_BINS, LAT_LABELS = (
    [0, 23.5, 40, 55, 91],
    ["0-23.5 (tropics)", "23.5-40", "40-55", ">55 (subpolar)"],
)


def _products(args):
    i, path, cfg, drop = args
    try:
        counts, c_all, c_kept, daily = cube.channel_products(path, cfg, drop)
        return i, dict(counts, status="ok"), c_all, c_kept, daily
    except Exception as e:  # recorded, the channel is removed by the cascade
        return i, dict(status=f"failed: {type(e).__name__}: {e}"), None, None, None


def channel_cubes(table, cfg, out, workers=8, drops=None, tag="kept"):
    """Histograms of every channel (memory-mapped) and daily band levels; returns window counts and daily levels."""
    edges = db_edges(cfg)
    (out / "cubes").mkdir(parents=True, exist_ok=True)
    shape = (len(table), 12, len(GRID), len(edges) - 1)
    c_all = (
        np.lib.format.open_memmap(
            out / "cubes" / "all.npy", mode="w+", dtype=np.float32, shape=shape
        )
        if tag == "kept"
        else None
    )
    c_kept = np.lib.format.open_memmap(
        out / "cubes" / f"{tag}.npy", mode="w+", dtype=np.float32, shape=shape
    )
    rows, daily = [None] * len(table), {}
    tasks = [(i, p, cfg, (drops or {}).get(k)) for i, (p, k) in enumerate(zip(table["product"], table.seed_id))
             if isinstance(p, str) and Path(p).exists()]  # fmt: skip
    with ProcessPoolExecutor(workers) as ex:
        for i, cnt, a, k, dl in ex.map(_products, tasks, chunksize=2):
            rows[i] = dict(seed_id=table.seed_id.iloc[i], **cnt)
            if a is not None:
                if c_all is not None:
                    c_all[i] = a
                c_kept[i] = k
                daily[table.seed_id.iloc[i]] = dl
    for i in range(len(rows)):
        rows[i] = rows[i] or dict(seed_id=table.seed_id.iloc[i], status="no product")
    for c in (c_all, c_kept):
        if c is not None:
            c.flush()
    return pd.DataFrame(rows), daily


def _median_level(cube_row, band, edges):
    """Median window level [dB] of one channel over a period band, from its histogram."""
    sel = (GRID >= band[0]) & (GRID <= band[1])
    h = np.asarray(cube_row).sum(axis=0)[sel]
    tot = h.sum(axis=1)
    if not (tot > 0).any():
        return np.nan
    centres = 0.5 * (edges[:-1] + edges[1:])
    c = h.cumsum(axis=1)
    return float(
        np.median(
            [centres[np.searchsorted(x, 0.5 * n)] for x, n in zip(c, tot) if n > 0]
        )
    )


def usable_limits(t, cfg):
    """Usable frequency band of each channel [Hz]: the table's limits, the Nyquist fraction, and the 1-Hz rule."""
    p = cfg["population"]
    f_lo = t.get("f_lo_hz", pd.Series(1.0 / GRID.max(), index=t.index)).fillna(
        1.0 / GRID.max()
    )
    f_hi = t.get("f_hi_hz", pd.Series(np.inf, index=t.index)).fillna(np.inf)
    f_hi = np.minimum(f_hi, p["nyquist_fraction"] * t.fs)
    f_hi = np.where(t.fs <= 1.0, np.minimum(f_hi, 1.0 / p["one_hz_min_period_s"]), f_hi)
    return f_lo.values, np.asarray(f_hi, float)


def select(t, counts, cfg, edges, cube_kept):
    """Channel cascade; returns the table with ``drop_reason``, ``in_atlas`` and gate offsets, and the cascade table."""
    p = cfg["population"]
    c = Cascade(t)
    d = c.d
    d["n_windows"] = (
        counts.set_index("seed_id").n_kept.reindex(d.seed_id).fillna(0).values
    )
    c.step(counts.status.eq("ok").values, "no native-rate product")
    if "deployment" in d:
        c.step(d.deployment.eq("seafloor"), "moored or unverified deployment")
    c.step(d.type.ne("Other"), "logger internal or unknown channel")
    if "pa_units" in d:
        c.step(d.pa_units.astype(bool), "StationXML input units not pressure")
    if "exclude_reason" in d:
        for why in [w for w in d.exclude_reason.dropna().unique() if str(w).strip()]:
            c.step(d.exclude_reason.fillna("").ne(why), str(why))
    c.step(
        d.n_windows >= p["min_windows"],
        f"fewer than {p['min_windows']} windows (7 days)",
    )
    if "gauge" in d:
        d["band_oct"] = np.log2(d.f_hi / d.f_lo)
        d["is_dh"] = d.cha.str[1:].eq("DH")
        live = d[c.live].sort_values(
            ["gauge", "is_dh", "band_oct"], ascending=[True, False, False]
        )
        c.step(
            d.seed_id.isin(set(live.groupby("gauge").head(1).seed_id)),
            "second channel of the same gauge",
        )
    for a, b in p["level_gate_bands_s"]:
        lev = np.array(
            [
                _median_level(cube_kept[i], (a, b), edges) if c.live.iloc[i] else np.nan
                for i in range(len(d))
            ]
        )
        col = f"level_{a:g}_{b:g}s_db"
        d[col] = lev + d.gain_db.values
        off, ok = screens.level_gate(d[col], d.basin, c.live, p["level_gate_db"])
        d[col.replace("_db", "_offset_db")] = off
        c.step(
            ok,
            f"{a:g}-{b:g} s level more than {p['level_gate_db']:g} dB from basin median",
        )
    if "pz_level_db" in d:
        deep = d.depth_m >= cfg["screens"]["deep_m"]
        off = screens.pz_offsets(d.loc[c.live & deep].set_index("seed_id").pz_level_db,
                                 d.loc[c.live & deep].rename(columns={"seed_id": "key"}))  # fmt: skip
        cal = screens.calibration_screen(d.rename(columns={"seed_id": "key"}), off, cfg)
        for col in cal:
            d[col] = cal[col].reindex(d.seed_id).values
        c.step(
            ~d.cal_excluded.fillna(False).astype(bool),
            f"P-Z calibration level more than {cfg['screens']['cal_db']:g} dB from same-depth DPG/APG",
        )
    d["in_atlas"] = c.live
    return d, c.table()


Cascade = screens.Cascade


def standard_models(A):
    """The groups of the atlas: global, sensor, depth, sensor x depth, basin, distances, latitude, brand, seasons."""
    d = A.d
    allu = list(range(len(d)))
    deep = d.depth_m >= 1000
    A.add("all", "all", allu)
    A.add("all", "deep_ge200m", list(np.where(d.depth_m >= 200)[0]))
    A.by("sensor_type", "type")
    A.by("depth", "depth_bin")
    for (x, y), s in d.groupby(["type", "depth_bin"]):
        A.add("type_x_depth", f"{x}_{y}", list(s.index))
    d["lat_bin"] = classify.binned(d.lat.abs(), LAT_BINS, LAT_LABELS)
    for col, axis in (("brand", "sensor_brand"), ("basin", "basin"), ("distance_bin", "coast_distance"),
                      ("shelf_bin", "shelf_break_distance"), ("lat_bin", "latitude")):  # fmt: skip
        if col not in d:
            continue
        for g, s in d.groupby(col):
            A.add(axis, g, list(s.index))
            A.add(axis + "_deep", g, list(np.where(d.index.isin(s.index) & deep)[0]))
    for s_ in classify.SEASONS:
        A.add("season", s_, allu, lambda u, s_=s_: classify.season_months(s_, d.lat[u]))
    for hemi, m in (("north", d.lat >= 0), ("south", d.lat < 0)):
        for s_ in classify.SEASONS:
            A.add(
                "hemisphere_season",
                f"{hemi}_{s_}",
                list(np.where(m)[0]),
                lambda u, s_=s_: classify.season_months(s_, d.lat[u]),
            )
        for k in range(12):
            A.add(
                "month_" + hemi,
                f"{k + 1:02d}",
                list(np.where(m)[0]),
                lambda u, k=k: np.array([k]),
            )


def build(table_path, out, cfg, workers=8, make_figures=True):
    """Run the atlas: cubes, selection, screens, models, figures. Returns the atlas object."""
    out = Path(out)
    tdir = out / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(table_path, low_memory=False)
    t = t.drop(
        columns=[c for c in ("key", "in_atlas", "drop_reason") if c in t]
    )  # the pipeline's key is seed_id
    t["gain_db"] = t.get("gain_db", pd.Series(0.0, index=t.index)).fillna(0.0)
    t = classify.add_classes(t, cfg)
    t["f_lo"], t["f_hi"] = usable_limits(t, cfg)
    edges = db_edges(cfg)
    # 1. window screens and histograms of every channel
    counts, daily = channel_cubes(t, cfg, out, workers)
    counts.to_csv(tdir / "channel_windows.csv", index=False)
    kept = np.load(out / "cubes" / "kept.npy", mmap_mode="r")
    # 2. channel cascade
    d, casc = select(t, counts, cfg, edges, kept)
    # 3. day screen, then the histograms again without the dropped days
    a = d[d.in_atlas].rename(columns={"seed_id": "key"})
    drops, dsum = screens.day_screen({k: daily.get(k) for k in a.key}, a, cfg)
    pd.DataFrame(
        [(k, x) for k, s in drops.items() for x in sorted(s)], columns=["key", "day"]
    ).to_csv(tdir / "screen_days.csv", index=False)
    dsum.to_csv(tdir / "screen_days_summary.csv", index=False)
    a_tab = d[d.in_atlas].reset_index(drop=True)
    cnt2, daily2 = channel_cubes(a_tab, cfg, out, workers, drops, tag="screened")
    cnt2.to_csv(tdir / "channel_windows_screened.csv", index=False)
    # 4. band screen and the usable-band taper
    a_tab = a_tab.rename(columns={"seed_id": "key"})
    bands = screens.band_screen(daily2, drops, a_tab, cfg)
    bands.to_csv(tdir / "screen_gauge_band.csv", index=False)
    U = np.vstack([band_taper(lo, hi) for lo, hi in zip(a_tab.f_lo, a_tab.f_hi)])
    pos = {k: i for i, k in enumerate(a_tab.key)}
    for r in bands[bands.removed].itertuples():
        lo, hi = cfg["bands"][r.band]
        U[pos[r.key], (GRID >= lo) & (GRID <= hi)] = 0
    casc = pd.concat([casc, pd.DataFrame([
        dict(step="windows: no level, clipped, dead, deployment edges", removed=0, remaining=len(a_tab),
             windows_removed=int(counts.set_index("seed_id").reindex(a_tab.key)[["n_no_level", "n_clipped", "n_dead", "n_edge"]].sum().sum())),
        dict(step=f"day screen: day level more than {cfg['screens']['day_db']:g} dB from network peers", removed=0, remaining=len(a_tab),
             windows_removed=int(cnt2.n_day_screen.sum())),
        dict(step=f"band screen: gauge band more than {cfg['screens']['band_db']:g} dB from peers (bands removed)",
             removed=0, remaining=len(a_tab), windows_removed=int(bands.removed.sum()))])], ignore_index=True)  # fmt: skip
    casc.to_csv(tdir / "selection_cascade.csv", index=False)
    d.to_csv(tdir / "population.csv", index=False)
    np.save(out / "cubes" / "taper.npy", U)
    np.save(out / "cubes" / "periods.npy", GRID)
    np.save(out / "cubes" / "db_edges.npy", edges)
    # 5. models
    scr = np.load(out / "cubes" / "screened.npy", mmap_mode="r")
    A = Atlas(a_tab, scr, U, cfg, out / "models")
    standard_models(A)
    A.index_table().to_csv(tdir / "model_index.csv", index=False)
    A.season_paired().to_csv(tdir / "season_paired.csv", index=False)
    (tdir / "settings.json").write_text(
        json.dumps(dict(cfg, product=native.PRODUCT_VERSION), indent=1)
    )
    if make_figures:
        standard_figures(A, out / "figures")
    return A


def standard_figures(A, fig_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab = pd.read_csv(A.out / "all" / "all.csv")
    P, used = A.densities[("all", "all")]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    figures.density_panel(ax, P, A.edges, f"All gauges ({int(used.max())})",
                          {"5 %": (tab.p5_db, dict(color="w", lw=1.2)), "median": (tab.p50_db, dict(color="w", lw=2)),
                           "95 %": (tab.p95_db, dict(color="w", lw=1.2))})  # fmt: skip
    ax.set_ylim(-60, 90)
    figures.save(fig, fig_dir / "global_density.png")
    for axis, col in (("depth", "depth_bin"), ("sensor_type", "type"), ("coast_distance", "distance_bin"),
                      ("shelf_break_distance", "shelf_bin")):  # fmt: skip
        if col not in A.d:
            continue
        groups = sorted(A.d[col].dropna().unique())
        tabs = {}
        for g in groups:
            f = A.out / axis / f"{slug(g)}.csv"
            if f.exists():
                tabs[g] = pd.read_csv(f)
        cmap = plt.get_cmap("viridis", max(len(tabs), 2))
        colours = (
            figures.TYPE_COLOURS
            if axis == "sensor_type"
            else {g: cmap(i) for i, g in enumerate(tabs)}
        )
        fig, ax = plt.subplots(figsize=(9, 5.5))
        figures.model_lines(ax, tabs, colours)
        ax.legend(fontsize=8, loc="upper left", title=axis.replace("_", " "))
        ax.set_ylim(-60, 90)
        figures.save(fig, fig_dir / f"models_{axis}.png")


def main_env_threads():
    for v in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "MKL_NUM_THREADS",
    ):
        os.environ.setdefault(v, "1")
