"""Station-balanced pressure noise percentiles by group, from a scored table.

Reads a channel table (one row per pressure epoch, e.g. the completed QC and
calibration table) and selects channels by one column value, by default the
provisional final class "good". One channel is kept per station and sensor
type, the one with the widest usable band. Only PSD bins lying entirely inside
that channel's usable response band are used, and windows far below the
channel's own median (dead or clipped data) are dropped. Each station's windows
are pooled into one histogram, then stations are weighted equally at each
frequency (compy_pressure.plot_network_ppsd.station_distribution).

Groupings: all, basin, type, depth_bin and season. Season is local: months are
shifted by six for stations south of the equator, so "winter" is DJF in the
north and JJA in the south. Every grouping writes a percentile CSV, a
distribution NPZ and three figures: percentile bands, a heat map of the
station-balanced probability density, and the median PSD of every station
used as a line. A JSON summary lists the inputs. No gain is
applied and the inputs are not changed. This is a descriptive model: it is
provisional until the selection is approved, after which the same command can
be rerun on the approved table.
"""

import argparse
import csv
import json
from collections import defaultdict, namedtuple
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from obspy.imaging.cm import pqlx
from obspy.signal import PPSD

from compy_pressure.plot_network_ppsd import accumulate_histogram, station_distribution
from compy_pressure.pressure_test import load_references

PERCENTILES = (5, 10, 50, 90, 95)
GROUPINGS = ("all", "basin", "type", "depth_bin", "season")
SEASONS = ("winter", "spring", "summer", "autumn")
# Display only: frequencies covered by fewer stations are left out of figures.
MIN_STATIONS_PLOTTED = 3
# A window whose band median is this far below the channel median is dropped.
DEAD_WINDOW_DB = 30.0
DEPTH_ORDER = ("< 1 km", "1–2 km", "2–3 km", "3–4 km", "4–5 km", "> 5 km")
Model = namedtuple(
    "Model",
    "frequency curves coverage stations windows frequency_edges db_edges "
    "probability lines",
)


def local_season(month, latitude):
    """Meteorological season at the station: DJF winter in the north."""
    if latitude < 0:
        month = (month + 5) % 12 + 1
    return SEASONS[(month % 12) // 3]


def one_per_station(table):
    """Keep, per station and sensor type, the channel with the widest band."""
    width = np.log(table.band_max_hz / table.band_min_hz)
    ranked = table.assign(
        _width=width, _station=table.Seed_id.str.rsplit(".", n=2).str[0]
    )
    ranked = ranked.sort_values(
        ["_station", "type", "_width"], ascending=[True, True, False]
    )
    keep = ranked.drop_duplicates(["_station", "type"])
    return keep.drop(columns=["_width"]).rename(columns={"_station": "station"})


def band_selection(ppsd, low, high):
    """Increasing frequencies, display edges and columns of bins inside [low, high].

    A bin is kept only when its whole smoothing interval lies in the band; the
    histogram edges are the PPSD's non-overlapping display edges.
    """
    f = 1 / np.asarray(ppsd.period_bin_centers, dtype=float)
    order = np.argsort(f)
    left = (1 / np.asarray(ppsd.period_bin_right_edges))[order]
    right = (1 / np.asarray(ppsd.period_bin_left_edges))[order]
    indices = np.flatnonzero((left >= low) & (right <= high))
    if not len(indices):
        raise ValueError("No PSD bin lies inside the usable band")
    if np.any(np.diff(indices) != 1):
        raise ValueError("Usable band is not contiguous")
    edges = (1 / np.asarray(ppsd.period_xedges))[::-1][indices[0] : indices[-1] + 2]
    return f[order][indices], edges, order[indices]


def channel_histograms(row):
    """Histogram states of one channel: the whole record and each local season."""
    ppsd = PPSD.load_npz(row["ppsd_npz_file"])
    f, edges, columns = band_selection(ppsd, row["band_min_hz"], row["band_max_hz"])
    values = np.asarray(ppsd.psd_values, dtype=float)[:, columns]
    level = np.median(values, axis=1)
    live = level >= np.median(level) - DEAD_WINDOW_DB
    values = values[live]
    months = np.array([t.month for t in ppsd.times_processed])[live]
    seasons = np.array([local_season(m, row["Lat"]) for m in months])
    states = {"all": {}}
    medians = {"all": (f, np.median(values, axis=0))}
    accumulate_histogram(states["all"], f, values, edges)
    for season in SEASONS:
        if (seasons == season).any():
            states[season] = {}
            accumulate_histogram(states[season], f, values[seasons == season], edges)
            medians[season] = (f, np.median(values[seasons == season], axis=0))
    return {
        "Seed_id": row["Seed_id"],
        "station": row["station"],
        "states": states,
        "medians": medians,
        "windows": int(live.sum()),
        "dropped_windows": int((~live).sum()),
    }


def safe_histograms(row):
    try:
        return channel_histograms(row)
    except Exception as exc:
        return {"Seed_id": row["Seed_id"], "error": f"{type(exc).__name__}: {exc}"}


def percentile_curves(probability, db_edges, coverage):
    cdf = probability.cumsum(axis=0)
    levels = (db_edges[:-1] + db_edges[1:]) / 2
    curves = {}
    for p in PERCENTILES:
        values = levels[np.argmax(cdf >= p / 100, axis=0)].astype(float)
        values[coverage == 0] = np.nan
        curves[p] = values
    return curves


def group_order(grouping, names):
    order = {"season": SEASONS, "depth_bin": DEPTH_ORDER}.get(grouping)
    return [n for n in order if n in names] if order else sorted(names)


def build_models(results, rows, grouping):
    """{group: Model}: percentiles, distribution and station medians per group."""
    by_seed = rows.set_index("Seed_id")
    states, lines = defaultdict(list), defaultdict(list)
    for result in results:
        row = by_seed.loc[result["Seed_id"]]
        medians = result.get("medians", {})
        if grouping == "season":
            for season in SEASONS:
                if season in result["states"]:
                    states[season].append(result["states"][season])
                    if season in medians:
                        lines[season].append(medians[season])
        else:
            group = "all" if grouping == "all" else str(row[grouping])
            states[group].append(result["states"]["all"])
            if "all" in medians:
                lines[group].append(medians["all"])
    models = {}
    for group in group_order(grouping, states):
        f, edges, db, probability, coverage = station_distribution(states[group])
        models[group] = Model(
            f,
            percentile_curves(probability, db, coverage),
            coverage,
            len(states[group]),
            int(sum(s["windows"] for s in states[group])),
            edges,
            db,
            probability,
            lines[group],
        )
    return models


def plot_references(axis):
    metadata, brown = load_references()
    brown = brown[brown[:, 0] >= 0.2]  # used at and above 0.2 Hz, as in Test 1
    axis.fill_between(
        brown[:, 0],
        brown[:, 1],
        brown[:, 2],
        color="#4c9a5b",
        alpha=0.18,
        lw=0,
        label="Brown et al. 2014 (≥0.2 Hz)",
    )
    for i, point in enumerate(metadata["webb"]["points"]):
        axis.plot(
            [point["frequency_hz"]] * 2,
            [point["low_db_pa2_hz"], point["high_db_pa2_hz"]],
            color="#c0392b",
            lw=4,
            alpha=0.8,
            solid_capstyle="butt",
            label="Webb 1998 Fig. 2 (≤0.2 Hz)" if i == 0 else None,
        )


def plot_grouping(grouping, models, path, title):
    n = len(models)
    cols = min(n, 3)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(8.5, 5.6) if n == 1 else (5.2 * cols, 3.9 * rows),
        sharex=True,
        sharey=True,
        layout="constrained",
        squeeze=False,
    )
    for axis, (group, m) in zip(axes.flat, models.items(), strict=False):
        f, curves, coverage, stations = m.frequency, m.curves, m.coverage, m.stations
        shown = coverage >= MIN_STATIONS_PLOTTED
        plot_references(axis)
        c = {p: np.where(shown, v, np.nan) for p, v in curves.items()}
        axis.fill_between(
            f, c[5], c[95], color="#1f5f8b", alpha=0.18, lw=0, label="5–95 %"
        )
        axis.fill_between(
            f, c[10], c[90], color="#1f5f8b", alpha=0.30, lw=0, label="10–90 %"
        )
        axis.semilogx(f, c[50], color="#0b2f4a", lw=2, label="median")
        axis.set_title(f"{group}: {stations} stations")
        axis.grid(which="both", alpha=0.2)
    for axis in axes.flat[n:]:
        axis.set_visible(False)
    for axis in axes[-1]:
        axis.set_xlabel("Frequency (Hz)")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"dB re 1 Pa$^2$/Hz")
    axes.flat[0].set(xlim=(1e-3, 50), ylim=(-70, 90))
    axes.flat[0].legend(loc="lower left", fontsize=7.5)
    fig.suptitle(title)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_medians(grouping, models, path, title):
    fig, axis = plt.subplots(figsize=(9, 5.5), layout="constrained")
    plot_references(axis)
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(models)))
    for color, (group, m) in zip(colors, models.items(), strict=True):
        median = np.where(m.coverage >= MIN_STATIONS_PLOTTED, m.curves[50], np.nan)
        axis.semilogx(
            m.frequency, median, color=color, lw=2, label=f"{group} ({m.stations})"
        )
    axis.set(
        xlim=(1e-3, 50),
        ylim=(-60, 70),
        xlabel="Frequency (Hz)",
        ylabel=r"dB re 1 Pa$^2$/Hz",
        title=title,
    )
    axis.grid(which="both", alpha=0.2)
    axis.legend(fontsize=8, ncol=2)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def panel_grid(n):
    cols = min(n, 3)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(8.5, 5.6) if n == 1 else (5.2 * cols, 3.9 * rows),
        sharex=True,
        sharey=True,
        layout="constrained",
        squeeze=False,
    )
    for axis in axes.flat[n:]:
        axis.set_visible(False)
    for axis in axes[-1]:
        axis.set_xlabel("Frequency (Hz)")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"dB re 1 Pa$^2$/Hz")
    return fig, axes


def reference_outlines(axis):
    """Webb bars and Brown limits as lines, so they do not hide a heat map."""
    metadata, brown = load_references()
    brown = brown[brown[:, 0] >= 0.2]
    for column in (1, 2):
        axis.semilogx(
            brown[:, 0],
            brown[:, column],
            color="#2e7d32",
            lw=1.2,
            label="Brown et al. 2014 (≥0.2 Hz)" if column == 1 else None,
        )
    for i, point in enumerate(metadata["webb"]["points"]):
        axis.plot(
            [point["frequency_hz"]] * 2,
            [point["low_db_pa2_hz"], point["high_db_pa2_hz"]],
            color="#c0392b",
            lw=3,
            solid_capstyle="butt",
            label="Webb 1998 Fig. 2 (≤0.2 Hz)" if i == 0 else None,
        )


def plot_heatmaps(models, path, title):
    """Station-balanced probability per 1 dB bin, with the percentile curves."""
    fig, axes = panel_grid(len(models))
    shown = {g: m.coverage >= MIN_STATIONS_PLOTTED for g, m in models.items()}
    peak = max(
        np.percentile(100 * m.probability[:, shown[g]], 99.5)
        for g, m in models.items()
        if shown[g].any()
    )
    for axis, (group, m) in zip(axes.flat, models.items(), strict=False):
        density = np.ma.masked_less_equal(100 * m.probability, 0)
        density[:, ~shown[group]] = np.ma.masked
        mesh = axis.pcolormesh(
            m.frequency_edges,
            m.db_edges,
            density,
            cmap=pqlx,
            vmin=0,
            vmax=peak,
            shading="flat",
            rasterized=True,
        )
        curves = {p: np.where(shown[group], v, np.nan) for p, v in m.curves.items()}
        axis.semilogx(m.frequency, curves[50], color="k", lw=1.6, label="median")
        for p in (5, 95):
            axis.semilogx(
                m.frequency,
                curves[p],
                color="k",
                lw=0.9,
                ls="--",
                label="5th / 95th percentile" if p == 5 else None,
            )
        reference_outlines(axis)
        axis.set_title(f"{group}: {m.stations} stations")
        axis.grid(which="both", alpha=0.15)
    axes.flat[0].set(xlim=(1e-3, 50), ylim=(-70, 90))
    axes.flat[0].legend(loc="lower left", fontsize=7.5)
    fig.colorbar(
        mesh, ax=axes, shrink=0.8, label="Probability (% per dB, stations equal)"
    )
    fig.suptitle(title)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_lines(models, path, title):
    """Median PSD of every station used, with the group median on top."""
    fig, axes = panel_grid(len(models))
    for axis, (group, m) in zip(axes.flat, models.items(), strict=False):
        alpha = float(np.clip(20 / max(len(m.lines), 1), 0.18, 0.7))
        for i, (f, median) in enumerate(m.lines):
            axis.semilogx(
                f,
                median,
                color="#1f5f8b",
                lw=0.6,
                alpha=alpha,
                label="one station (median PSD)" if i == 0 else None,
            )
        shown = m.coverage >= MIN_STATIONS_PLOTTED
        axis.semilogx(
            m.frequency,
            np.where(shown, m.curves[50], np.nan),
            color="k",
            lw=2,
            label="group median",
        )
        reference_outlines(axis)
        axis.set_title(f"{group}: {m.stations} stations")
        axis.grid(which="both", alpha=0.2)
    axes.flat[0].set(xlim=(1e-3, 50), ylim=(-70, 90))
    axes.flat[0].legend(loc="lower left", fontsize=7.5)
    fig.suptitle(title)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_curves(path, models):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "group",
                "frequency_hz",
                *[f"p{p}_db_pa2_hz" for p in PERCENTILES],
                "stations",
            ]
        )
        for group, m in models.items():
            f, curves, coverage = m.frequency, m.curves, m.coverage
            for i, hz in enumerate(f):
                if coverage[i]:
                    writer.writerow(
                        [
                            group,
                            f"{hz:.6g}",
                            *[f"{curves[p][i]:.1f}" for p in PERCENTILES],
                            int(coverage[i]),
                        ]
                    )


def run(
    table_path, output, *, select_column, select_value, band_columns, workers, label
):
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"{output} is not empty; choose a new folder")
    output.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(table_path, low_memory=False)
    chosen = table[table[select_column].astype(str).eq(select_value)].rename(
        columns={band_columns[0]: "band_min_hz", band_columns[1]: "band_max_hz"}
    )
    chosen = chosen.dropna(subset=["band_min_hz", "band_max_hz", "ppsd_npz_file"])
    rows = one_per_station(chosen)
    records = rows[
        ["Seed_id", "station", "ppsd_npz_file", "band_min_hz", "band_max_hz", "Lat"]
    ].to_dict("records")
    with ProcessPoolExecutor(workers, mp_context=get_context("spawn")) as pool:
        results = list(pool.map(safe_histograms, records, chunksize=4))
    failed = [r for r in results if "error" in r]
    results = [r for r in results if "error" not in r]
    summary = {
        "table": str(Path(table_path).resolve()),
        "selection": {select_column: select_value},
        "status": label,
        "channels_selected": int(len(chosen)),
        "stations_used": len(results),
        "windows_used": int(sum(r["windows"] for r in results)),
        "windows_dropped_dead": int(sum(r["dropped_windows"] for r in results)),
        "failed": failed,
        "weighting": "equal stations per frequency; windows pooled within station",
        "dead_window_rule_db": DEAD_WINDOW_DB,
        "percentiles": list(PERCENTILES),
        "gain_applied": False,
        "groupings": {},
    }
    pd.DataFrame(
        [
            {k: r[k] for k in ("Seed_id", "station", "windows", "dropped_windows")}
            for r in results
        ]
    ).to_csv(output / "channels_used.csv", index=False)
    labels = {
        "all": "All stations",
        "basin": "Ocean basin",
        "type": "Sensor type",
        "depth_bin": "Water depth",
        "season": "Local season",
    }
    for grouping in GROUPINGS:
        models = build_models(results, rows, grouping)
        scope = (
            "of all stations" if grouping == "all" else f"by {labels[grouping].lower()}"
        )
        write_curves(output / f"noise_model_{grouping}.csv", models)
        np.savez_compressed(
            output / f"noise_model_{grouping}.npz",
            **{
                f"{g}__{key}": value
                for g, m in models.items()
                for key, value in (
                    ("frequency_hz", m.frequency),
                    ("frequency_edges_hz", m.frequency_edges),
                    ("db_edges", m.db_edges),
                    ("station_balanced_probability", m.probability),
                    ("coverage", m.coverage),
                    *((f"p{p}", m.curves[p]) for p in PERCENTILES),
                )
            },
        )
        title = f"Pressure noise model {scope} ({summary['status']})"
        plot_grouping(grouping, models, output / f"noise_model_{grouping}.png", title)
        plot_heatmaps(
            models,
            output / f"noise_model_{grouping}_heatmap.png",
            f"Station-balanced pressure PSD density {scope} ({summary['status']})",
        )
        plot_lines(
            models,
            output / f"noise_model_{grouping}_lines.png",
            f"Median PSD of every station used, {scope} ({summary['status']})",
        )
        if len(models) > 1:
            plot_medians(
                grouping,
                models,
                output / f"noise_model_{grouping}_medians.png",
                f"Median pressure noise {scope}",
            )
        summary["groupings"][grouping] = {
            g: {"stations": m.stations, "windows": m.windows} for g, m in models.items()
        }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=1, default=str) + "\n"
    )
    print(json.dumps({k: v for k, v in summary.items() if k != "failed"}, indent=1))
    return summary


def main():
    from compy_pressure.network import worker_count

    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--select-column", default="final_class")
    parser.add_argument("--select-value", default="good")
    parser.add_argument(
        "--band-columns",
        nargs=2,
        default=("pressure_expected_band_min_hz", "pressure_expected_band_max_hz"),
    )
    parser.add_argument("--workers", type=worker_count, default=4)
    parser.add_argument(
        "--label",
        default="provisional",
        help="Status written on the outputs; use 'approved' only for an approved table",
    )
    args = parser.parse_args()
    run(
        args.table,
        args.output_dir,
        select_column=args.select_column,
        select_value=args.select_value,
        band_columns=args.band_columns,
        workers=args.workers,
        label=args.label,
    )


if __name__ == "__main__":
    main()
