"""Overlay saved network PPSDs by station, in frequency, without recalculation.

One figure per channel keeps pressure and acceleration units separate. By
default every saved window contributes to exact per-epoch 5/50/95 percentiles.
--all-curves additionally draws every window, without subsampling (can be slow).
--curves-only draws only individual windows, with no percentile computation.
No gain, raw waveform processing, master update or station selection is applied.
"""

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from obspy.imaging.cm import obspy_sequential, pqlx
from obspy.signal import PPSD
from obspy.signal.spectral_estimation import get_nhnm, get_nlnm

from compy_pressure.plot_ppsd import MASTER_CSV, PROJECT_DIRECTORY
from compy_pressure.pressure_test import REFERENCE_COLORS, load_references

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

CHANNELS = ("BDH", "BHZ", "BH1", "BH2")


def accumulate_histogram(state, frequencies, values, frequency_edges):
    """Pool native display bins on their edge union, without PSD interpolation."""
    if (
        len(frequency_edges) != len(frequencies) + 1
        or np.any(np.diff(frequency_edges) <= 0)
        or np.any(frequencies <= frequency_edges[:-1])
        or np.any(frequencies >= frequency_edges[1:])
    ):
        raise ValueError("Invalid native frequency bins")
    lo = int(np.floor(values.min()))
    hi = int(np.floor(values.max())) + 1
    if not state:
        state.update(
            edges=frequency_edges.copy(),
            lo=lo,
            hi=hi,
            counts=np.zeros((hi - lo, len(frequencies)), dtype=np.uint64),
            windows=0,
        )
    new_lo, new_hi = min(lo, state["lo"]), max(hi, state["hi"])
    old_counts = np.pad(
        state["counts"], ((state["lo"] - new_lo, new_hi - state["hi"]), (0, 0))
    )
    db_edges = np.arange(new_lo, new_hi + 1, dtype=float)
    native = np.zeros((new_hi - new_lo, len(frequencies)), dtype=np.uint64)
    for column in range(len(frequencies)):
        native[:, column] = np.histogram(values[:, column], db_edges)[0].astype(
            np.uint64
        )
    union = np.union1d(state["edges"], frequency_edges)
    centers = np.sqrt(union[:-1] * union[1:])
    combined = np.zeros((new_hi - new_lo, len(centers)), dtype=np.uint64)
    for source_edges, source in [
        (state["edges"], old_counts),
        (frequency_edges, native),
    ]:
        index = np.searchsorted(source_edges, centers, side="right") - 1
        keep = (index >= 0) & (index < source.shape[1])
        combined[:, keep] += source[:, index[keep]]
    state.update(lo=new_lo, hi=new_hi, counts=combined, edges=union)
    state["windows"] += len(values)


def station_distribution(states):
    """Normalize each station first, then weight available stations equally."""
    if not states:
        raise ValueError("No station histograms")
    edges = np.unique(np.concatenate([s["edges"] for s in states]))
    centers = np.sqrt(edges[:-1] * edges[1:])
    lo, hi = min(s["lo"] for s in states), max(s["hi"] for s in states)
    total = np.zeros((hi - lo, len(centers)))
    coverage = np.zeros(len(centers), dtype=int)
    for state in states:
        index = np.searchsorted(state["edges"], centers, side="right") - 1
        keep = (index >= 0) & (index < state["counts"].shape[1])
        bins = np.flatnonzero(keep)
        denominators = state["counts"].sum(axis=0)
        bins = bins[denominators[index[bins]] > 0]
        total[state["lo"] - lo : state["hi"] - lo, bins] += (
            state["counts"][:, index[bins]] / denominators[index[bins]]
        )
        coverage[bins] += 1
    probability = total / np.maximum(coverage, 1)
    if not np.allclose(probability.sum(axis=0)[coverage > 0], 1):
        raise RuntimeError("Station-balanced probability does not normalize")
    return centers, edges, np.arange(lo, hi + 1), probability, coverage


def sha256(path):
    # hashlib.file_digest needs Python 3.11; ComPy supports 3.10
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def saved_rows(master, network, channels):
    """Keep full epoch identities, but never count one NPZ twice."""
    with Path(master).open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"Seed_id", "ppsd_npz_file", "ppsd_status", "ppsd_units"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Master lacks saved PPSD paths/status/units")
        rows = []
        seen = {}
        for row in reader:
            parts = row["Seed_id"].split(".")
            if (
                len(parts) != 4
                or parts[0] != network
                or parts[3] not in channels
                or row["ppsd_status"] not in {"complete", "partial"}
                or not row["ppsd_npz_file"]
            ):
                continue
            path = Path(row["ppsd_npz_file"]).expanduser()
            path = (PROJECT_DIRECTORY / path).resolve()
            identity = tuple(
                row.get(k, "")
                for k in (
                    "Seed_id",
                    "channel_start_time",
                    "channel_end_time",
                    "inventory_file",
                    "ppsd_units",
                )
            )
            if path in seen:
                if seen[path] != identity:
                    raise ValueError(f"One PPSD has conflicting master epochs: {path}")
                continue
            seen[path] = identity
            rows.append(row | {"ppsd_npz_file": str(path)})
    return sorted(rows, key=lambda r: (r["Seed_id"], r.get("channel_start_time", "")))


def spectrum(ppsd, row):
    """Validate identity/units and return increasing frequencies and all windows."""
    if ppsd.id != row["Seed_id"]:
        raise ValueError("PPSD identifier does not match master row")
    pressure = ppsd.special_handling == "hydrophone"
    expected = "Pa^2/Hz" if pressure else "m^2/s^4/Hz"
    if ppsd.special_handling not in (None, "hydrophone"):
        raise ValueError("Unsupported PPSD response handling")
    if row["ppsd_units"] != expected:
        raise ValueError(f"PPSD/master unit mismatch: expected {expected}")
    f = 1 / np.asarray(ppsd.period_bin_centers, dtype=float)
    values = np.asarray(ppsd.psd_values)
    if (
        values.ndim != 2
        or values.shape != (len(ppsd.times_processed), len(f))
        or not values.size
        or not np.isfinite(values).all()
        or not np.isfinite(f).all()
        or np.any(f <= 0)
        or len(np.unique(f)) != len(f)
    ):
        raise ValueError("Invalid or empty saved PPSD frequency/window array")
    if len({str(t) for t in ppsd.times_processed}) != len(values):
        raise ValueError("Duplicate PPSD times")
    order = np.argsort(f)
    return f[order], values[:, order], pressure


def add_references(axis, pressure, lower, upper, *, include_unverified_points=None):
    """Do not extrapolate digitized pressure bounds or mix pressure/seismic units."""
    if pressure:
        metadata, brown = load_references()
        # Brown above, Webb (1998) at and below the 0.2 Hz split.
        brown_lower = max(lower, metadata["brown"].get("minimum_comparison_hz", 0))
        for col, label, style in [(1, "low", "--"), (2, "high", "-")]:
            keep = (brown[:, 0] >= brown_lower) & (brown[:, 0] <= upper)
            axis.plot(
                brown[keep, 0],
                brown[keep, col],
                color=REFERENCE_COLORS["brown_figure8"],
                ls=style,
                lw=2.2,
                label=f"Brown et al. 2014 {label} (approx. digitization)",
                zorder=6,
            )
        include_points = (
            metadata["webb"].get("use_in_plot", False)
            if include_unverified_points is None
            else include_unverified_points
        )
        points = metadata["webb"]["points"] if include_points else []
        for bound, marker in [("low", "v"), ("high", "^")]:
            pts = [p for p in points if lower <= p["frequency_hz"] <= upper]
            if not pts:
                continue
            axis.plot(
                [p["frequency_hz"] for p in pts],
                [p[f"{bound}_db_pa2_hz"] for p in pts],
                ls="none",
                marker=marker,
                color=REFERENCE_COLORS["webb_instructions"],
                ms=7,
                zorder=7,
                label=f"Webb 1998 Fig. 2 {bound}",
            )
    else:
        for getter, label, style in [
            (get_nlnm, "Peterson NLNM", "--"),
            (get_nhnm, "Peterson NHNM", "-"),
        ]:
            periods, db = getter()
            f = 1 / np.asarray(periods)
            keep = (f >= lower) & (f <= upper)
            order = np.argsort(f[keep])
            axis.plot(
                f[keep][order],
                np.asarray(db)[keep][order],
                color="black",
                ls=style,
                lw=2.2,
                label=label,
                zorder=6,
            )


def plot_network(
    master=MASTER_CSV,
    *,
    network="YV",
    channels=CHANNELS,
    output=None,
    all_curves=False,
    show=False,
    pressure_limits=(-100.0, 120.0),
    curves_only=False,
    single_color=False,
    probability=False,
    cmap="pqlx",
    max_percentage=30.0,
    min_percentage=0.0,
):
    if (
        not np.isfinite([min_percentage, max_percentage]).all()
        or not 0 <= min_percentage < max_percentage <= 100
    ):
        raise ValueError("Require 0 <= minimum occurrence < maximum <= 100 percent")
    heat_cmap = {"pqlx": pqlx, "obspy-default": obspy_sequential}.get(cmap)
    if heat_cmap is None:
        heat_cmap = plt.get_cmap(cmap)
    if probability and (all_curves or curves_only):
        raise ValueError("Choose probability heat map or individual curves, not both")
    single_color = single_color or probability
    all_curves = all_curves or curves_only
    if (
        len(pressure_limits) != 2
        or not np.isfinite(pressure_limits).all()
        or pressure_limits[0] >= pressure_limits[1]
    ):
        raise ValueError("Pressure limits require two increasing finite dB values")
    master = Path(master).expanduser().resolve()
    master_before = sha256(master)
    rows = saved_rows(master, network, channels)
    if not rows:
        raise ValueError("No registered saved PPSDs match the requested channels")
    output = (
        Path(output)
        if output
        else (
            PROJECT_DIRECTORY
            / "results/REVIEWS/PPSD_NETWORK"
            / f"run_{datetime.now(UTC):%Y%m%dT%H%M%SZ}_{uuid4().hex[:8]}"
        )
    )
    output = output.expanduser().resolve()
    raw = (Path.home() / "Desktop/Data_Bruit_FM").resolve()
    exports = (PROJECT_DIRECTORY / "exports").resolve()
    if any(output == p or p in output.parents for p in [raw, exports]):
        raise ValueError("Output must not be inside raw inputs or production exports")
    output.mkdir(parents=True, exist_ok=False)
    stations = sorted({r["Seed_id"].split(".")[1] for r in rows})
    colors = dict(
        zip(
            stations,
            plt.get_cmap("turbo")(np.linspace(0, 1, len(stations))),
            strict=True,
        )
    )
    report = {
        "network": network,
        "master": str(master),
        "mode": (
            "probability_heatmap"
            if probability
            else "all_windows_only"
            if curves_only
            else "all_windows_plus_percentiles"
            if all_curves
            else "percentiles"
        ),
        "percentiles": [] if curves_only or probability else [5, 50, 95],
        "probability_weighting": "equal saved windows, not equal stations"
        if probability
        else None,
        "single_color": single_color,
        "inputs": [],
        "figures": [],
        "gain_applied": False,
        "master_updated": False,
        "pressure_detail_limits_db": list(pressure_limits),
        "probability_display": {
            "colormap": heat_cmap.name,
            "minimum_shown_percent": min_percentage,
            "colour_scale_max_percent": max_percentage,
            "histogram_counts_modified": False,
        }
        if probability
        else None,
    }
    for channel in dict.fromkeys(channels):
        selected = [r for r in rows if r["Seed_id"].split(".")[-1] == channel]
        if not selected:
            continue
        figure, axis = plt.subplots(figsize=(16, 8.5))
        figure.subplots_adjust(left=0.075, right=0.77, bottom=0.19, top=0.88)
        present, lower, upper, pressure_kind, window_count = set(), np.inf, 0, None, 0
        histogram = {}
        try:
            for row in selected:
                path = Path(row["ppsd_npz_file"])
                record = {
                    k: row.get(k, "")
                    for k in (
                        "Seed_id",
                        "channel_start_time",
                        "channel_end_time",
                        "inventory_file",
                        "ppsd_npz_file",
                        "ppsd_status",
                        "ppsd_units",
                    )
                }
                if not path.is_file():
                    report["inputs"].append(record | {"plot_status": "missing_npz"})
                    print(f"Missing saved file: {path}", flush=True)
                    continue
                before = sha256(path)
                ppsd = PPSD.load_npz(str(path))
                f, values, pressure = spectrum(ppsd, row)
                if pressure_kind is not None and pressure != pressure_kind:
                    raise ValueError("Cannot mix physical units on one channel axis")
                pressure_kind = pressure
                station = row["Seed_id"].split(".")[1]
                color = "#4477AA" if single_color else colors[station]
                if all_curves:
                    curves = np.empty((len(values), len(f), 2), dtype=np.float32)
                    curves[:, :, 0], curves[:, :, 1] = f, values
                    axis.add_collection(
                        LineCollection(
                            curves,
                            colors=[color],
                            linewidths=0.25,
                            alpha=0.025,
                            rasterized=True,
                        )
                    )
                if probability:
                    accumulate_histogram(
                        histogram, f, values, 1 / np.asarray(ppsd.period_xedges)[::-1]
                    )
                elif not curves_only:
                    low, median, high = np.percentile(values, [5, 50, 95], axis=0)
                    axis.plot(f, median, color=color, lw=1.2)
                    axis.plot(f, low, color=color, lw=0.5, ls="--", alpha=0.5)
                    axis.plot(f, high, color=color, lw=0.5, ls="--", alpha=0.5)
                lower, upper = min(lower, f[0]), max(upper, f[-1])
                present.add(station)
                window_count += len(values)
                if sha256(path) != before:
                    raise RuntimeError(f"Saved PPSD changed during plotting: {path}")
                report["inputs"].append(
                    record
                    | {
                        "plot_status": "plotted",
                        "sha256": before,
                        "windows": len(values),
                        "first_window": str(ppsd.times_processed[0]),
                        "last_window": str(ppsd.times_processed[-1]),
                        "frequency_min_hz": float(f[0]),
                        "frequency_max_hz": float(f[-1]),
                    }
                )
                print(f"{ppsd.id}: {len(values):,} saved windows", flush=True)
            if not present:
                plt.close(figure)
                continue
            if probability:
                counts = histogram["counts"]
                coverage = counts.sum(axis=0)
                if np.any(coverage > window_count):
                    raise RuntimeError("Probability histogram double-counted windows")
                db_edges = np.arange(histogram["lo"], histogram["hi"] + 1)
                occurrence = np.ma.masked_equal(
                    100.0 * counts / np.maximum(coverage, 1), 0
                )
                occurrence = np.ma.masked_less(occurrence, min_percentage)
                mesh = axis.pcolormesh(
                    histogram["edges"],
                    db_edges,
                    occurrence,
                    cmap=heat_cmap,
                    shading="flat",
                    rasterized=True,
                    vmin=0,
                    vmax=max_percentage,
                )
                figure.colorbar(
                    mesh, ax=axis, label="Occurrence [% per 1-dB bin]", pad=0.02
                )
                np.savez_compressed(
                    output / f"{network}_{channel}_probability.npz",
                    counts=counts,
                    frequency_edges_hz=histogram["edges"],
                    db_edges=db_edges,
                    windows=window_count,
                    windows_per_frequency=coverage,
                )
            add_references(axis, pressure_kind, lower, upper)
            axis.set(
                xscale="log",
                xlim=(lower, upper),
                xlabel="Frequency [Hz]",
                ylabel=(
                    r"Pressure PSD [dB re 1 Pa$^2$/Hz]"
                    if pressure_kind
                    else r"Acceleration PSD [dB re 1 (m/s$^2$)$^2$/Hz]"
                ),
            )
            axis.grid(which="both", alpha=0.2)
            axis.autoscale_view(scalex=False)
            axis.set_title(
                f"{network} {channel}: {len(present)} stations, "
                f"{window_count:,} saved PSD windows\n"
                + (
                    f"Occurrence [%]; colour scale 0–{max_percentage:g}%; "
                    f"display floor {min_percentage:g}% (counts unchanged)"
                    if probability
                    else "Every saved window; no averaging or subsampling"
                    if curves_only
                    else "Solid: epoch median; thin dashed: 5th/95th percentiles"
                )
            )
            models = axis.legend(loc="lower left", fontsize=8, framealpha=0.9)
            if not single_color:
                axis.add_artist(models)
                axis.legend(
                    handles=[
                        Line2D([0], [0], color=colors[s], label=s)
                        for s in sorted(present)
                    ],
                    title="Station",
                    loc="upper left",
                    bbox_to_anchor=(1.01, 1),
                    ncol=2,
                    fontsize=8,
                    frameon=False,
                )
            else:
                figure.subplots_adjust(right=0.96)
            figure.text(
                0.075,
                0.055,
                "Saved nominal-response spectra only. No gain or acceptance filter. "
                "Full saved frequency range is not a verified sensor passband.\n"
                "Observations are not a new validated noise model. "
                "Pressure bounds: Webb 1998 Fig. 2 points (≤0.2 Hz) / "
                "Brown digitization (≥0.2 Hz); no extrapolation.",
                fontsize=8,
            )
            suffix = "probability_ppsd" if probability else "station_ppsd"
            path = output / f"{network}_{channel}_{suffix}.png"
            if pressure_kind:
                full_range = output / f"{network}_{channel}_{suffix}_full_range.png"
                figure.savefig(full_range, dpi=160)
                report["figures"].append(str(full_range))
                axis.set_ylim(*pressure_limits)
                axis.set_title(
                    axis.get_title() + f"\nDetail view: {pressure_limits[0]:g} to "
                    f"{pressure_limits[1]:g} dB; full-range PNG also saved"
                )
            figure.savefig(path, dpi=160)
            report["figures"].append(str(path))
        finally:
            if not show:
                plt.close(figure)
    if sha256(master) != master_before:
        raise RuntimeError("Master changed during plotting")
    for record in report["inputs"]:
        if record["plot_status"] == "plotted" and (
            sha256(record["ppsd_npz_file"]) != record["sha256"]
        ):
            raise RuntimeError("Saved PPSD changed before plotting finished")
    if not report["figures"]:
        raise ValueError("No readable saved PPSD files; no figure generated")
    report["master_sha256"] = master_before
    report["protected_inputs_unchanged"] = True
    with (output / "plot_manifest.json").open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"Plots and read-only manifest: {output}", flush=True)
    if show:
        plt.show()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-csv", type=Path, default=MASTER_CSV)
    parser.add_argument("--network", default="YV")
    parser.add_argument("--channels", nargs="+", default=list(CHANNELS))
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--all-curves",
        action="store_true",
        help="Draw every individual window as well (slow/memory-heavy)",
    )
    parser.add_argument("--show", action="store_true", help="Open figures after saving")
    parser.add_argument(
        "--probability",
        action="store_true",
        help="PPSD occurrence heat map; pool all windows, no summary lines",
    )
    parser.add_argument(
        "--cmap",
        default="pqlx",
        help="Default: classic ObsPy pqlx; obspy-default uses installed default",
    )
    parser.add_argument(
        "--max-percentage",
        type=float,
        default=30.0,
        help="Colour scale maximum, percent (ObsPy convention: 30)",
    )
    parser.add_argument(
        "--min-percentage",
        type=float,
        default=0.0,
        help="Hide rarer bins in display only; counts/normalization unchanged",
    )
    parser.add_argument(
        "--curves-only",
        action="store_true",
        help="Every individual saved PSD; no medians, means or percentiles",
    )
    parser.add_argument(
        "--single-color",
        action="store_true",
        help="Use one blue colour for all stations; omit the station legend",
    )
    parser.add_argument(
        "--pressure-limits",
        type=float,
        nargs=2,
        default=(-100.0, 120.0),
        metavar=("LOW_DB", "HIGH_DB"),
        help="Pressure detail y limits; separate full-range PNG retains extremes",
    )
    args = parser.parse_args()
    plot_network(
        args.master_csv,
        network=args.network,
        channels=args.channels,
        output=args.output_dir,
        all_curves=args.all_curves,
        show=args.show,
        pressure_limits=args.pressure_limits,
        curves_only=args.curves_only,
        single_color=args.single_color,
        probability=args.probability,
        cmap=args.cmap,
        max_percentage=args.max_percentage,
        min_percentage=args.min_percentage,
    )


if __name__ == "__main__":
    main()
