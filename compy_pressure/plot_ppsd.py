"""Plot saved PSD distributions, their time evolution, and all individual curves."""

import argparse
import csv
import math
from datetime import timezone
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from obspy.signal import PPSD
from obspy.signal.spectral_estimation import get_nhnm, get_nlnm

from compy_pressure.common import PROJECT_DIRECTORY
from compy_pressure.screening import dropout_windows

MASTER_CSV = PROJECT_DIRECTORY / "exports" / "master_inventory.csv"
PLOT_TYPES = ("histogram", "spectrogram", "lines")


def density_label(ppsd):
    pressure = ppsd.special_handling in {"hydrophone", "infrasound"}
    units = r"Pa$^2$/Hz" if pressure else r"m$^2$/s$^4$/Hz"
    return f"Power spectral density [dB re 1 {units}]"


def decorate_density_axis(axis, ppsd, values):
    """Overlay Peterson acceleration models only on seismic PSD axes."""
    frequencies = 1 / ppsd.period_bin_centers
    lower, upper = float(min(frequencies)), float(max(frequencies))
    amplitude_min, amplitude_max = float(np.nanmin(values)), float(np.nanmax(values))
    if ppsd.special_handling is None:
        for label, getter, color, style in [
            ("Low noise (NLNM)", get_nlnm, "black", "-"),
            ("High noise (NHNM)", get_nhnm, "#D55E00", "--"),
        ]:
            periods, levels = getter()
            model_frequencies = 1 / np.asarray(periods)
            keep = (model_frequencies >= lower) & (model_frequencies <= upper)
            model_frequencies, levels = (
                model_frequencies[keep],
                np.asarray(levels)[keep],
            )
            order = np.argsort(model_frequencies)
            axis.plot(
                model_frequencies[order],
                levels[order],
                color=color,
                linestyle=style,
                linewidth=1.8,
                label=label,
                zorder=5,
                path_effects=[
                    path_effects.withStroke(linewidth=3.2, foreground="white")
                ],
            )
            if len(levels):
                amplitude_min = min(amplitude_min, float(levels.min()))
                amplitude_max = max(amplitude_max, float(levels.max()))
        axis.legend(loc="upper right", fontsize=9, framealpha=0.95)
    axis.set_xscale("log")
    axis.set_xlim(lower, upper)
    axis.set_xlabel("Frequency [Hz]")
    axis.set_ylabel(density_label(ppsd))
    axis.set_ylim(
        math.floor(amplitude_min / 10) * 10 - 10,
        math.ceil(amplitude_max / 10) * 10 + 10,
    )


def plot_histogram(ppsd, values):
    frequencies = 1 / ppsd.period_bin_centers
    figure = ppsd.plot(
        show=False,
        xaxis_frequency=True,
        period_lim=(float(min(frequencies)), float(max(frequencies))),
        show_noise_models=False,
    )
    decorate_density_axis(figure.axes[0], ppsd, values)
    return figure


def plot_spectrogram(ppsd, values):
    """Frequency-versus-UTC plot; colors encode PSD amplitude, not probability.

    Each column represents a saved PSD window and spans one processing step.
    Real gaps are left blank; fractional-sample clock phases are tolerated.
    """
    figure, axis = plt.subplots(figsize=(10, 5.5), layout="constrained")
    dates = np.array([time.matplotlib_date for time in ppsd.times_processed])
    timestamps = np.array([float(time) for time in ppsd.times_processed])
    frequency_edges = 1 / ppsd.period_xedges[::-1]
    frequencies = 1 / ppsd.period_bin_centers
    order = np.argsort(frequencies)
    minimum, maximum = float(np.nanmin(values)), float(np.nanmax(values))
    breaks = (
        np.flatnonzero(np.diff(timestamps) > ppsd.step + 1 / ppsd.sampling_rate + 1e-6)
        + 1
    )
    starts, ends = np.r_[0, breaks], np.r_[breaks, len(dates)]
    for start, end in zip(starts, ends, strict=True):
        time_edges = np.r_[dates[start:end], dates[end - 1] + ppsd.step / 86400]
        mesh = axis.pcolormesh(
            time_edges,
            frequency_edges,
            values[start:end, order].T,
            cmap="viridis",
            vmin=minimum,
            vmax=maximum,
            shading="flat",
            rasterized=True,
        )
    figure.colorbar(mesh, ax=axis, label=density_label(ppsd))
    axis.set_facecolor("#E6E6E6")
    axis.set_yscale("log")
    axis.set_ylim(float(min(frequencies)), float(max(frequencies)))
    axis.set_xlim(dates[0], dates[-1] + ppsd.step / 86400)
    locator = mdates.AutoDateLocator(minticks=4, maxticks=8, tz=timezone.utc)
    axis.xaxis.set_major_locator(locator)
    axis.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(locator, tz=timezone.utc)
    )
    axis.set_xlabel("PSD window start time [UTC]")
    axis.set_ylabel("Frequency [Hz]")
    axis.set_title(f"{ppsd.id}: PSD evolution over time")
    return figure


def plot_all_lines(ppsd, values):
    """Draw every saved smoothed PSD curve, with no subsampling or averaging."""
    figure, axis = plt.subplots(figsize=(10, 6), layout="constrained")
    dropped = int(np.isnan(values).all(axis=1).sum())
    values = values[~np.isnan(values).all(axis=1)]
    frequencies = 1 / ppsd.period_bin_centers
    order = np.argsort(frequencies)
    curves = np.empty((len(values), len(frequencies), 2))
    curves[:, :, 0] = frequencies[order]
    curves[:, :, 1] = values[:, order]
    axis.add_collection(
        LineCollection(
            curves,
            colors="#0072B2",
            linewidths=0.45,
            alpha=0.045,
            rasterized=True,
        )
    )
    decorate_density_axis(axis, ppsd, values)
    axis.grid(which="both", alpha=0.25)
    title = f"{ppsd.id}: all {len(values):,} individual PSD curves"
    if dropped:
        title += f" ({dropped:,} constant-data windows left out)"
    axis.set_title(title)
    return figure


def plot_station(
    master_csv=MASTER_CSV,
    *,
    network="YV",
    station="RR52",
    output_directory=None,
    show=True,
    plot_types=PLOT_TYPES,
    selected_rows=None,
):
    """Create up to three plot types per channel from the master CSV's NPZ paths."""
    unknown = set(plot_types) - set(PLOT_TYPES)
    if unknown or not plot_types:
        raise ValueError(
            f"Choose plot types from {PLOT_TYPES}; unknown: {sorted(unknown)}"
        )
    master_csv = Path(master_csv).expanduser().resolve()
    if selected_rows is not None:
        rows = [
            row
            for row in selected_rows.to_dict("records")
            if row.get("ppsd_status") in {"complete", "partial"}
            and row.get("ppsd_npz_file")
        ]
    else:
        with master_csv.open(newline="") as stream:
            rows = _saved_rows(stream, network, station)
    if not rows:
        raise ValueError(f"No saved PPSD results found for {network}.{station}")
    output_directory = (
        Path(output_directory).expanduser().resolve()
        if output_directory
        else PROJECT_DIRECTORY / "results" / "TEST" / network / station / "plots"
    )
    output_directory.mkdir(parents=True, exist_ok=True)
    figures = []
    for row in rows:
        npz_path = Path(row["ppsd_npz_file"])
        ppsd = PPSD.load_npz(str(npz_path))
        if not ppsd.times_processed:
            raise ValueError(f"No PSD segments in {npz_path}")
        values = np.asarray(ppsd.psd_values, dtype=float)
        dropouts = dropout_windows(ppsd)
        if dropouts.any() and not dropouts.all():
            # Constant-data windows stay in the NPZ; plots show them as gaps.
            values[dropouts] = np.nan
        for kind in dict.fromkeys(plot_types):
            factory, suffix = {
                "histogram": (plot_histogram, "ppsd"),
                "spectrogram": (plot_spectrogram, "spectrogram"),
                "lines": (plot_all_lines, "all_psd_lines"),
            }[kind]
            figure = factory(ppsd, values)
            path = output_directory / f"{npz_path.parent.name}__{suffix}.png"
            figure.savefig(path, dpi=160)
            print(f"{ppsd.id}: {kind} plot saved to {path}", flush=True)
            figures.append(figure)
            if not show:
                plt.close(figure)
    if show:
        plt.show()
    return figures


def _saved_rows(stream, network, station):
    """Select a station's available results from a CSV stream."""
    reader = csv.DictReader(stream)
    required = {"Seed_id", "ppsd_npz_file", "ppsd_status"}
    if not required.issubset(reader.fieldnames or []):
        raise ValueError("Run ppsd.py first to populate the master PPSD results")
    rows = [
        row
        for row in reader
        if row["Seed_id"].startswith(f"{network}.{station}.")
        and row["ppsd_status"] in {"complete", "partial"}
        and row["ppsd_npz_file"]
    ]
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-csv", type=Path, default=MASTER_CSV)
    parser.add_argument("--network", default="YV")
    parser.add_argument("--station", default="RR52")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--plot-type", action="append", choices=PLOT_TYPES)
    parser.add_argument(
        "--no-show", action="store_true", help="Save PNGs without windows"
    )
    args = parser.parse_args()
    plot_station(
        args.master_csv,
        network=args.network,
        station=args.station,
        output_directory=args.output_dir,
        show=not args.no_show,
        plot_types=args.plot_type or PLOT_TYPES,
    )


if __name__ == "__main__":
    main()
