"""One figure: the median saved PSD of every channel in a network, with limits.

Reads the registered saved PPSD NPZ files only (no recalculation, no gain, no
master update). Each curve is the median over all saved windows of one
channel epoch, the same median as ``compy_pressure.plot_network_ppsd``. Pressure
channels get the project's reference limits from ``reference_data``:

* Brown et al. (2014) IDC hydroacoustic low/high (green, at and above 0.2 Hz):
  approximate digitization of in-water IMS hydrophone curves, drawn only
  inside the digitized support;
* Webb (1998) Figure 2 low/high points (red, at and below 0.2 Hz): drawn at
  their four frequencies only, never interpolated.

These are display references, not seafloor acceptance bounds. Seismic channels
get Peterson NLNM/NHNM instead. A channel whose median lies entirely outside
the plotted dB range (for example RR12's flat -3077 dB pressure archive) is
not drawn and is named as off scale in the legend and the report.

Example::

    python scripts/plot_network_median_ppsd.py --network YV --channel BDH --show
"""

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patheffects
from matplotlib.lines import Line2D
from obspy.signal import PPSD
from obspy.signal.spectral_estimation import get_nhnm, get_nlnm

from compy_pressure.plot_network_ppsd import saved_rows, sha256, spectrum
from compy_pressure.plot_ppsd import MASTER_CSV
from compy_pressure.pressure_test import REFERENCE_COLORS, load_references

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

DEFAULT_YLIM = {True: (-80.0, 100.0), False: (-200.0, -40.0)}
HALO = [patheffects.withStroke(linewidth=5.5, foreground="white")]
WEBB_COLOR = REFERENCE_COLORS["webb_instructions"]


def median_curves(rows, stations=None):
    """Median over saved windows per channel epoch, with input hashes."""
    curves = []
    for row in rows:
        station = row["Seed_id"].split(".")[1]
        if stations and station not in stations:
            continue
        path = Path(row["ppsd_npz_file"])
        if not path.is_file():
            print(f"Missing saved PPSD, skipped: {path}")
            continue
        before = sha256(path)
        ppsd = PPSD.load_npz(str(path))
        f, values, pressure = spectrum(ppsd, row)
        if sha256(path) != before:
            raise RuntimeError(f"Saved PPSD changed while reading: {path}")
        curves.append(
            {
                "station": station,
                "seed_id": row["Seed_id"],
                "channel_start_time": row.get("channel_start_time", ""),
                "ppsd_npz_file": str(path),
                "sha256": before,
                "windows": len(values),
                "pressure": pressure,
                "frequency_hz": f,
                "median_db": np.median(values, axis=0),
            }
        )
    if len({c["pressure"] for c in curves}) > 1:
        raise ValueError("Cannot mix pressure and seismic units in one figure")
    return curves


def off_scale(median, ylim):
    """True when no part of the curve is inside the plotted dB range."""
    median = np.asarray(median, dtype=float)
    return bool(np.all((median < ylim[0]) | (median > ylim[1])))


def draw_pressure_limits(axis, lower, upper):
    """Brown lines and Webb (1998) low-high brackets, on top of the data."""
    metadata, brown = load_references()
    brown_lower = max(lower, metadata["brown"].get("minimum_comparison_hz", 0))
    keep = (brown[:, 0] >= brown_lower) & (brown[:, 0] <= upper)
    handles = []
    for column, name, style in [(1, "low", "--"), (2, "high", "-")]:
        (line,) = axis.plot(
            brown[keep, 0],
            brown[keep, column],
            color=REFERENCE_COLORS["brown_figure8"],
            ls=style,
            lw=3,
            zorder=20,
            path_effects=HALO,
            label=f"Brown et al. (2014) {name}: in-water IMS ref., approx. digitized",
        )
        handles.append(line)
    webb = metadata["webb"]
    points = [
        p
        for p in webb["points"]
        if webb.get("use_in_plot") and lower <= p["frequency_hz"] <= upper
    ]
    if points:
        label = webb.get("display_label", "Webb 1998 Fig. 2")
        for p in points:
            axis.plot(
                [p["frequency_hz"]] * 2,
                [p["low_db_pa2_hz"], p["high_db_pa2_hz"]],
                color=WEBB_COLOR,
                lw=2.5,
                zorder=21,
                path_effects=HALO,
            )
        for bound, marker in [("low", "v"), ("high", "^")]:
            (mark,) = axis.plot(
                [p["frequency_hz"] for p in points],
                [p[f"{bound}_db_pa2_hz"] for p in points],
                ls="none",
                marker=marker,
                ms=12,
                color=WEBB_COLOR,
                markeredgecolor="white",
                markeredgewidth=1.5,
                zorder=22,
                label=f"{label} {bound} (points only)",
            )
            handles.append(mark)
    return handles


def draw_seismic_limits(axis, lower, upper):
    handles = []
    for getter, label, style in [
        (get_nlnm, "Peterson NLNM", "--"),
        (get_nhnm, "Peterson NHNM", "-"),
    ]:
        periods, db = getter()
        f = 1 / np.asarray(periods)
        keep = (f >= lower) & (f <= upper)
        order = np.argsort(f[keep])
        (line,) = axis.plot(
            f[keep][order],
            np.asarray(db)[keep][order],
            color="black",
            ls=style,
            lw=3,
            zorder=20,
            path_effects=HALO,
            label=label,
        )
        handles.append(line)
    return handles


def plot(curves, network, channel, ylim, output, show=False):
    pressure = curves[0]["pressure"]
    ylim = tuple(ylim) if ylim else DEFAULT_YLIM[pressure]
    lower = min(c["frequency_hz"][0] for c in curves)
    upper = max(c["frequency_hz"][-1] for c in curves)
    colors = plt.get_cmap("turbo")(np.linspace(0.02, 0.98, len(curves)))
    figure, axis = plt.subplots(figsize=(16, 11))
    station_handles, flagged = [], []
    for curve, color in zip(curves, colors, strict=True):
        name = curve["station"]
        if sum(c["station"] == name for c in curves) > 1:
            name += f" ({curve['channel_start_time'][:10]})"
        if off_scale(curve["median_db"], ylim):
            flagged.append(curve)
            level = float(np.median(curve["median_db"]))
            station_handles.append(
                Line2D(
                    [],
                    [],
                    color=color,
                    lw=1.6,
                    ls=":",
                    label=f"{name}: off scale ({level:.0f} dB)",
                )
            )
            continue
        (line,) = axis.plot(
            curve["frequency_hz"],
            curve["median_db"],
            color=color,
            lw=1.2,
            alpha=0.9,
            label=name,
        )
        station_handles.append(line)
    draw = draw_pressure_limits if pressure else draw_seismic_limits
    limit_handles = draw(axis, lower, upper)
    axis.set(
        xscale="log",
        xlim=(lower, upper),
        ylim=ylim,
        xlabel="Frequency [Hz]",
        ylabel=(
            r"Pressure PSD [dB re 1 Pa$^2$/Hz]"
            if pressure
            else r"Acceleration PSD [dB re 1 (m/s$^2$)$^2$/Hz]"
        ),
    )
    axis.grid(which="both", alpha=0.25)
    windows = sum(c["windows"] for c in curves)
    axis.set_title(
        f"{network} {channel}: median of saved PSD windows per channel "
        f"({len(curves)} channels, {windows:,} windows; no gain applied)",
        loc="left",
        fontsize=12,
    )
    # Legends sit above the axes, anchored top-left, so no curve is covered.
    # The axes top is then set from the measured legend heights.
    left = 0.07
    first = figure.legend(
        handles=limit_handles,
        loc="upper left",
        bbox_to_anchor=(left, 0.995),
        ncol=2,
        fontsize=10,
        title="Limits",
        title_fontproperties={"weight": "bold"},
        frameon=True,
    )
    figure.add_artist(first)
    inverse = figure.transFigure.inverted()
    renderer = figure.canvas.get_renderer()
    below = first.get_window_extent(renderer).transformed(inverse).y0
    second = figure.legend(
        handles=station_handles,
        loc="upper left",
        bbox_to_anchor=(left, below - 0.008),
        ncol=10,
        fontsize=8.5,
        title="Channels (median PSD)",
        frameon=True,
    )
    bottom = second.get_window_extent(renderer).transformed(inverse).y0
    figure.subplots_adjust(left=left, right=0.98, bottom=0.07, top=bottom - 0.045)
    output.mkdir(parents=True, exist_ok=False)
    stem = f"{network}_{channel}_median_ppsd"
    figure.savefig(output / f"{stem}.png", dpi=200)
    figure.savefig(output / f"{stem}.pdf")
    if show:
        plt.show()
    plt.close(figure)
    return ylim, flagged, stem


def write_outputs(curves, output, stem, report):
    with open(output / f"{stem}.csv", "x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["seed_id", "channel_start_time", "frequency_hz", "median_db"])
        for c in curves:
            for f, db in zip(c["frequency_hz"], c["median_db"], strict=True):
                writer.writerow(
                    [c["seed_id"], c["channel_start_time"], f"{f:.8g}", f"{db:.4f}"]
                )
    with open(output / f"{stem}.json", "x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--network", default="YV")
    parser.add_argument(
        "--channel", default="BDH", help="BDH (pressure), BHZ, BH1 or BH2"
    )
    parser.add_argument(
        "--master", default=str(MASTER_CSV), help="master_inventory.csv"
    )
    parser.add_argument("--stations", nargs="+", help="optional subset, e.g. RR17 RR52")
    parser.add_argument("--ylim", nargs=2, type=float, metavar=("LOW", "HIGH"))
    parser.add_argument("--out", help="new folder (default under results/REVIEWS)")
    parser.add_argument("--show", action="store_true", help="also open the figure")
    args = parser.parse_args()
    master = Path(args.master).expanduser().resolve()
    if not master.is_file():
        raise SystemExit(f"Master inventory not found: {master} (pass --master)")
    project = master.parent.parent
    output = (
        Path(args.out)
        if args.out
        else (
            project
            / "results/REVIEWS/PPSD_MEDIAN"
            / f"run_{datetime.now(UTC):%Y%m%dT%H%M%SZ}_{uuid4().hex[:8]}"
        )
    )
    output = output.expanduser().resolve()
    if master.parent == output or master.parent in output.parents:
        raise SystemExit("Output must not be inside the production exports folder")
    master_before = sha256(master)
    rows = saved_rows(master, args.network, (args.channel,))
    curves = median_curves(rows, set(args.stations or []))
    if not curves:
        raise SystemExit(f"No saved {args.network} {args.channel} PPSDs found")
    ylim, flagged, stem = plot(
        curves, args.network, args.channel, args.ylim, output, args.show
    )
    if sha256(master) != master_before:
        raise RuntimeError("Master inventory changed while plotting")
    report = {
        "network": args.network,
        "channel": args.channel,
        "created_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "statistic": "median over all saved windows per channel epoch",
        "master": str(master),
        "master_sha256": master_before,
        "ylim_db": list(ylim),
        "gain_applied": False,
        "off_scale": [c["seed_id"] for c in flagged],
        "inputs": [
            {
                k: c[k]
                for k in (
                    "seed_id",
                    "channel_start_time",
                    "ppsd_npz_file",
                    "sha256",
                    "windows",
                )
            }
            for c in curves
        ],
    }
    write_outputs(curves, output, stem, report)
    for c in flagged:
        print(f"Off scale (not drawn): {c['seed_id']}")
    print(f"Figure: {output / (stem + '.png')}")


if __name__ == "__main__":
    main()
