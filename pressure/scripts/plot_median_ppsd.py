"""Median pressure PPSD of every station, nominal and with a provisional gain.

Read-only display: saved PPSD NPZs and the saved gain table are never changed.
A gain g multiplies pressure, so the PSD shifts by 20*log10(g) dB at every
frequency. The gains were fitted in 0.03-0.07 Hz only; outside that band the
shifted curve assumes a frequency-independent gain, which is not established.

    conda run -n bruitfm python scripts/plot_median_ppsd.py --gain pole_trial --show
"""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from obspy import read_inventory
from obspy.signal import PPSD

ROOT = Path(__file__).resolve().parents[1]
GAIN_COLUMNS = {
    "pole_trial": "pole_trial_g",
    "baseline": "previous_g",
    "weighted": "weighted_trial_g",
}
CALIBRATION_BAND_HZ = (0.03, 0.07)
# Below the PPSD histogram floor: ObsPy's dtiny (-3076.5 dB) marks zero power,
# e.g. a railed/flat channel. Such medians are dead-channel values, not noise.
FLOOR_DB = -250.0
Y_LIMITS_DB = (-80.0, 110.0)


def pressure_gain_db(gain):
    """An amplitude multiplier changes PSD by 10 log10(g²)."""
    return 20 * np.log10(gain)


def load_gains(table, network, column):
    gains = {}
    with Path(table).open() as stream:
        for row in csv.DictReader(stream):
            if row["network"] == network and row[column]:
                gains[row["station"]] = float(row[column])
    return gains


def load_station_medians(root, network):
    """Exact median over all saved windows; dead channels returned separately."""
    stations = []
    for path in sorted(Path(root, network).glob("*/*DH__*/ppsd.npz")):
        ppsd = PPSD.load_npz(str(path))
        if ppsd.special_handling != "hydrophone":
            continue
        values = np.asarray(ppsd.psd_values, dtype=float)
        frequency = 1 / np.asarray(ppsd.period_bin_centers)
        order = np.argsort(frequency)
        channel = read_inventory(str(path.with_name("inventory.xml")))[0][0][0]
        sensor = (channel.sensor.description if channel.sensor else "") or ""
        stations.append(
            {
                "station": path.parts[-3],
                "seed_id": ppsd.id,
                "family": "DPG" if "DPG" in sensor else "Hydrophone",
                "frequency": frequency[order],
                "median_db": np.median(values, axis=0)[order],
                "windows": len(values),
            }
        )
    if not stations:
        raise ValueError(f"No saved pressure PPSDs under {root}/{network}")
    dead = [s for s in stations if np.median(s["median_db"]) <= FLOOR_DB]
    return [s for s in stations if s not in dead], dead


def on_grid(grid, station, values):
    """Interpolate in log frequency; NaN outside the station's own support."""
    return np.interp(
        np.log(grid),
        np.log(station["frequency"]),
        values,
        left=np.nan,
        right=np.nan,
    )


def colours(stations):
    """Blues for hydrophones, oranges for DPGs; one shade per station."""
    result = {}
    for family, cmap in (("Hydrophone", plt.cm.Blues), ("DPG", plt.cm.Oranges)):
        names = [s["station"] for s in stations if s["family"] == family]
        for name, level in zip(names, np.linspace(0.45, 0.95, max(len(names), 1))):
            result[name] = cmap(level)
    return result


def overview_figure(stations, gains, gain_name, network, dead):
    colour = colours(stations)
    figure = plt.figure(figsize=(17, 13), layout="constrained")
    grid = figure.add_gridspec(3, 2, height_ratios=[3, 1.1, 1.1])
    nominal = figure.add_subplot(grid[0, 0])
    calibrated = figure.add_subplot(grid[0, 1], sharex=nominal, sharey=nominal)
    spread = figure.add_subplot(grid[1, :], sharex=nominal)
    gap = figure.add_subplot(grid[2, :], sharex=nominal)
    for s in stations:
        style = "-" if s["family"] == "DPG" else "--"
        nominal.semilogx(
            s["frequency"],
            s["median_db"],
            style,
            color=colour[s["station"]],
            lw=1.1,
            label=f"{s['station']} ({s['family']})",
        )
        gain = gains.get(s["station"])
        if gain is None:
            calibrated.semilogx(
                s["frequency"], s["median_db"], ":", color="0.7", lw=0.8
            )
        else:
            calibrated.semilogx(
                s["frequency"],
                s["median_db"] + pressure_gain_db(gain),
                style,
                color=colour[s["station"]],
                lw=1.1,
            )
    # Spread and family gap: stations with a gain, on frequencies all support.
    common = [s for s in stations if s["station"] in gains]
    low = max(s["frequency"].min() for s in common)
    high = min(s["frequency"].max() for s in common)
    frequency = np.geomspace(low, high, 400)
    before = np.array([on_grid(frequency, s, s["median_db"]) for s in common])
    shift = np.array([[pressure_gain_db(gains[s["station"]])] for s in common])
    dpg = np.array([s["family"] == "DPG" for s in common])
    for values, name, tone in (
        (before, "Nominal", "0.3"),
        (before + shift, f"With {gain_name} gain", "tab:green"),
    ):
        spread.semilogx(frequency, values.std(axis=0, ddof=1), color=tone, label=name)
        if dpg.any() and (~dpg).any():
            difference = np.median(values[dpg], axis=0) - np.median(
                values[~dpg], axis=0
            )
            gap.semilogx(frequency, difference, color=tone, label=name)
    spread.legend(loc="upper right")
    gap.axhline(0, color="black", lw=0.8)
    gap.legend(loc="upper right")
    spread.set_ylabel("Station-to-station\nstd of medians [dB]")
    spread.set_title(f"Spread across the {len(common)} stations that have a gain")
    gap.set_ylabel("DPG median minus\nhydrophone median [dB]")
    gap.set_xlabel("Frequency [Hz]")
    gap.set_title(
        f"Instrument-family gap ({int(dpg.sum())} DPG vs "
        f"{int((~dpg).sum())} hydrophone stations with gains)"
    )
    for axis in (nominal, calibrated, spread, gap):
        axis.axvspan(*CALIBRATION_BAND_HZ, color="green", alpha=0.08)
        axis.grid(which="both", alpha=0.2)
    nominal.set_ylim(*Y_LIMITS_DB)
    nominal.set_ylabel(r"Median pressure PSD [dB re 1 Pa$^2$/Hz]")
    excluded = ", ".join(s["station"] for s in dead)
    nominal.set_title(
        f"{network}: nominal response ({len(stations)} stations"
        + (f"; dead channel excluded: {excluded}" if dead else "")
        + ")"
    )
    calibrated.set_title(
        f"{network}: provisional {gain_name} gain applied ({len(common)} stations); "
        "grey dotted = no gain"
    )
    nominal.legend(
        fontsize=6.5,
        ncol=3,
        loc="lower left",
        title="dashed = hydrophone, solid = DPG",
        title_fontsize=7,
    )
    figure.suptitle(
        "Median of all saved PPSD windows per station. Green band: Rayleigh fit band "
        "0.03-0.07 Hz. Display only; no gain is applied to saved products."
    )
    return figure


def grid_figure(stations, gains, gain_name, network):
    columns = 7
    rows = int(np.ceil(len(stations) / columns))
    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(3.0 * columns, 2.4 * rows),
        sharex=True,
        sharey=True,
        layout="constrained",
        squeeze=False,
    )
    for axis in axes.flat[len(stations) :]:
        axis.axis("off")
    for axis, s in zip(axes.flat, stations):
        axis.semilogx(
            s["frequency"], s["median_db"], color="black", lw=1, label="Nominal"
        )
        gain = gains.get(s["station"])
        if gain is not None:
            axis.semilogx(
                s["frequency"],
                s["median_db"] + pressure_gain_db(gain),
                color="tab:orange" if s["family"] == "DPG" else "tab:blue",
                lw=1,
                label=f"g={gain:.3f} ({pressure_gain_db(gain):+.1f} dB)",
            )
        axis.axvspan(*CALIBRATION_BAND_HZ, color="green", alpha=0.08)
        axis.set_ylim(*Y_LIMITS_DB)
        axis.set_title(f"{s['station']} {s['family']}", fontsize=9)
        axis.legend(fontsize=6, loc="lower left")
        axis.grid(which="both", alpha=0.2)
    figure.supxlabel("Frequency [Hz]")
    figure.supylabel(r"Median pressure PSD [dB re 1 Pa$^2$/Hz]")
    figure.suptitle(f"{network}: per-station median PPSD, nominal and {gain_name} gain")
    return figure


def band_summary(stations, gains):
    """Hydrophone-vs-DPG level gap and spread in the calibration band."""
    lines = []
    for label, shifted in (("nominal", False), ("calibrated", True)):
        level = {}
        for s in stations:
            if s["station"] not in gains:
                continue
            band = (s["frequency"] >= CALIBRATION_BAND_HZ[0]) & (
                s["frequency"] <= CALIBRATION_BAND_HZ[1]
            )
            value = np.median(s["median_db"][band])
            if shifted:
                value += pressure_gain_db(gains[s["station"]])
            level.setdefault(s["family"], []).append(value)
        values = sum(level.values(), [])
        gap = (
            np.median(level["DPG"]) - np.median(level["Hydrophone"])
            if {"DPG", "Hydrophone"} <= level.keys()
            else float("nan")
        )
        lines.append(
            f"{label:10s}: DPG minus hydrophone {gap:+.2f} dB; "
            f"station std {np.std(values, ddof=1):.2f} dB (n={len(values)})"
        )
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", default="YV")
    parser.add_argument("--results-root", default=ROOT / "results/NETWORKS/TEST")
    parser.add_argument(
        "--gain-table", default=ROOT / "docs/socrates_all_station_comparison.csv"
    )
    parser.add_argument("--gain", choices=GAIN_COLUMNS, default="pole_trial")
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "results/REVIEWS/MEDIAN_PPSD"
    )
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    gains = load_gains(args.gain_table, args.network, GAIN_COLUMNS[args.gain])
    stations, dead = load_station_medians(args.results_root, args.network)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{args.network}_{args.gain}"
    overview = overview_figure(stations, gains, args.gain, args.network, dead)
    overview.savefig(args.output_dir / f"{tag}_median_ppsd_overview.png", dpi=150)
    grid = grid_figure(stations, gains, args.gain, args.network)
    grid.savefig(args.output_dir / f"{tag}_median_ppsd_per_station.png", dpi=120)
    with (args.output_dir / f"{tag}_station_gains.csv").open("w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["station", "family", "windows", "gain", "gain_db"])
        for s in stations + dead:
            gain = gains.get(s["station"])
            writer.writerow(
                [
                    s["station"],
                    s["family"],
                    s["windows"],
                    gain or "",
                    f"{pressure_gain_db(gain):.3f}" if gain else "",
                ]
            )
    print(f"{len(stations)} stations, {len(gains)} with {args.gain} gain")
    if dead:
        print("Excluded dead channels:", ", ".join(s["station"] for s in dead))
    print("\n".join(band_summary(stations, gains)))
    print(f"Saved to {args.output_dir.resolve()}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
