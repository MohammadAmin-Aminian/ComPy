"""Histogram of provisional Rayleigh pressure gains; read-only display.

g multiplies recorded pressure (PSD shift 20*log10 g). Published YL values are
1/g in the Zha & Webb convention and are not plotted here.

    conda run -n bruitfm python scripts/plot_gain_histogram.py --gain pole_trial --show
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from obspy import read_inventory

ROOT = Path(__file__).resolve().parents[1]
GAIN_COLUMNS = {
    "pole_trial": "pole_trial_g",
    "baseline": "previous_g",
    "weighted": "weighted_trial_g",
}
# Validated categorical slots in fixed order (light mode reference palette).
NETWORK_COLOURS = {"YV": "#2a78d6", "YL": "#eb6834", "7D": "#1baf7a", "XO": "#eda100"}
FAMILY_COLOURS = {"Hydrophone (DEPAS)": "#e87ba4", "DPG (INSU)": "#008300"}
BIN_WIDTH = 0.05


def yv_family(station, root):
    paths = sorted(Path(root, "YV", station).glob("*DH__*/inventory.xml"))
    if not paths:
        return "unknown"
    channel = read_inventory(str(paths[0]))[0][0][0]
    text = (channel.sensor.description if channel.sensor else "") or ""
    return "DPG (INSU)" if "DPG" in text else "Hydrophone (DEPAS)"


def stacked(axis, groups, colours, bins):
    bottom = np.zeros(len(bins) - 1)
    for name, values in groups.items():
        counts, _ = np.histogram(values, bins)
        axis.bar(
            bins[:-1],
            counts,
            width=BIN_WIDTH,
            bottom=bottom,
            align="edge",
            color=colours[name],
            edgecolor="white",
            linewidth=1.0,
            label=f"{name}: n={len(values)}, median {np.median(values):.3f}",
        )
        bottom += counts
    for name, values in groups.items():
        axis.axvline(np.median(values), color=colours[name], lw=1.5, ls="--")
    axis.axvline(1.0, color="black", lw=1, ls=":", label="g = 1 (nominal response)")
    axis.set_ylabel("Stations")
    axis.yaxis.get_major_locator().set_params(integer=True)
    axis.grid(axis="y", alpha=0.25)
    axis.legend(fontsize=8, loc="upper right", frameon=False)
    top = axis.secondary_xaxis(
        "top",
        functions=(
            lambda g: 20 * np.log10(np.clip(g, 1e-6, None)),
            lambda db: 10 ** (db / 20),
        ),
    )
    top.set_xlabel("PSD shift 20·log10(g) [dB]", fontsize=8)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gain", choices=GAIN_COLUMNS, default="pole_trial")
    parser.add_argument(
        "--table", default=ROOT / "docs/socrates_all_station_comparison.csv"
    )
    parser.add_argument("--results-root", default=ROOT / "results/NETWORKS/TEST")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    output = args.output or (
        ROOT / f"results/REVIEWS/MEDIAN_PPSD/gain_histogram_{args.gain}.png"
    )
    column = GAIN_COLUMNS[args.gain]
    table = pd.read_csv(args.table).dropna(subset=[column])
    gains = table[column].to_numpy()
    bins = (
        np.arange(
            np.floor(gains.min() / BIN_WIDTH), np.ceil(gains.max() / BIN_WIDTH) + 1
        )
        * BIN_WIDTH
    )
    yv = table[table.network == "YV"].copy()
    yv["family"] = [yv_family(s, args.results_root) for s in yv.station]
    families = {
        name: yv.loc[yv.family == name, column].to_numpy()
        for name in FAMILY_COLOURS
        if (yv.family == name).any()
    }
    networks = {
        name: table.loc[table.network == name, column].to_numpy()
        for name in NETWORK_COLOURS
        if (table.network == name).any()
    }
    figure, (top, bottom) = plt.subplots(
        2, 1, figsize=(10, 8.5), sharex=True, layout="constrained"
    )
    stacked(top, families, FAMILY_COLOURS, bins)
    top.set_title(f"YV: {len(yv)} stations by pressure-sensor family")
    stacked(bottom, networks, NETWORK_COLOURS, bins)
    bottom.set_title(f"All networks: {len(table)} stations")
    bottom.set_xlabel(
        "Pressure multiplier g (corrected pressure = g × recorded pressure)"
    )
    figure.suptitle(
        f"Provisional Rayleigh gains ({args.gain}); "
        "not applied, not accepted calibrations",
        fontsize=11,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=150)
    summary = (
        table.groupby("network")[column]
        .describe(percentiles=[0.25, 0.5, 0.75])[
            ["count", "min", "25%", "50%", "75%", "max"]
        ]
        .round(3)
    )
    print(summary.to_string())
    for name, values in families.items():
        print(
            f"YV {name}: n={len(values)} median {np.median(values):.3f} "
            f"range {values.min():.3f}-{values.max():.3f}"
        )
    print(f"Saved {output.resolve()}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
