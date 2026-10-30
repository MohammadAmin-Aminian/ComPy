"""One figure, one panel per channel: saved PPSD distributions of a station.

Reads registered saved PPSD NPZ files only; occurrence per 1-dB bin per
frequency (equal windows), ObsPy PQLX colours. Pressure panels show the
Brown/Wayne references, seismic panels Peterson NLNM/NHNM.

    python scripts/plot_station_ppsd_panels.py --station RR52 --out fig.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from obspy.imaging.cm import pqlx
from obspy.signal import PPSD

from compy_pressure.plot_network_ppsd import add_references, saved_rows, spectrum
from compy_pressure.plot_ppsd import MASTER_CSV

CHANNELS = ("BDH", "BHZ", "BH1", "BH2")


def occurrence(values, low, high):
    edges = np.arange(np.floor(low), np.ceil(high) + 1)
    counts = np.stack([np.histogram(col, edges)[0] for col in values.T], axis=1)
    return edges, np.ma.masked_equal(100.0 * counts / len(values), 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--network", default="YV")
    parser.add_argument("--station", default="RR52")
    parser.add_argument("--master", default=str(MASTER_CSV))
    parser.add_argument("--out", required=True, help="PNG path; must not exist")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    out = Path(args.out).expanduser()
    if out.exists():
        raise SystemExit(f"Refusing to overwrite {out}")
    rows = [
        r
        for r in saved_rows(args.master, args.network, CHANNELS)
        if r["Seed_id"].split(".")[1] == args.station
    ]
    figure, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
    for axis, channel in zip(axes.flat, CHANNELS, strict=True):
        row = next((r for r in rows if r["Seed_id"].endswith("." + channel)), None)
        if row is None:
            axis.set_title(f"{args.network}.{args.station} {channel}: no saved PPSD")
            axis.axis("off")
            continue
        ppsd = PPSD.load_npz(row["ppsd_npz_file"])
        f, values, pressure = spectrum(ppsd, row)
        low, high = (-100, 120) if pressure else (-200, -50)
        edges, percent = occurrence(np.clip(values, low, high - 1e-9), low, high)
        f_edges = np.concatenate(
            [
                [f[0] ** 1.5 / f[1] ** 0.5],
                np.sqrt(f[:-1] * f[1:]),
                [f[-1] ** 1.5 / f[-2] ** 0.5],
            ]
        )
        mesh = axis.pcolormesh(
            f_edges, edges, percent, cmap=pqlx, vmin=0, vmax=30, rasterized=True
        )
        add_references(axis, pressure, f[0], f[-1])
        axis.set(
            xscale="log",
            xlim=(f[0], f[-1]),
            ylim=(low, high),
            xlabel="Frequency [Hz]",
            ylabel=r"dB re 1 Pa$^2$/Hz" if pressure else r"dB re 1 (m/s$^2$)$^2$/Hz",
            title=f"{ppsd.id}  ({len(values):,} windows)",
        )
        axis.grid(which="both", alpha=0.2)
        axis.legend(loc="lower left", fontsize=7)
        figure.colorbar(mesh, ax=axis, label="Occurrence [%]")
    out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(out, dpi=200)
    print(f"Saved {out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
