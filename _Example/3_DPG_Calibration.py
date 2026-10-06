"""Estimate gauge gain from raw data and a local StationXML response inventory."""

import argparse

from obspy import read, read_inventory
import Pressure_calibration as dpg


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input", help="raw-count MiniSEED, not displacement output from example 0"
    )
    parser.add_argument(
        "inventory", help="StationXML containing vertical and pressure responses"
    )
    parser.add_argument("--minmag", type=float, default=7)
    parser.add_argument("--band", type=float, nargs=2, default=[0.03, 0.07])
    args = parser.parse_args()
    gain = dpg.calculate_spectral_ratio(
        read(args.input),
        inventory=read_inventory(args.inventory),
        mag=args.minmag,
        f_min=args.band[0],
        f_max=args.band[1],
    )
    print(f"Pressure gain: {gain:.6g}")


if __name__ == "__main__":
    main()
