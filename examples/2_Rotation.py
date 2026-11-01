"""Correct tilt on response-corrected displacement data."""

import argparse
from pathlib import Path

import numpy as np
from obspy import read
import compy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("output", type=Path)
    parser.add_argument("--hours", type=float, default=1)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    rotated, azimuth, angle, variance_ratio = compy.Rotate(
        read(args.input), args.hours, plot=False
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rotated.write(str(args.output), format="MSEED")
    print(
        f"Cleaned {len(angle)} windows; median after/before variance: {np.median(variance_ratio):.3g}"
    )


if __name__ == "__main__":
    main()
