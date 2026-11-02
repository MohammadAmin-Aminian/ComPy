"""Correct tilt and horizontal-coherent noise on response-corrected displacement.

Follows Aminian et al. (2025, sections 4.1.1 and 4.1.4): hourly rotation of the
vertical (TiSKitPy CleanRotator), then removal of the remaining
horizontal-coherent noise with transfer functions estimated over windows of up
to ``--clean-hours`` hours. Hours that cannot be rotated (tilt above
TiSKitPy's limit, gaps) are kept uncorrected and reported.

Example: python examples/2_Rotation.py CC11.2018.160.mseed CC11_clean.mseed
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
from obspy import read

import compy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input", help="Z, 1, 2 (and pressure) traces in displacement [m]"
    )
    parser.add_argument("output", type=Path)
    parser.add_argument("--hours", type=float, default=1, help="rotation window [h]")
    parser.add_argument(
        "--clean-hours", type=float, default=24, help="coherent-noise window [h]"
    )
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    stream = read(args.input)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rotated, azimuth, angle, variance_ratio = compy.correct_tilt(
            stream, args.hours, plot=False, clean_horizontals=False, skip_failed=True
        )
        span_h = (rotated[0].stats.endtime - rotated[0].stats.starttime) / 3600
        cleaned = compy.remove_coherent_noise(
            rotated, min(args.clean_hours, span_h), cover_all=True
        )
    for w in caught:
        print(f"warning: {w.message}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cleaned.write(str(args.output), format="MSEED")
    ok = np.isfinite(angle)
    print(
        f"Rotated {ok.sum()} of {len(angle)} windows (median tilt {np.nanmedian(angle):.2f} deg, "
        f"median variance after/before {np.nanmedian(variance_ratio):.3g}); wrote {args.output}"
    )


if __name__ == "__main__":
    main()
