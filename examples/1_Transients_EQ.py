"""Fit/remove one deployment-specific periodic transient from the vertical trace.

Input must be in units matching the clipping thresholds. Review the fitted
transient before using its cleaned waveform for compliance.
"""

import argparse
from pathlib import Path

from obspy import read
from tiskitpy import PeriodicTransient, TimeSpans


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("output", type=Path)
    parser.add_argument("--period", type=float, required=True, help="period in seconds")
    parser.add_argument("--clip", type=float, nargs=2, required=True)
    parser.add_argument("--period-step", type=float, default=0.05)
    parser.add_argument("--minmag", type=float, default=5.5)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    stream = read(args.input)
    vertical = stream.select(component="Z")
    if len(vertical) != 1:
        raise ValueError("expected exactly one vertical trace")
    z = vertical[0]
    spans = TimeSpans.from_eqs(
        (z.stats.starttime, z.stats.endtime),
        minmag=args.minmag,
        days_per_magnitude=0.5,
        save_eq_file=False,
    )
    transient = PeriodicTransient(
        "periodic", args.period, args.period_step, args.clip, z.stats.starttime
    )
    transient.calc_timing(z, eq_spans=spans)
    transient.calc_transient(z, eq_spans=spans, plots=False)
    cleaned = transient.remove_transient(z, plots=False)
    stream[stream.traces.index(z)] = cleaned
    args.output.parent.mkdir(parents=True, exist_ok=True)
    stream.write(str(args.output), format="MSEED")


if __name__ == "__main__":
    main()
