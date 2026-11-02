"""Remove a periodic transient (glitch) from each channel, avoiding earthquakes.

Unattended by default, with the method of pipeline stage B: earthquake spans
from the USGS catalogue (Deen et al. 2017 rule, 1.5 d per magnitude unit above
``--minmag``) are excluded. The transient period is found by folding the data at
trial periods around ``--period``. A stacked template is subtracted only where
the fold strength reaches ``--min-strength``; otherwise the channel is left
unchanged. Input: response-corrected traces (e.g. the output of
0_Download_data.py); at least several cycles of the transient are needed.

``--interactive`` uses TiSKitPy's PeriodicTransient instead. It asks for the
clip levels at the keyboard, so it needs ``--clip`` and a terminal.

Example: python examples/1_Transients_EQ.py RR52_day.mseed RR52_noglitch.mseed
"""

import argparse
from pathlib import Path

import numpy as np
from obspy import read

import compy.cleaning as cleaning


def earthquake_spans(start, end, minmag, days_per_magnitude):
    from tiskitpy import TimeSpans

    spans = TimeSpans.from_eqs(
        (start, end),
        minmag=minmag,
        days_per_magnitude=days_per_magnitude,
        save_eq_file=False,
    )
    return list(zip(spans.start_times, spans.end_times))


def remove_unattended(stream, args, eq_spans):
    for tr in stream:
        tr.data = tr.data.astype(np.float64)
        detect = tr.copy()
        detect.detrend("linear")
        if not tr.stats.channel.endswith("H"):  # seismometer: detect on acceleration
            detect.differentiate()
            detect.differentiate()
        detect.filter("bandpass", freqmin=0.002, freqmax=0.3, zerophase=True)
        found = cleaning.periodic_template(
            detect,
            period_range=(args.period - args.search, args.period + args.search),
            period_step=args.period_step,
            exclude=eq_spans,
        )
        if found.strength < args.min_strength:
            print(
                f"{tr.id}: no transient (fold strength {found.strength:.2f}); unchanged"
            )
            continue
        disp = tr.copy()
        disp.detrend("linear")
        disp.filter("highpass", freq=0.001, zerophase=True)
        template = cleaning.periodic_template(
            disp,
            period_range=(found.period - 0.05, found.period + 0.05),
            period_step=0.01,
            exclude=eq_spans,
        )
        times = (
            tr.stats.starttime.timestamp
            + np.arange(tr.stats.npts) / tr.stats.sampling_rate
        )
        tr.data -= template.evaluate(times)
        print(
            f"{tr.id}: removed transient, period {template.period:.3f} s, strength {found.strength:.2f}"
        )


def remove_interactive(stream, args):
    from tiskitpy import PeriodicTransient, TimeSpans

    if args.clip is None:
        raise SystemExit("--interactive needs --clip LOW HIGH")
    z = stream.select(component="Z")[0]
    spans = TimeSpans.from_eqs((z.stats.starttime, z.stats.endtime), minmag=args.minmag,
                               days_per_magnitude=0.5, save_eq_file=False)  # fmt: skip
    transient = PeriodicTransient(
        "periodic", args.period, args.period_step, args.clip, z.stats.starttime
    )
    transient.calc_timing(z, eq_spans=spans)
    transient.calc_transient(z, eq_spans=spans, plots=False)
    stream[stream.traces.index(z)] = transient.remove_transient(z, plots=False)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input")
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--period", type=float, default=3600.0, help="expected period [s]"
    )
    parser.add_argument(
        "--search", type=float, default=4.0, help="period search half-width [s]"
    )
    parser.add_argument("--period-step", type=float, default=0.02)
    parser.add_argument(
        "--min-strength", type=float, default=2.0, help="fold-strength threshold"
    )
    parser.add_argument("--minmag", type=float, default=5.85)
    parser.add_argument("--days-per-magnitude", type=float, default=1.5)
    parser.add_argument(
        "--no-eqs", action="store_true", help="skip the catalogue download"
    )
    parser.add_argument(
        "--interactive", action="store_true", help="TiSKitPy PeriodicTransient"
    )
    parser.add_argument(
        "--clip", type=float, nargs=2, help="clip levels for --interactive"
    )
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    stream = read(args.input)
    if args.interactive:
        remove_interactive(stream, args)
    else:
        t0 = min(tr.stats.starttime for tr in stream)
        t1 = max(tr.stats.endtime for tr in stream)
        eqs = (
            []
            if args.no_eqs
            else earthquake_spans(t0, t1, args.minmag, args.days_per_magnitude)
        )
        remove_unattended(stream, args, eqs)
    for tr in stream:
        tr.data = tr.data.astype(np.float32)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    stream.write(str(args.output), format="MSEED", encoding="FLOAT32")


if __name__ == "__main__":
    main()
