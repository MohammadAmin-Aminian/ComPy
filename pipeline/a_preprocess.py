"""Stage A: decimate raw OBS data to 1.25 Hz and remove instrument responses.

For each station and day this writes two MiniSEED files under work/<STA>/:
  counts/<STA>.<YYYY.JJJ>.mseed   decimated raw counts (BHZ, BH1, BH2, BDH)
  phys/<STA>.<YYYY.JJJ>.mseed     BHZ/BH1/BH2 displacement [m], BDH pressure [Pa]

Choices (documented in the report):
  * Decimation to 1.25 Hz with TiSKitPy SAC FIR filters (config.FACTORS;
    5, 5, 2 from 62.5 Hz as in Aminian 2024, thesis section 3.1.1). Channel codes are restored after
    decimation so the StationXML responses still match.
  * Each raw segment is trimmed to start on the 0.8 s output grid (residual
    offset < one raw sample, identical on all channels of the digitizer).
  * Responses are removed after decimation with water_level=None and an
    explicit pre-filter. With obspy's default 60 dB water level, a displacement
    response removed at 62.5 Hz is clipped below ~0.03 Hz (tested: 22 dB error).
  * Chunks of CHUNK_DAYS days are padded by PAD_S seconds on each side so FIR
    and response edge effects fall outside the kept day.
  * Saturated or dead channels (Aminian et al. 2025, section 4: windows in
    which any seismometer component is saturated are eliminated). Raw counts
    are scanned in QC_BLOCK_S blocks; a block with any sample at or above
    CLIP_COUNTS, or with a standard deviation below DEAD_STD_COUNTS, is bad on
    every channel, widened by QC_MARGIN_S. All channels are cut there before
    decimation, so response ringing stays inside the trimmed segment edges.
    The cut spans are listed in work/<STA>/qc/<STA>.<YYYY.JJJ>.json.
  * Raw traces whose sampling rate is not OUT_RATE times the product of
    config.FACTORS are left out and listed under "wrong_rate" in the QC file
    (ComPy 3.2.1; before, a day at another rate was appended to its neighbour
    sample by sample and decimated with the same factors).
  * A segment that crosses a StationXML epoch of its channel is split there, so
    that each part is corrected with its own response; the times are listed
    under "response_epochs" (ComPy 3.2.1; before, the response valid at the
    start of the segment was applied to all of it).
  * A station without any day file in the archive ends the run with an error.
"""

import argparse
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction

import numpy as np
from obspy import Stream, UTCDateTime, read, read_inventory

from config import (
    CHANNELS,
    FACTORS,
    INVENTORY,
    NETWORK,
    OUT_RATE,
    P,
    SDS,
    SDS_LOCATION,
    WORK,
    YEARS,
    Z,
)

#: Sampling rate [Hz] that config.FACTORS decimates to OUT_RATE.
EXPECTED_RATE = OUT_RATE * float(np.prod(FACTORS))
#: Relative tolerance on the rate of a raw trace (clock drift is a few 1e-8).
RATE_TOLERANCE = 1e-3
#: Frequencies [Hz] at which two response epochs are compared.
EPOCH_TEST_FREQUENCIES = np.array([0.002, 0.01, 0.05, 0.2, 0.5])

GRID = 1 / OUT_RATE
CHUNK_DAYS = 6
PAD_S = 1800.0
MIN_SEGMENT_S = 2 * PAD_S
PRE_FILT = (0.0008, 0.0015, 0.45, 0.6)
DAY = 86400.0
CLIP_COUNTS = 0.99 * 2**23  # 24-bit digitizer full scale
CLIP_COUNTS_32 = 0.99 * 2**31  # 32-bit digitizers and absolute pressure gauges
DEAD_STD_COUNTS = 5.0
QC_BLOCK_S = 600.0
QC_MARGIN_S = 300.0


def day_files(station, channel, day):
    """Path of one channel-day file in the SDS archive."""
    path = SDS / str(day.year) / NETWORK / station / f"{channel}.D"
    return (
        path
        / f"{NETWORK}.{station}.{SDS_LOCATION}.{channel}.D.{day.year}.{day.julday:03d}"
    )


def station_days(station):
    """Valid days with a vertical file for a station."""
    days = set()
    for year in YEARS:
        folder = SDS / str(year) / NETWORK / station / f"{Z}.D"
        if folder.is_dir():
            for f in folder.iterdir():
                if f.name.startswith("."):
                    continue
                try:
                    UTCDateTime(year=year, julday=int(f.name[-3:]))
                except ValueError:  # e.g. day 366 of a non-leap year, or a partial copy
                    print(f"SKIP invalid day file {f}", file=sys.stderr)
                    continue
                days.add((year, int(f.name[-3:])))
    return [UTCDateTime(year=y, julday=j) for y, j in sorted(days)]


def grid_ceil(t):
    """First point of the output sample grid at or after a time."""
    return UTCDateTime(np.ceil(t.timestamp / GRID) * GRID)


def join_day_files(traces, tolerance=1.5):
    """Join consecutive day-file traces of one channel into gap-free segments.

    Day files overlap by ~1 min and their sample grids slide by ~2.4 ms/day
    (clock correction). When the offset wraps past a whole sample, obspy's
    merge sees a one-sample conflict and masks the overlap. Here overlapping
    samples are dropped and traces whose junction is within ``tolerance``
    samples are concatenated (timing slip < 1 sample = 16 ms, removed again
    when snapping to the 0.8 s output grid); larger gaps start a new segment,
    and so does a trace with another sampling rate.
    """
    traces = sorted(traces, key=lambda tr: tr.stats.starttime)
    segments = []
    for tr in traces:
        tr = tr.copy()
        if not segments:
            segments.append(tr)
            continue
        prev = segments[-1]
        delta = prev.stats.delta
        expected = prev.stats.endtime + delta
        if tr.stats.endtime < expected - delta / 2:
            continue  # wholly inside the previous trace
        if tr.stats.starttime < expected - delta / 2:
            tr = tr.slice(expected - delta / 2, tr.stats.endtime)
        same_rate = abs(tr.stats.delta - delta) <= 1e-6 * delta
        if same_rate and abs(tr.stats.starttime - expected) <= tolerance * delta:
            prev.data = np.concatenate([prev.data, tr.data])
        else:
            segments.append(tr)
    return segments


def clip_threshold(data, width):
    """Clip level in counts of one raw channel.

    CLIP_COUNTS assumes 24-bit counts. A channel that exceeds 2^23 in more than
    half of its blocks cannot be 24-bit: it comes from a 32-bit digitizer or
    carries a large DC level (absolute pressure gauge), so CLIP_COUNTS_32 applies.
    A 24-bit channel stuck at full scale stays at or below 2^23 and is still clipped.
    """
    n = data.size // width
    if not n:
        return CLIP_COUNTS
    peaks = np.abs(data[: n * width]).reshape(n, width).max(axis=1)
    return CLIP_COUNTS_32 if np.mean(peaks > 2**23) > 0.5 else CLIP_COUNTS


def bad_spans(trace):
    """Clipped or flat QC_BLOCK_S blocks of one raw trace, widened by QC_MARGIN_S."""
    width = int(round(QC_BLOCK_S * trace.stats.sampling_rate))
    n_blocks = -(-trace.stats.npts // width)
    data = np.asarray(trace.data, dtype=np.float64)
    limit = clip_threshold(data, width)
    spans = []
    for i in range(n_blocks):
        block = data[i * width : (i + 1) * width]
        if block.size < 2:
            continue
        clipped = np.abs(block).max() >= limit
        dead = block.std() < DEAD_STD_COUNTS
        if not (clipped or dead):
            continue
        start = trace.stats.starttime + i * QC_BLOCK_S - QC_MARGIN_S
        end = trace.stats.starttime + (i + 1) * QC_BLOCK_S + QC_MARGIN_S
        reason = f"{trace.stats.channel}:{'clip' if clipped else 'dead'}"
        if spans and start <= spans[-1][1]:
            spans[-1][1] = max(spans[-1][1], end)
            spans[-1][2].add(reason)
        else:
            spans.append([start, end, {reason}])
    return spans


def merge_spans(spans):
    """Merge overlapping bad spans and their reasons."""
    merged = []
    for start, end, reasons in sorted(spans, key=lambda s: s[0]):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
            merged[-1][2] |= reasons
        else:
            merged.append([start, end, set(reasons)])
    return merged


def cut_spans(trace, spans):
    """Pieces of ``trace`` outside the sorted, merged ``spans``."""
    pieces, cursor, end = [], trace.stats.starttime, trace.stats.endtime
    for start, stop, _ in spans:
        if stop <= cursor or start >= end:
            continue
        # Spans of different channels can start a few milliseconds apart: a gap
        # shorter than one sample holds no piece (slicing it raised ValueError
        # and lost the whole chunk, XE.CC10 and XE.EC05).
        if start - trace.stats.delta > cursor:
            pieces.append(trace.slice(cursor, start - trace.stats.delta).copy())
        cursor = max(cursor, stop)
    if cursor < end:
        pieces.append(trace.slice(cursor, end).copy())
    return [p for p in pieces if p.stats.npts > 1]


class WrongRate(ValueError):
    """A raw trace whose sampling rate config.FACTORS does not bring to OUT_RATE."""


def has_expected_rate(trace):
    """True when ``trace`` has the sampling rate that config.FACTORS expects."""
    return (
        abs(trace.stats.sampling_rate - EXPECTED_RATE) <= RATE_TOLERANCE * EXPECTED_RATE
    )


def _response_at(inventory, trace, time):
    """Complex response of the channel of ``trace`` at ``time``, or None outside every epoch."""
    try:
        response = inventory.get_response(trace.id, time)
    except Exception:  # obspy raises a bare Exception when no epoch matches
        return None
    output = "DEF" if trace.stats.channel == P else "DISP"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return response.get_evalresp_response_for_frequencies(
            EPOCH_TEST_FREQUENCIES, output=output
        )


def epoch_boundaries(inventory, trace):
    """Times inside ``trace`` at which the response of its channel changes.

    A StationXML epoch that starts or ends inside the trace counts when the
    response differs across it, or when one side has no response at all.
    Epochs that repeat the same response (administrative splits) do not.
    """
    stats = trace.stats
    times = set()
    for network in inventory.select(
        network=stats.network, station=stats.station, channel=stats.channel
    ):
        for station in network:
            for channel in station:
                if channel.location_code != stats.location:
                    continue
                for t in (channel.start_date, channel.end_date):
                    if t is not None and stats.starttime < t <= stats.endtime:
                        times.add(t.timestamp)
    changes = []
    for timestamp in sorted(times):
        t = UTCDateTime(timestamp)
        before = _response_at(inventory, trace, t - stats.delta)
        after = _response_at(inventory, trace, t + stats.delta)
        same = (
            before is not None
            and after is not None
            and np.allclose(before, after, rtol=1e-9, atol=0.0)
        )
        if not same:
            changes.append(t)
    return changes


def split_at(trace, times):
    """Pieces of ``trace`` between ``times``.

    A sample exactly at a boundary stays in the earlier piece: where one epoch
    ends at the instant the next starts, obspy returns the earlier response.
    """
    stats = trace.stats
    edges = [0]
    for t in times:
        edges.append(int(np.floor((t - stats.starttime) / stats.delta + 1e-6)) + 1)
    edges.append(stats.npts)
    pieces = []
    for first, last in zip(edges[:-1], edges[1:]):
        if last <= first:
            continue
        piece = trace.copy()
        piece.data = trace.data[first:last].copy()
        piece.stats.starttime = stats.starttime + first * stats.delta
        pieces.append(piece)
    return pieces


def decimate_segment(trace):
    """Decimate one gap-free raw trace on the 0.8 s grid; returns None if too short.

    Raises:
        WrongRate: The trace is not sampled at OUT_RATE times the product of
            config.FACTORS. A driver that decimates several rates replaces this
            function and decides for itself.
    """
    from tiskitpy import Decimator

    if not has_expected_rate(trace):
        raise WrongRate(
            f"{trace.stats.sampling_rate:g} Hz, expected {EXPECTED_RATE:g} Hz "
            f"(config.FACTORS {tuple(FACTORS)} to {OUT_RATE:g} Hz)"
        )

    # Rates a few 1e-8 off nominal would give segments that obspy refuses to merge.
    nominal = float(Fraction(trace.stats.sampling_rate).limit_denominator(100))
    if abs(nominal - trace.stats.sampling_rate) < 1e-5 * nominal:
        trace.stats.sampling_rate = nominal
    start = grid_ceil(trace.stats.starttime)
    if start > trace.stats.endtime:  # a piece that ends before the next grid point
        return None
    trace = trace.slice(start, trace.stats.endtime)
    if trace.stats.npts * trace.stats.delta < MIN_SEGMENT_S:
        return None
    offset = trace.stats.starttime - start
    trace.data = trace.data.astype(np.float64)
    out = Decimator(list(FACTORS)).decimate(Stream([trace]), keep_dtype=False)[0]
    out.stats.channel = trace.stats.channel
    out.stats.starttime = out.stats.starttime - offset  # snap (< 1 raw sample) to grid
    return out


def process_chunk(station, chunk_days):
    """Process one chunk of consecutive days and write its daily files."""
    import logging

    logging.disable(logging.INFO)
    inv = read_inventory(str(INVENTORY)).select(station=station)
    t0 = chunk_days[0] - PAD_S
    t1 = chunk_days[-1] + DAY + PAD_S
    n_days = int(round((chunk_days[-1] - chunk_days[0]) / DAY)) + 1
    span = [chunk_days[0] + i * DAY for i in range(-1, n_days + 1)]
    segments = []
    for channel in CHANNELS:
        raw = Stream()
        for day in span:
            path = day_files(station, channel, day)
            if path.exists():
                try:
                    raw += read(str(path), starttime=t0, endtime=t1)
                except Exception as exc:  # corrupt file: report and continue
                    print(f"READ-FAIL {path}: {exc}", file=sys.stderr)
        raw = raw.split()  # masked gaps inside a file become separate traces
        segments.extend(join_day_files(list(raw)))
    if not segments:
        return station, [], "no data"
    bad = merge_spans([span for tr in segments for span in bad_spans(tr)])
    segments = [piece for tr in segments for piece in cut_spans(tr, bad)]
    qc_file = (
        WORK
        / station
        / "qc"
        / f"{station}.{chunk_days[0].year}.{chunk_days[0].julday:03d}.json"
    )
    qc_file.parent.mkdir(parents=True, exist_ok=True)
    qc = {
        "days": [str(d) for d in chunk_days],
        "cut": [[str(a), str(b), sorted(r)] for a, b, r in bad],
        "cut_hours": sum(b - a for a, b, _ in bad) / 3600,
        "rule": f"|counts| >= {CLIP_COUNTS:.0f} (24-bit) or {CLIP_COUNTS_32:.0f} (other channels), "
                f"or std < {DEAD_STD_COUNTS} in {QC_BLOCK_S:.0f}-s blocks, +/- {QC_MARGIN_S:.0f} s",
    }  # fmt: skip
    qc_file.write_text(json.dumps(qc, indent=1))
    counts, phys = Stream(), Stream()
    for trace in segments:
        try:
            decimated = decimate_segment(trace)
        except WrongRate as exc:
            # Leave the trace out and say so: decimating it with factors meant
            # for another rate would write data on a wrong time axis.
            print(f"WRONG-RATE {trace.id} {trace.stats.starttime} - "
                  f"{trace.stats.endtime}: {exc}; left out", file=sys.stderr)  # fmt: skip
            qc.setdefault("wrong_rate", []).append(
                [trace.stats.channel, str(trace.stats.starttime),
                 str(trace.stats.endtime), float(trace.stats.sampling_rate)]
            )  # fmt: skip
            qc["expected_rate_hz"] = EXPECTED_RATE
            qc_file.write_text(json.dumps(qc, indent=1))
            continue
        if decimated is None:
            continue
        counts += decimated.copy()
        boundaries = epoch_boundaries(inv, decimated)
        if boundaries:
            # Each part gets the response of its own epoch (obspy takes the
            # response valid at the start time of the trace it is given).
            qc.setdefault("response_epochs", []).extend(
                [decimated.stats.channel, str(t)] for t in boundaries
            )
            qc_file.write_text(json.dumps(qc, indent=1))
        for out in split_at(decimated, boundaries):
            if out.stats.npts * out.stats.delta < MIN_SEGMENT_S:
                continue
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    out.remove_response(
                        inventory=inv,
                        output="DEF" if out.stats.channel == P else "DISP",
                        pre_filt=PRE_FILT,
                        water_level=None,
                        taper=True,
                        taper_fraction=min(0.05, 600 / (out.stats.npts * GRID)),
                    )
            except ValueError as exc:
                if "No matching response" not in str(exc):
                    raise
                # Data outside every StationXML epoch: lose this part, not the chunk.
                qc.setdefault("no_response", []).append(
                    [
                        out.stats.channel,
                        str(out.stats.starttime),
                        str(out.stats.endtime),
                    ]
                )
                qc_file.write_text(json.dumps(qc, indent=1))
                continue
            # Drop the response/FIR edge zones of every segment.
            out.trim(out.stats.starttime + PAD_S / 2, out.stats.endtime - PAD_S / 2)
            phys += out
    written = []
    for kind, stream in (("counts", counts), ("phys", phys)):
        stream.merge(method=0)
        for day in chunk_days:
            part = stream.slice(day, day + DAY - GRID / 2, nearest_sample=False).copy()
            part = part.split()
            part.traces = [tr for tr in part if tr.stats.npts > 0]
            if not part:
                continue
            for tr in part:
                tr.data = tr.data.astype(np.float32)
            out = WORK / station / kind / f"{station}.{day.year}.{day.julday:03d}.mseed"
            out.parent.mkdir(parents=True, exist_ok=True)
            part.write(str(out), format="MSEED", encoding="FLOAT32")
            written.append(str(out))
    if not written and qc.get("wrong_rate"):
        return station, written, "no data at the expected rate"
    return station, written, "ok"


def main():
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stations", nargs="+")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    tasks = []
    empty = []
    for station in args.stations:
        days = station_days(station)
        if not days:
            empty.append(station)
        todo = [
            d
            for d in days
            if not (
                WORK / station / "phys" / f"{station}.{d.year}.{d.julday:03d}.mseed"
            ).exists()
        ]
        runs, run = [], []
        for day in todo:
            if run and day - run[-1] > 1.5 * DAY:
                runs.append(run)
                run = []
            run.append(day)
        runs += [run] if run else []
        for run in runs:
            for i in range(0, len(run), CHUNK_DAYS):
                tasks.append((station, run[i : i + CHUNK_DAYS]))
        print(f"{station}: {len(days)} days, {len(todo)} to do", flush=True)
    if empty:
        # Before 3.2.1 such a station was passed over and the run ended with exit 0.
        folder = SDS / "<year>" / NETWORK / "<station>" / f"{Z}.D"
        raise SystemExit(
            f"Stage A: no day file for {', '.join(empty)} in {folder} "
            f"(years {', '.join(str(y) for y in YEARS)}). Check config.SDS (or COMPY_SDS), "
            "NETWORK, YEARS and the station codes."
        )
    done = 0
    failed = 0
    with ProcessPoolExecutor(args.workers) as pool:
        futures = [pool.submit(process_chunk, *t) for t in tasks]
        for future in as_completed(futures):
            done += 1
            try:
                station, written, status = future.result()
                if status != "ok":
                    failed += 1
                print(
                    f"[{done}/{len(tasks)}] {station} {status} {len(written)} files",
                    flush=True,
                )
            except Exception as exc:
                failed += 1
                print(f"[{done}/{len(tasks)}] FAILED {exc!r}", flush=True)
    if failed:
        raise SystemExit(f"Stage A failed for {failed} preprocessing chunks")


if __name__ == "__main__":
    main()
