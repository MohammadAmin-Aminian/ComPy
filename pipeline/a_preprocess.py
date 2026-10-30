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
"""

import argparse
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

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

GRID = 1 / OUT_RATE
CHUNK_DAYS = 6
PAD_S = 1800.0
MIN_SEGMENT_S = 2 * PAD_S
PRE_FILT = (0.0008, 0.0015, 0.45, 0.6)
DAY = 86400.0
CLIP_COUNTS = 0.99 * 2**23  # 24-bit digitizer full scale
DEAD_STD_COUNTS = 5.0
QC_BLOCK_S = 600.0
QC_MARGIN_S = 300.0


def day_files(station, channel, day):
    path = SDS / str(day.year) / NETWORK / station / f"{channel}.D"
    return (
        path
        / f"{NETWORK}.{station}.{SDS_LOCATION}.{channel}.D.{day.year}.{day.julday:03d}"
    )


def station_days(station):
    days = set()
    for year in YEARS:
        folder = SDS / str(year) / NETWORK / station / f"{Z}.D"
        if folder.is_dir():
            for f in folder.iterdir():
                if not f.name.startswith("."):
                    days.add((year, int(f.name[-3:])))
    return [UTCDateTime(year=y, julday=j) for y, j in sorted(days)]


def grid_ceil(t):
    return UTCDateTime(np.ceil(t.timestamp / GRID) * GRID)


def join_day_files(traces, tolerance=1.5):
    """Join consecutive day-file traces of one channel into gap-free segments.

    Day files overlap by ~1 min and their sample grids slide by ~2.4 ms/day
    (clock correction). When the offset wraps past a whole sample, obspy's
    merge sees a one-sample conflict and masks the overlap. Here overlapping
    samples are dropped and traces whose junction is within ``tolerance``
    samples are concatenated (timing slip < 1 sample = 16 ms, removed again
    when snapping to the 0.8 s output grid); larger gaps start a new segment.
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
        if abs(tr.stats.starttime - expected) <= tolerance * delta:
            prev.data = np.concatenate([prev.data, tr.data])
        else:
            segments.append(tr)
    return segments


def bad_spans(trace):
    """Clipped or flat QC_BLOCK_S blocks of one raw trace, widened by QC_MARGIN_S."""
    width = int(round(QC_BLOCK_S * trace.stats.sampling_rate))
    n_blocks = -(-trace.stats.npts // width)
    data = np.asarray(trace.data, dtype=np.float64)
    spans = []
    for i in range(n_blocks):
        block = data[i * width : (i + 1) * width]
        if block.size < 2:
            continue
        clipped = np.abs(block).max() >= CLIP_COUNTS
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
        if start > cursor:
            pieces.append(trace.slice(cursor, start - trace.stats.delta).copy())
        cursor = max(cursor, stop)
    if cursor < end:
        pieces.append(trace.slice(cursor, end).copy())
    return [p for p in pieces if p.stats.npts > 1]


def decimate_segment(trace):
    """Decimate one gap-free raw trace on the 0.8 s grid; returns None if too short."""
    from tiskitpy import Decimator

    start = grid_ceil(trace.stats.starttime)
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
    qc_file.write_text(json.dumps({
        "days": [str(d) for d in chunk_days],
        "cut": [[str(a), str(b), sorted(r)] for a, b, r in bad],
        "cut_hours": sum(b - a for a, b, _ in bad) / 3600,
        "rule": f"|counts| >= {CLIP_COUNTS:.0f} or std < {DEAD_STD_COUNTS} in {QC_BLOCK_S:.0f}-s blocks, "
                f"+/- {QC_MARGIN_S:.0f} s",
    }, indent=1))  # fmt: skip
    counts, phys = Stream(), Stream()
    for trace in segments:
        out = decimate_segment(trace)
        if out is None:
            continue
        counts += out.copy()
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
    return station, written, "ok"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stations", nargs="+")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    tasks = []
    for station in args.stations:
        days = station_days(station)
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
