"""Stage A: decimate raw YV OBS data to 1.25 Hz and remove instrument responses.

For each station and day this writes two MiniSEED files under work/<STA>/:
  counts/<STA>.<YYYY.JJJ>.mseed   decimated raw counts (BHZ, BH1, BH2, BDH)
  phys/<STA>.<YYYY.JJJ>.mseed     BHZ/BH1/BH2 displacement [m], BDH pressure [Pa]

Choices (documented in the report):
  * Decimation 62.5 -> 1.25 Hz with TiSKitPy SAC FIR filters, factors 5, 5, 2,
    as in Aminian (2024, thesis section 3.1.1). Channel codes are restored after
    decimation so the StationXML responses still match.
  * Each raw segment is trimmed to start on the 0.8 s output grid (residual
    offset < one raw sample, 16 ms, identical on all channels of the digitizer).
  * Responses are removed after decimation with water_level=None and an
    explicit pre-filter. With obspy's default 60 dB water level, a displacement
    response removed at 62.5 Hz is clipped below ~0.03 Hz (tested: 22 dB error).
  * Chunks of CHUNK_DAYS days are padded by PAD_S seconds on each side so FIR
    and response edge effects fall outside the kept day.
"""

import argparse
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from obspy import Stream, UTCDateTime, read, read_inventory

SDS = Path("/Volumes/DATA_bruitfm/SDS")
INVENTORY = Path("/Volumes/DATA_bruitfm/INVENTORIES/YV.xml")
WORK = Path(__file__).resolve().parents[1] / "work"
CHANNELS = ("BHZ", "BH1", "BH2", "BDH")
FACTORS = (5, 5, 2)
OUT_RATE = 1.25
GRID = 1 / OUT_RATE
CHUNK_DAYS = 6
PAD_S = 1800.0
MIN_SEGMENT_S = 2 * PAD_S
PRE_FILT = (0.0008, 0.0015, 0.45, 0.6)
DAY = 86400.0


def day_files(station, channel, day):
    path = SDS / str(day.year) / "YV" / station / f"{channel}.D"
    return path / f"YV.{station}.00.{channel}.D.{day.year}.{day.julday:03d}"


def station_days(station):
    days = set()
    for year in (2012, 2013):
        folder = SDS / str(year) / "YV" / station / "BHZ.D"
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
    out.stats.starttime = out.stats.starttime - offset  # snap (<16 ms) to grid
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
                output="DEF" if out.stats.channel == "BDH" else "DISP",
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
    with ProcessPoolExecutor(args.workers) as pool:
        futures = [pool.submit(process_chunk, *t) for t in tasks]
        for future in as_completed(futures):
            done += 1
            try:
                station, written, status = future.result()
                print(
                    f"[{done}/{len(tasks)}] {station} {status} {len(written)} files",
                    flush=True,
                )
            except Exception as exc:
                print(f"[{done}/{len(tasks)}] FAILED {exc!r}", flush=True)


if __name__ == "__main__":
    main()
