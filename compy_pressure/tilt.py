"""Tilt-corrected vertical PPSDs with ComPy (TiSKitPy CleanRotator + DataCleaner).

For one station location, a vertical and two horizontal channels of one band
are read day by day, low-passed and resampled to 1 Hz, corrected to
acceleration in chunks of several days (tapered ends discarded) and cleaned of
tilt and horizontal-coherent noise by ``compy.Rotate`` in WINDOW_HOURS windows.
A window whose rotation fails is dropped and counted. Two PPSDs with identical
windows are saved: the 1 Hz vertical before and after cleaning, in ``tilt/``
beside the vertical channel's saved PPSD. Raw data and saved PPSDs are only
read.

Run after the seismic ``ppsd`` stage, which saves the response snapshots:

    python -m compy_pressure.tilt --master-dir M --data-dir SDS --results-dir R --workers max
"""

import argparse
import csv
import json
import logging
import math
import os
import traceback
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from multiprocessing import get_context
from pathlib import Path

import numpy as np
from obspy import Stream, UTCDateTime, read, read_inventory
from obspy.signal import PPSD

from compy_pressure.inventory import load_master_dataset

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

TARGET_RATE = 1.0
WINDOW_HOURS = 6
CHUNK_DAYS = 7
MARGIN_SECONDS = 3600
PRE_FILT = (0.0005, 0.001, 0.4, 0.45)
PPSD_LENGTH = 10800
PPSD_OVERLAP = 0.5
PERIOD_LIMITS = (2.5, 1000.0)
DB_BINS = (-250.0, -40.0, 1.0)
BAND_PREFERENCE = "LBHESM"
BANDS_HZ = ((0.002, 0.01), (0.01, 0.03), (0.03, 0.1), (0.1, 0.3))
HORIZONTAL = {"1": "1", "2": "2", "N": "1", "E": "2"}
# Some networks (1T) code the vertical "3"; its StationXML dip is +-90.
VERTICAL = "Z3"
THREAD_VARIABLES = (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def is_seismometer(code):
    return len(code) == 3 and code[1] == "H" and code[2] in VERTICAL + "12NE"


def choose_triplets(master):
    """One vertical with two horizontals per station location and Z epoch.

    Only channels with a completed PPSD (and so a saved response snapshot) are
    used. The band is chosen in BAND_PREFERENCE order: 1 Hz L channels first,
    since they are cheap to read and already at the target rate.
    """
    rows = master[master.ppsd_status.isin(["complete", "partial"])].copy()
    rows["code"] = rows.Seed_id.str.split(".").str[3]
    rows = rows[rows.code.map(is_seismometer)]
    rows["location"] = rows.Seed_id.str.rsplit(".", n=1).str[0]
    triplets = []
    for location, group in rows.groupby("location"):
        for band in BAND_PREFERENCE:
            members = group[group.code.str[0].eq(band)]
            verticals = members[members.code.str[2].isin(list(VERTICAL))]
            if verticals.empty:
                continue
            found = []
            for _, z in verticals.iterrows():
                start, end = z.channel_start_time, z.channel_end_time or "9999"
                overlap = members[
                    members.code.str[2].isin(list(HORIZONTAL))
                    & (members.channel_start_time < end)
                    & (
                        (members.channel_end_time.eq(""))
                        | (members.channel_end_time > start)
                    )
                ]
                one = overlap[overlap.code.str[2].map(HORIZONTAL).eq("1")]
                two = overlap[overlap.code.str[2].map(HORIZONTAL).eq("2")]
                if len(one) and len(two):
                    found.append((z, one.iloc[0], two.iloc[0]))
            if found:
                triplets.extend(found)
                break
    return triplets


def day_key(path):
    parts = Path(path).name.split(".")
    return f"{int(parts[5]):04d}.{int(parts[6]):03d}"


def to_target_rate(stream):
    """Detrend, low-pass and resample each contiguous trace to 1 Hz on whole seconds."""
    out = Stream()
    for tr in stream.split():
        rate = tr.stats.sampling_rate
        if tr.stats.npts < 600 * rate:
            continue
        tr.data = tr.data.astype(np.float64)
        tr.detrend("linear")
        if abs(rate - TARGET_RATE) > 1e-9:
            tr.filter("lowpass", freq=0.4 * TARGET_RATE, corners=4, zerophase=True)
            start = UTCDateTime(math.ceil(tr.stats.starttime.timestamp))
            tr.interpolate(TARGET_RATE, method="lanczos", a=12, starttime=start)
        else:
            offset = tr.stats.starttime.timestamp % 1
            if offset > 1e-3:
                continue
            tr.stats.starttime = UTCDateTime(round(tr.stats.starttime.timestamp))
        out += tr
    return out


def read_day(paths, station, code):
    stream = Stream()
    for path in paths:
        if Path(path).stat().st_size:
            stream += read(str(path)).select(station=station, channel=code)
    if not stream:
        return stream
    stream.merge(method=1, fill_value=None)
    return to_target_rate(stream)


def corrected_segments(stream, inventory):
    """Response-corrected (acceleration) contiguous traces with tapered ends removed."""
    out = Stream()
    stream.merge(method=1, fill_value=None)
    for tr in stream.split():
        duration = tr.stats.endtime - tr.stats.starttime
        if duration < WINDOW_HOURS * 3600 + 2 * MARGIN_SECONDS:
            continue
        tr.detrend("linear")
        fraction = min(0.5, 1.8 * MARGIN_SECONDS / duration)
        tr.remove_response(
            inventory=inventory,
            output="ACC",
            pre_filt=PRE_FILT,
            taper=True,
            taper_fraction=fraction,
        )
        tr.trim(tr.stats.starttime + MARGIN_SECONDS, tr.stats.endtime - MARGIN_SECONDS)
        out += tr
    return out


def common_intervals(z, one, two):
    """Slices of Z, 1 and 2 covering the same samples, renamed to Z/1/2."""
    triples = []
    for tz in z:
        for t1 in one:
            for t2 in two:
                start = max(tz.stats.starttime, t1.stats.starttime, t2.stats.starttime)
                end = min(tz.stats.endtime, t1.stats.endtime, t2.stats.endtime)
                if end - start < WINDOW_HOURS * 3600:
                    continue
                parts = [t.slice(start, end).copy() for t in (tz, t1, t2)]
                n = min(len(p) for p in parts)
                for p, suffix in zip(parts, ("Z", "1", "2"), strict=True):
                    p.data = p.data[:n]
                    p.stats.channel = p.stats.channel[:2] + suffix
                triples.append(Stream(parts))
    return triples


def clean_windows(stream, rotate):
    """Apply compy.Rotate window by window; failed windows are dropped."""
    from compy_streams import split_stream

    raw, clean, rows = Stream(), Stream(), []
    for window in split_stream(stream, WINDOW_HOURS * 3600):
        start = window[0].stats.starttime
        try:
            rotated, azimuth, angle, ratio = rotate(window, WINDOW_HOURS, plot=False)
            z = rotated.select(component="Z")[0]
            if not np.all(np.isfinite(z.data)):
                raise ValueError("non-finite cleaned samples")
        except Exception as exc:
            rows.append(
                {"start": str(start), "status": "failed", "error": str(exc)[:200]}
            )
            continue
        raw += window.select(component="Z")[0].copy()
        clean += z
        rows.append(
            {
                "start": str(start),
                "status": "cleaned",
                "azimuth_deg": float(azimuth[0]),
                "tilt_deg": float(angle[0]),
                "variance_after_over_before": float(ratio[0]),
            }
        )
    return raw, clean, rows


def new_ppsd(trace):
    return PPSD(
        trace.stats,
        # Data are already acceleration: "ringlaser" handling only divides by
        # the sensitivity, so a unit, pole-free response leaves them unchanged.
        metadata={"poles": [], "zeros": [], "gain": 1.0, "sensitivity": 1.0},
        special_handling="ringlaser",
        ppsd_length=PPSD_LENGTH,
        overlap=PPSD_OVERLAP,
        period_limits=PERIOD_LIMITS,
        db_bins=DB_BINS,
    )


def band_differences(raw, clean):
    """Median over common windows of cleaned minus raw PSD (dB) in BANDS_HZ."""
    f = 1 / np.asarray(raw.period_bin_centers)
    a, b = np.asarray(raw.psd_values), np.asarray(clean.psd_values)
    n = min(len(a), len(b))
    result = {}
    for low, high in BANDS_HZ:
        keep = (f >= low) & (f < high)
        if n and keep.any():
            result[f"{low:g}-{high:g}"] = float(
                np.median(b[:n][:, keep] - a[:n][:, keep])
            )
        else:
            result[f"{low:g}-{high:g}"] = None
    return result


def triplet_job(z, one, two, files):
    """Process one station triplet; return a summary dictionary."""
    logging.disable(logging.INFO)
    import compy

    output = Path(z["ppsd_npz_file"]).parent / "tilt"
    summary_path = output / "summary.json"
    seeds = [z["Seed_id"], one["Seed_id"], two["Seed_id"]]
    if summary_path.exists():
        previous = json.loads(summary_path.read_text())
        if previous.get("status") == "complete" and previous.get("channels") == seeds:
            return previous
    output.mkdir(exist_ok=True)
    inventory = read_inventory(str(Path(z["ppsd_npz_file"]).parent / "inventory.xml"))
    for row in (one, two):
        inventory += read_inventory(
            str(Path(row["ppsd_npz_file"]).parent / "inventory.xml")
        )
    start = UTCDateTime(max(r["channel_start_time"] for r in (z, one, two)))
    ends = [r["channel_end_time"] for r in (z, one, two) if r["channel_end_time"]]
    end = UTCDateTime(min(ends)) if ends else UTCDateTime(2100, 1, 1)
    by_day = [defaultdict(list) for _ in seeds]
    for index, seed in enumerate(seeds):
        for path in files.get(seed, []):
            by_day[index][day_key(path)].append(path)
    days = sorted(set.intersection(*(set(d) for d in by_day)))
    days = [
        d
        for d in days
        if start - 86400
        <= UTCDateTime(f"{d[:4]}-01-01") + (int(d[5:]) - 1) * 86400
        < end
    ]
    ppsds, rows, used_days = {}, [], 0
    for begin in range(0, len(days), CHUNK_DAYS):
        chunk = days[begin : begin + CHUNK_DAYS]
        streams = []
        for index, seed in enumerate(seeds):
            _, station, _, code = seed.split(".")
            stream = Stream()
            for day in chunk:
                stream += read_day(by_day[index][day], station, code)
            streams.append(
                corrected_segments(stream, inventory) if stream else Stream()
            )
        if not all(streams):
            continue
        used_days += len(chunk)
        for triple in common_intervals(*streams):
            raw, clean, window_rows = clean_windows(triple, compy.Rotate)
            rows.extend(window_rows)
            for kind, traces in (("raw", raw), ("clean", clean)):
                traces.merge(method=1, fill_value=None)
                for tr in traces.split():
                    if kind not in ppsds:
                        ppsds[kind] = new_ppsd(tr)
                    ppsds[kind].add(tr)
    cleaned = sum(r["status"] == "cleaned" for r in rows)
    summary = {
        "status": "complete" if ppsds else "no_common_data",
        "channels": seeds,
        "method": (
            f"1 Hz, response removed to acceleration (pre_filt {PRE_FILT}), "
            f"compy.Rotate in {WINDOW_HOURS} h windows; failed windows dropped"
        ),
        "compy_version": getattr(compy, "__version__", ""),
        "days_with_all_three_channels": len(days),
        "days_processed": used_days,
        "windows_cleaned": cleaned,
        "windows_failed": len(rows) - cleaned,
        "updated_utc": datetime.now(UTC).isoformat(),
    }
    if ppsds:
        raw_path = output / "ppsd_vertical_1hz.npz"
        clean_path = output / "ppsd_vertical_tilt_corrected.npz"
        ppsds["raw"].save_npz(str(raw_path))
        ppsds["clean"].save_npz(str(clean_path))
        angles = [r["tilt_deg"] for r in rows if r["status"] == "cleaned"]
        summary |= {
            "psd_windows": len(ppsds["clean"].times_processed),
            "median_tilt_deg": float(np.median(angles)) if angles else None,
            "median_cleaned_minus_raw_db": band_differences(
                ppsds["raw"], ppsds["clean"]
            ),
            "files": {
                "raw_ppsd": str(raw_path),
                "tilt_corrected_ppsd": str(clean_path),
            },
        }
    with (output / "windows.csv").open("w", newline="") as stream:
        names = ["start", "status", "azimuth_deg", "tilt_deg"]
        names += ["variance_after_over_before", "error"]
        writer = csv.DictWriter(stream, fieldnames=names, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def safe_job(z, one, two, files):
    try:
        return triplet_job(z, one, two, files)
    except Exception as exc:
        return {
            "status": "failed",
            "channels": [z["Seed_id"], one["Seed_id"], two["Seed_id"]],
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc()[-2000:],
        }


def run_tilt(master_dir, data_dir, results_dir, *, networks=None, workers=4):
    """Run every triplet in parallel; write <results_dir>/_workflow/tilt.json."""
    from compy_pressure.network import apply_location_aliases, scan_files

    master = load_master_dataset(master_dir)
    if networks:
        master = master[master.Seed_id.str.split(".").str[0].isin(networks)]
    triplets = choose_triplets(master)
    scope = sorted({z["Seed_id"].split(".")[0] for z, _, _ in triplets})
    files = apply_location_aliases(scan_files(data_dir, scope), set(master.Seed_id))
    report_path = Path(results_dir) / "_workflow" / "tilt.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report = {"entries": {}}
    print(
        f"tilt: {len(triplets)} station triplets in {len(scope)} networks", flush=True
    )
    for name in THREAD_VARIABLES:
        os.environ.setdefault(name, "1")
    with ProcessPoolExecutor(
        max_workers=max(1, min(workers, len(triplets) or 1)),
        mp_context=get_context("spawn"),
    ) as pool:
        futures = {}
        for z, one, two in triplets:
            seeds = (z["Seed_id"], one["Seed_id"], two["Seed_id"])
            subset = {s: files.get(s, []) for s in seeds}
            key = f"{z['Seed_id']}|{z['channel_start_time']}"
            futures[
                pool.submit(safe_job, z.to_dict(), one.to_dict(), two.to_dict(), subset)
            ] = key
        for count, future in enumerate(as_completed(futures), 1):
            summary = future.result()
            report["entries"][futures[future]] = summary
            print(
                f"tilt {count}/{len(futures)} {futures[future]}: {summary['status']}"
                + (f" ({summary.get('error')})" if summary.get("error") else ""),
                flush=True,
            )
            if count % 10 == 0:
                report_path.write_text(json.dumps(report, indent=1) + "\n")
    counts = defaultdict(int)
    for entry in report["entries"].values():
        counts[entry["status"]] += 1
    report["status_counts"] = dict(counts)
    report_path.write_text(json.dumps(report, indent=1) + "\n")
    print(f"tilt: {dict(counts)}; details: {report_path}", flush=True)
    return report


def main():
    from compy_pressure.network import worker_count

    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--master-dir", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--network", action="append")
    parser.add_argument("--workers", type=worker_count, default=4)
    args = parser.parse_args()
    run_tilt(
        args.master_dir,
        args.data_dir,
        args.results_dir,
        networks=args.network,
        workers=args.workers,
    )


if __name__ == "__main__":
    main()
