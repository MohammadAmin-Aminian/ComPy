"""Calculate response-corrected PPSDs and update the master inventory exports."""

import argparse
import hashlib
import json
import math
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import obspy
from obspy import UTCDateTime, read_inventory
from obspy.core import Stats
from obspy.signal import PPSD

from compy_pressure.archive import discover_waveforms, read_window
from compy_pressure.common import PROJECT_DIRECTORY
from compy_pressure.inventory import (
    EPOCH_COLUMNS,
    OUTPUT_DIRECTORY,
    load_master_dataset,
    update_processing_results,
)

RESULT_DIRECTORY = PROJECT_DIRECTORY / "results" / "TEST"


def is_pressure_unit(units):
    """Recognize the same declared pressure units in processing and selection."""
    return str(units).upper().replace(" ", "") in {"PA", "PASCAL", "PASCALS"}


def response_handling(channel):
    """Use declared physical units, never channel letters, to select correction."""
    response = channel.response
    if (
        not response
        or not response.response_stages
        or not response.instrument_sensitivity
    ):
        raise ValueError("Full response stages and instrument sensitivity are required")
    sensitivity = response.instrument_sensitivity
    # Some supplied pressure responses encode polarity with negative gains.
    # Keep that sign for complex correction; PSD correction uses |R| squared.
    if not np.isfinite(sensitivity.value) or sensitivity.value == 0:
        raise ValueError("Instrument sensitivity must be finite and nonzero")
    if sensitivity.output_units.upper().strip() not in {"COUNTS", "COUNT"}:
        raise ValueError("Expected raw digital counts at the response output")
    units = sensitivity.input_units.upper().replace(" ", "")
    if is_pressure_unit(units):
        return "hydrophone", "Pa^2/Hz"
    if units in {"M", "M/S", "M/S**2", "M/S^2", "M/S/S"}:
        return None, "m^2/s^4/Hz"
    raise ValueError(f"Unsupported response input units: {sensitivity.input_units}")


def select_epoch(row):
    """Select exactly the source channel epoch represented by a master row."""
    network, station, location, code = row["Seed_id"].split(".")
    inventory = read_inventory(row["inventory_file"], format="STATIONXML").select(
        network=network,
        station=station,
        location=location,
        channel=code,
    )
    channels = [
        channel
        for net in inventory
        for sta in net
        for channel in sta
        if str(channel.start_date or "") == row["channel_start_time"]
        and str(channel.end_date or "") == row["channel_end_time"]
    ]
    # Some StationXML files repeat a station block verbatim; identical copies
    # describe one epoch. Differing copies remain ambiguous and are refused.
    if len(channels) > 1 and all(other == channels[0] for other in channels[1:]):
        channels = channels[:1]
    if len(channels) != 1:
        raise ValueError(f"Expected one exact StationXML epoch, found {len(channels)}")
    channel = channels[0]
    for net in inventory:
        for sta in net:
            sta.channels = [candidate for candidate in sta if candidate is channel]
    return inventory, channel


def atomic_json(path, content):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(content, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def save_checkpoint(ppsd, path, manifest):
    """Commit spectra before progress, so an interrupted run can safely retry."""
    temporary = path.with_name("ppsd.tmp.npz")
    ppsd.save_npz(str(temporary))
    temporary.replace(path)
    atomic_json(path.with_name("processing.json"), manifest)


def process_epoch(
    row,
    records,
    results_directory,
    *,
    start=None,
    end=None,
    ppsd_length=10800.0,
    overlap=0.5,
    minimum_frequency=0.001,
    chunk_seconds=86400.0,
    rebuild=False,
    cancel_event=None,
):
    """Process one epoch in bounded, overlapping reads and resume checkpoints.

    The time interval is half-open; only complete PSD windows within it are used.
    FFT support is checked separately from the instrument's physical usable band.
    No instrument cutoff or calibration pass/fail is inferred by this function.
    """
    if not (ppsd_length > 0 and 0 <= overlap < 1 and minimum_frequency > 0):
        raise ValueError("Require positive length/frequency and 0 <= overlap < 1")
    step = ppsd_length * (1 - overlap)
    if chunk_seconds < step:
        raise ValueError("Chunk duration must be at least one PSD step")
    chunk_seconds = math.floor(chunk_seconds / step) * step
    inventory, channel = select_epoch(row)
    handling, psd_units = response_handling(channel)
    rate = float(channel.sample_rate)
    if not math.isclose(rate, float(row["sample_rate_hz"]), rel_tol=1e-9):
        raise ValueError("Master sample rate differs from StationXML")
    if minimum_frequency >= rate / 2:
        raise ValueError("Minimum frequency must be below Nyquist")
    records = [record for record in records if record.seed_id == row["Seed_id"]]
    epoch_start = channel.start_date
    epoch_end = channel.end_date
    records = [
        record
        for record in records
        if (epoch_start is None or record.end >= epoch_start)
        and (epoch_end is None or record.start < epoch_end)
    ]
    if not records:
        return {
            "ppsd_status": "no_data",
            "ppsd_error": "No waveform data in this epoch",
        }
    if any(not math.isclose(record.sample_rate, rate) for record in records):
        raise ValueError("Waveform sample rate differs from the StationXML epoch")
    origin = min(record.start for record in records)
    if epoch_start is not None:
        origin = max(origin, epoch_start)
    available_end = max(record.end for record in records) + 1 / rate
    if epoch_end is not None:
        available_end = min(available_end, epoch_end)
    requested_start = max(origin, UTCDateTime(start)) if start else origin
    requested_end = min(available_end, UTCDateTime(end)) if end else available_end
    if requested_end <= requested_start:
        return {"ppsd_status": "no_data", "ppsd_error": "Empty requested time interval"}
    stats = Stats(
        dict(
            zip(
                ["network", "station", "location", "channel"],
                row["Seed_id"].split("."),
                strict=True,
            )
        )
    )
    stats.sampling_rate = rate
    ppsd = PPSD(
        stats,
        inventory,
        ppsd_length=ppsd_length,
        overlap=overlap,
        skip_on_gaps=True,
        special_handling=handling,
        db_bins=(-250.0, 100.0, 1.0),
        period_limits=(2 / rate, 1 / minimum_frequency),
    )
    # The full low-frequency smoothing bin must be supported by actual FFT bins.
    if ppsd.psd_periods[-1] < max(ppsd.period_bin_right_edges):
        raise ValueError(
            "PSD windows are too short for the requested low frequency; "
            "increase --ppsd-length"
        )
    response, frequencies = channel.response.get_evalresp_response(
        t_samp=1 / rate,
        nfft=ppsd.nfft,
        output="DEF" if handling else "VEL",
    )
    if not np.isfinite(response[1:]).all() or np.any(np.abs(response[1:]) == 0):
        raise ValueError("Instrument response is nonfinite or zero in the FFT band")
    inventory_path = Path(row["inventory_file"]).expanduser().resolve()
    inventory_hash = hashlib.sha256(inventory_path.read_bytes()).hexdigest()
    source_files = [
        {
            "path": str(path),
            "size": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns,
        }
        for path in sorted({record.path for record in records})
    ]
    config = {
        "seed_id": row["Seed_id"],
        "epoch_start": str(epoch_start or ""),
        "epoch_end": str(epoch_end or ""),
        "sample_rate_hz": rate,
        "ppsd_length_seconds": ppsd_length,
        "overlap": overlap,
        "minimum_frequency_hz": minimum_frequency,
        "chunk_seconds": chunk_seconds,
        "special_handling": handling,
        "psd_units": psd_units,
        "db_bins": [-250.0, 100.0, 1.0],
        "skip_on_gaps": True,
        "origin": str(origin),
        "inventory_sha256": inventory_hash,
        "obspy_version": obspy.__version__,
        "chunk_strategy": "sample_aligned_v2",
    }
    epoch_tag = epoch_start.strftime("%Y%m%dT%H%M%S") if epoch_start else "open"
    network, station, _, _ = row["Seed_id"].split(".")
    folder = (
        Path(results_directory).expanduser().resolve()
        / network
        / station
        / (f"{row['Seed_id']}__{epoch_tag}__{inventory_hash[:12]}")
    )
    folder.mkdir(parents=True, exist_ok=True)
    npz_path = folder / "ppsd.npz"
    manifest_path = folder / "processing.json"
    manifest = {
        "config": config,
        "source_files": source_files,
        "source_inventory": str(inventory_path),
        "completed_chunks": {},
        "errors": {},
        "warning_counts": {},
    }
    if manifest_path.exists() and not rebuild:
        previous = json.loads(manifest_path.read_text())
        if previous["config"] != config or previous["source_files"] != source_files:
            raise ValueError(
                "Saved settings/sources differ; use --rebuild to recalculate"
            )
        if not npz_path.exists():
            raise ValueError("Checkpoint NPZ is missing; use --rebuild")
        manifest = previous
        # Initialize with the real stats before loading spectra. ObsPy 1.5.1
        # load_npz(metadata=...) preloads responses before restoring saved stats.
        ppsd.add_npz(str(npz_path))
    elif npz_path.exists() and not rebuild:
        raise ValueError("NPZ has no provenance; use --rebuild")
    inventory.write(str(folder / "inventory.xml"), format="STATIONXML")
    # Resume from the next accepted PSD start. Small timestamp phase changes at
    # daily file boundaries must not cause ObsPy to reject a whole PSD window as
    # excessively overlapping. Keep original sample timestamps; move the next
    # read forward instead of relabeling samples or interpolating them.
    first = max(0, math.ceil((requested_start - origin) / step - 1e-10))
    cursor = origin + first * step
    attempted = 0
    interrupted = False
    while cursor + ppsd_length <= requested_end:
        if cancel_event is not None and cancel_event.is_set():
            interrupted = True
            break
        window_start = cursor
        key = str(window_start)
        if key in manifest["completed_chunks"]:
            cursor = UTCDateTime(manifest["completed_chunks"][key]["next_start"])
            continue
        # Two extra samples allow the last window despite fractional sample
        # phases between files; the requested/epoch end remains a hard limit.
        window_end = min(
            cursor + chunk_seconds + ppsd_length - step + 2 / rate,
            requested_end,
        )
        next_start = cursor + chunk_seconds
        try:
            stream = read_window(records, window_start, window_end)
            # PPSD handles response removal once; input stays in raw counts.
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                ppsd.add(stream)
            response_errors = []
            for warning in caught:
                message = str(warning.message)
                if message.startswith("Already covered time spans"):
                    continue
                category = message.split("\n")[0].split("Skipping trace:")[0]
                manifest["warning_counts"][category] = (
                    manifest["warning_counts"].get(category, 0) + 1
                )
                if "Error getting response" in message:
                    response_errors.append(message)
            if response_errors:
                raise ValueError("; ".join(response_errors))
            block_times = [
                time
                for time in ppsd.times_processed
                if window_start <= time and time + ppsd_length <= window_end
            ]
            if block_times:
                if window_end < cursor + chunk_seconds + ppsd_length - step:
                    # A pilot ended partway through this read. Its next PSD
                    # window must remain available when that pilot is extended.
                    next_start = block_times[-1] + step
                else:
                    next_start = max(next_start, block_times[-1] + step)
            manifest["completed_chunks"][key] = {
                "read_end": str(window_end),
                "next_start": str(next_start),
            }
            manifest["errors"].pop(key, None)
        except Exception as exc:
            manifest["errors"][key] = f"{type(exc).__name__}: {exc}"
        attempted += 1
        if attempted % 10 == 0:
            save_checkpoint(ppsd, npz_path, manifest)
            print(
                f"{row['Seed_id']}: through {window_end}, "
                f"{len(ppsd.times_processed)} PSD segments",
                flush=True,
            )
        cursor = next_start
    segments = len(ppsd.times_processed)
    status = "complete" if segments and not manifest["errors"] else "partial"
    if not segments:
        status = "failed" if manifest["errors"] else "insufficient_data"
    if interrupted:
        status = "partial"
    summary = {
        "ppsd_status": status,
        "ppsd_npz_file": str(npz_path),
        "ppsd_provenance_file": str(manifest_path),
        "ppsd_units": psd_units,
        "ppsd_segments": segments,
        "ppsd_source_files": len(source_files),
        "ppsd_length_seconds": ppsd_length,
        "ppsd_overlap": overlap,
        "ppsd_frequency_min_hz": float(1 / max(ppsd.period_bin_centers)),
        "ppsd_frequency_max_hz": float(1 / min(ppsd.period_bin_centers)),
        "ppsd_fft_min_frequency_hz": float(min(frequencies[1:])),
        "ppsd_nyquist_hz": rate / 2,
        "ppsd_histogram_min_db": float(ppsd.db_bin_edges[0]),
        "ppsd_histogram_max_db": float(ppsd.db_bin_edges[-1]),
        "ppsd_histogram_out_of_range_values": sum(
            int(
                np.count_nonzero(
                    (values < ppsd.db_bin_edges[0]) | (values > ppsd.db_bin_edges[-1])
                )
            )
            for values in ppsd.psd_values
        ),
        "ppsd_start_time": str(ppsd.times_processed[0]) if segments else None,
        "ppsd_end_time": str(ppsd.times_processed[-1] + ppsd_length)
        if segments
        else None,
        "ppsd_requested_start_time": str(requested_start),
        "ppsd_requested_end_time": str(requested_end),
        "ppsd_updated_at": datetime.now(timezone.utc).isoformat(),
        "ppsd_error": "Processing interrupted; rerun to resume"
        if interrupted
        else "; ".join(manifest["errors"].values())[:1000] or None,
    }
    manifest["summary"] = summary
    save_checkpoint(ppsd, npz_path, manifest)
    print(
        f"{row['Seed_id']}: {status}, {segments} segments saved to {npz_path}",
        flush=True,
    )
    return summary


def process_task(row, records, results_directory, options):
    """Isolate channel failures and return a master update to the parent process."""
    key = {name: row[name] for name in EPOCH_COLUMNS + ["inventory_file"]}
    try:
        summary = process_epoch(row, records, results_directory, **options)
    except Exception as exc:
        print(f"{row['Seed_id']}: failed: {exc}", flush=True)
        summary = {"ppsd_status": "failed", "ppsd_error": str(exc)}
    return key | summary


def run_processing(
    data_directory,
    master_directory,
    results_directory,
    *,
    workers=1,
    **options,
):
    """Process matching epochs; only the parent process writes the master files."""
    if workers < 1:
        raise ValueError("Worker count must be at least one")
    master = load_master_dataset(master_directory)
    seed_ids = options.pop("seed_ids", None)
    if seed_ids is None:
        seed_ids = set(master["Seed_id"])
    records = discover_waveforms(data_directory, seed_ids=set(seed_ids))
    available = {record.seed_id for record in records}
    selected = master[master.Seed_id.isin(available)].drop_duplicates(
        EPOCH_COLUMNS + ["inventory_file"],
    )
    rows = selected.to_dict("records")
    if workers == 1:
        updates = [
            process_task(row, records, results_directory, options) for row in rows
        ]
    else:
        with ProcessPoolExecutor(max_workers=min(workers, len(rows))) as pool:
            futures = [
                pool.submit(process_task, row, records, results_directory, options)
                for row in rows
            ]
            updates = [future.result() for future in futures]
    update_processing_results(master_directory, updates)
    return updates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--master-dir", type=Path, default=OUTPUT_DIRECTORY)
    parser.add_argument("--results-dir", type=Path, default=RESULT_DIRECTORY)
    parser.add_argument("--seed-id", action="append", dest="seed_ids")
    parser.add_argument("--start", help="Inclusive UTC start; e.g. 2013-01-01")
    parser.add_argument("--end", help="Exclusive UTC end; e.g. 2013-01-03")
    parser.add_argument("--ppsd-length", type=float, default=10800.0)
    parser.add_argument("--overlap", type=float, default=0.5)
    parser.add_argument("--minimum-frequency", type=float, default=0.001)
    parser.add_argument("--chunk-seconds", type=float, default=86400.0)
    parser.add_argument(
        "--workers", type=int, default=1, help="Independent channel workers"
    )
    parser.add_argument(
        "--rebuild", action="store_true", help="Recalculate saved results"
    )
    args = vars(parser.parse_args())
    data = args.pop("data_dir")
    master = args.pop("master_dir")
    results = args.pop("results_dir")
    updates = run_processing(data, master, results, **args)
    if any(update["ppsd_status"] in {"failed", "partial"} for update in updates):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
