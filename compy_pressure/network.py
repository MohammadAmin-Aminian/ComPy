"""Run the existing spectrum workflow for all networks found in a data folder."""

import argparse
import hashlib
import json
import os
import signal
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from fnmatch import fnmatchcase
from multiprocessing import get_context
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from compy_pressure.archive import alias_seed_id, index_waveform_files, waveform_seed_id
from compy_pressure.common import PROJECT_DIRECTORY
from compy_pressure.infragravity_test import run_infragravity_test
from compy_pressure.inventory import (
    EPOCH_COLUMNS,
    OUTPUT_DIRECTORY,
    RESULT_COLUMNS,
    XML_DIRECTORY,
    load_master_dataset,
    synchronize_inventory,
    update_infragravity_test_results,
    update_nearby_test_results,
    update_paired_results,
    update_pressure_test_results,
    update_processing_results,
)
from compy_pressure.nearby_test import (
    FINISHED as NEARBY_FINISHED,
)
from compy_pressure.nearby_test import (
    Criteria,
    add_criteria_arguments,
    criteria_from_args,
    ensure_master_idle,
    pool_fingerprint,
    pressure_rows,
    run_nearby_test,
)
from compy_pressure.paired_spectra import (
    FINISHED as PAIRED_FINISHED,
)
from compy_pressure.paired_spectra import (
    assert_workflows_finished,
    run_paired_spectra,
)
from compy_pressure.plot_ppsd import plot_station
from compy_pressure.ppsd import atomic_json, is_pressure_unit, process_task
from compy_pressure.pressure_test import run_pressure_test

DATA_DIRECTORY = Path.home() / "Desktop/Data_Bruit_FM/Data"
RESULT_DIRECTORY = PROJECT_DIRECTORY / "results/NETWORKS/TEST"
STAGES = (
    "inventory",
    "check",
    "ppsd",
    "plots",
    "pressure",
    "nearby",
    "infragravity",
    "paired",
)
KEYS = EPOCH_COLUMNS + ["inventory_file"]
SUCCESS = {"complete", "comparison_only", "no_comparable_bins", "screened"}
SUCCESS |= NEARBY_FINISHED | PAIRED_FINISHED
UPDATE_FUNCTIONS = {
    "ppsd": update_processing_results,
    "plots": update_processing_results,
    "pressure": update_pressure_test_results,
    "infragravity": update_infragravity_test_results,
    "nearby": update_nearby_test_results,
    "paired": update_paired_results,
}
PREFIXES = {
    "plots": "ppsd_plot",
    "pressure": "pressure_expected",
    "infragravity": "infragravity",
    "nearby": "nearby",
    "paired": "paired",
}
_CANCEL_EVENT = None
# Spawned workers inherit these before NumPy loads: one BLAS/FFT thread each.
THREAD_VARIABLES = (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def worker_count(value):
    """Parse --workers: a positive integer, or `max` for all cores but one."""
    if value == "max":
        return max(1, (os.cpu_count() or 2) - 1)
    count = int(value)
    if count < 1:
        raise argparse.ArgumentTypeError("Workers must be a positive integer or max")
    return count


def initialize_worker(cancel_event):
    """Let the parent handle Ctrl-C; workers finish their current bounded read."""
    global _CANCEL_EVENT
    _CANCEL_EVENT = cancel_event
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def _walk_files(root, networks=None):
    """Yield files in sorted order without one stat call per file.

    Directory entries carry their type, which matters on large spinning
    disks. In an SDS ``YEAR/NET/...`` tree, unselected network folders are
    skipped; filenames are still filtered, so this only saves time.
    """
    with os.scandir(root) as entries:
        entries = sorted(entries, key=lambda entry: entry.name)
    year = len(root.name) == 4 and root.name.isdigit()
    for entry in entries:
        if entry.is_dir(follow_symlinks=False):
            if networks and year and entry.name not in networks:
                continue
            yield from _walk_files(Path(entry.path), networks)
        elif entry.is_file():
            yield Path(entry.path)


def scan_files(data_directory, networks=None, stations=None, seed_ids=None):
    """Group daily filenames once; inspect trace headers later in channel jobs."""
    root = Path(data_directory).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Waveform directory does not exist: {root}")
    grouped = defaultdict(list)
    for path in sorted(_walk_files(root, networks)):
        seed_id = waveform_seed_id(path)
        if seed_id is None:
            continue
        parts = seed_id.split(".")
        if networks and parts[0] not in networks:
            continue
        if stations and parts[1] not in stations:
            continue
        if seed_ids and seed_id not in seed_ids:
            continue
        grouped[seed_id].append(path)
    if not grouped:
        raise ValueError("No matching SDS waveform files found for the selected scope")
    return dict(grouped)


def apply_location_aliases(files, inventory_ids):
    """Key aliased filenames by the inventory channel their headers must match.

    Only a filename identifier absent from the inventory is renamed, and only
    to an inventory identifier; headers are still checked in each job.
    """
    resolved = defaultdict(list)
    for seed_id, paths in files.items():
        alias = alias_seed_id(seed_id)
        if seed_id not in inventory_ids and alias in inventory_ids:
            seed_id = alias
        resolved[seed_id].extend(paths)
    return {seed_id: sorted(paths) for seed_id, paths in resolved.items()}


def drop_uninventoried(files, master_directory, results_directory):
    """Set aside waveform channels without StationXML rows, listing them."""
    known = set(load_master_dataset(master_directory).Seed_id)
    skipped = sorted(set(files) - known)
    if skipped:
        path = Path(results_directory) / "_workflow" / "uninventoried.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        atomic_json(path, {"seed_ids": skipped})
        print(f"Skipped {len(skipped)} channels without StationXML; see {path}")
    return {seed: paths for seed, paths in files.items() if seed in known}


def select_rows(master, files):
    return master[master.Seed_id.isin(files)].drop_duplicates(KEYS).copy()


def is_pressure(row):
    return is_pressure_unit(row.get("input_units", ""))


def build_plan(master, files, data_directory):
    """Report actual availability, absent inventory entries and source sizes."""
    ids = set(files)
    networks = sorted({seed.split(".")[0] for seed in ids})
    inventory_scope = master[master.Seed_id.str.split(".").str[0].isin(networks)]
    channels = []
    for seed, paths in sorted(files.items()):
        matches = master[master.Seed_id.eq(seed)].drop_duplicates(KEYS)
        channels.append(
            {
                "seed_id": seed,
                "files": len(paths),
                "bytes": sum(path.stat().st_size for path in paths),
                "inventory_epochs": len(matches),
                "pressure": bool(len(matches) and is_pressure(matches.iloc[0])),
                "sample_rates_hz": sorted(set(matches.sample_rate_hz)),
            }
        )
    stations = defaultdict(set)
    for seed in ids:
        net, sta, _, _ = seed.split(".")
        stations[net].add(sta)
    return {
        "data_directory": str(Path(data_directory).expanduser().resolve()),
        "networks": networks,
        "stations_by_network": {net: sorted(stas) for net, stas in stations.items()},
        "station_count": sum(map(len, stations.values())),
        "channel_count": len(ids),
        "channel_epochs": len(select_rows(master, files)),
        "waveform_files": sum(c["files"] for c in channels),
        "waveform_bytes": sum(c["bytes"] for c in channels),
        "unmatched_seed_ids": sorted(ids - set(master.Seed_id)),
        "inventory_only_seed_ids": sorted(set(inventory_scope.Seed_id) - ids),
        "channels": channels,
        "note": "Filename availability only; processing validates actual headers, "
        "response epochs and sample rates. No missing channel is fabricated.",
    }


def epoch_key(row):
    encoded = json.dumps([str(row[name]) for name in KEYS]).encode()
    return hashlib.sha256(encoded).hexdigest()


def code_digest():
    digest = hashlib.sha256()
    directory = Path(__file__).parent
    for path in sorted([*directory.glob("*.py"), *directory.glob("reference_data/*")]):
        if path.is_file():
            digest.update(path.name.encode())
            digest.update(path.read_bytes())
    return digest.hexdigest()


def file_digest(path):
    path = Path(path)
    if not path.is_file():
        return "missing"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fingerprint(stage, row, paths, options, code):
    """Resume only when code, input products, sources and scientific inputs agree."""
    inputs = {
        "code": code,
        "stage": stage,
        "options": options,
        "epoch": [str(row[name]) for name in KEYS],
        "inventory_sha256": file_digest(row["inventory_file"]),
        "sample_rate_hz": row["sample_rate_hz"],
    }
    if stage in {
        "ppsd",
        "paired",
    }:
        inputs["waveforms"] = [
            [str(p), p.stat().st_size, p.stat().st_mtime_ns] for p in paths
        ]
    if stage != "ppsd":
        source = row.get("ppsd_npz_file", "")
        if stage != "paired":
            inputs["ppsd_sha256"] = file_digest(source)
        inputs["ppsd_path"] = source
        inputs["response_sha256"] = file_digest(Path(source).parent / "inventory.xml")
        inputs["depth_m"] = row.get("Profondeur", "")
        inputs["scores"] = [row.get(name, "") for name in RESULT_COLUMNS]
    return hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()


def ppsd_job(row, paths, results_directory, options):
    """Validate only this channel's headers, containing failures to its epoch."""
    try:
        if _CANCEL_EVENT is not None:
            options = options | {"cancel_event": _CANCEL_EVENT}
        records = index_waveform_files(paths)
        return process_task(row, records, results_directory, options)
    except Exception as exc:
        return {name: row[name] for name in KEYS} | {
            "ppsd_status": "failed",
            "ppsd_error": f"{type(exc).__name__}: {exc}",
        }


def diagnostic_job(
    stage,
    row,
    paths,
    master_directory,
    options,
    reference_rows=None,
    dataset_cache=None,
):
    """Use existing scientific estimators one epoch at a time, without Excel writes."""
    net, sta, _, _ = row["Seed_id"].split(".")
    selected = pd.DataFrame([row])
    common = dict(
        network=net,
        station=sta,
        selected_rows=selected,
        update_master=False,
        return_updates=True,
        show=False,
    )
    try:
        if row.get("ppsd_status") not in {"complete", "partial"} or not row.get(
            "ppsd_npz_file"
        ):
            raise ValueError("No usable PPSD: run the ppsd stage for this epoch first")
        if stage == "plots":
            source = Path(row["ppsd_npz_file"])
            output = source.parent.parent / "plots"
            plot_station(
                Path(master_directory) / "master_inventory.csv",
                network=net,
                station=sta,
                selected_rows=selected,
                output_directory=output,
                show=False,
            )
            names = {
                "histogram": "ppsd",
                "spectrogram": "spectrogram",
                "lines": "all_psd_lines",
            }
            update = {name: row[name] for name in KEYS} | {
                "ppsd_plot_status": "complete",
                **{
                    f"ppsd_{name}_plot_file": str(
                        output / f"{source.parent.name}__{suffix}.png"
                    )
                    for name, suffix in names.items()
                },
            }
            return (
                "complete",
                [update],
                [v for k, v in update.items() if k.endswith("_plot_file")],
            )
        if stage == "pressure":
            summaries, updates = run_pressure_test(master_directory, **common)
            paths = list(summaries[0]["files"].values()) + [
                updates[0]["pressure_expected_result_file"]
            ]
        elif stage == "nearby":
            summaries, updates = run_nearby_test(
                master_directory,
                criteria=Criteria(**options["criteria"]),
                reference_rows=reference_rows,
                dataset_cache=dataset_cache,
                **common,
            )
            paths = list(summaries[0]["files"].values()) + [
                updates[0]["nearby_result_file"]
            ]
        elif stage == "paired":
            summary, updates = run_paired_spectra(
                master_directory,
                records=index_waveform_files(paths),
                reference_rows=reference_rows,
                **{k: v for k, v in options.items() if k != "paired_pool_signature"},
                **common,
            )
            summaries = [summary]
            paths = list(summary["files"].values())
        elif stage == "infragravity":
            summaries, updates = run_infragravity_test(master_directory, **common)
            paths = list(summaries[0]["files"].values()) + [
                updates[0]["infragravity_result_file"]
            ]
        else:
            raise ValueError(f"Unknown diagnostic stage: {stage}")
        return summaries[0]["status"], updates, paths
    finally:
        # A batch must not keep hundreds of closed figures or PSD curves alive.
        plt.close("all")


def run_stage(
    stage,
    master_directory,
    results_directory,
    files,
    *,
    workers=4,
    force=False,
    options=None,
    pressure_only=False,
    channels=None,
):
    """Checkpoint results, isolate failed epochs and let only the parent patch Excel.

    ``pressure_only`` restricts the PPSD and plot stages to pressure channels;
    the diagnostic stages always select pressure channels. ``channels`` keeps
    only channel codes matching one of these shell patterns (e.g. ``?H?``).
    """
    if stage not in UPDATE_FUNCTIONS:
        raise ValueError(f"Unknown processing stage: {stage}")
    if workers < 1:
        raise ValueError("Workers must be a positive integer")
    results_directory = Path(results_directory).expanduser().resolve()
    report_path = results_directory / "_workflow" / f"{stage}.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report = (
        json.loads(report_path.read_text())
        if report_path.exists()
        else {"stage": stage, "entries": {}}
    )
    master = load_master_dataset(master_directory)
    missing = set(files) - set(master.Seed_id)
    if missing:
        raise ValueError(
            f"Missing StationXML rows for {sorted(missing)}; "
            "run inventory stage after adding XML"
        )
    selected = select_rows(master, files)
    if pressure_only or stage in {
        "pressure",
        "infragravity",
        "nearby",
        "paired",
    }:
        selected = selected[[is_pressure(row) for _, row in selected.iterrows()]]
    if channels:
        codes = selected.Seed_id.str.rsplit(".", n=1).str[-1]
        selected = selected[[any(fnmatchcase(c, p) for p in channels) for c in codes]]
    rows = selected.to_dict("records")
    code = code_digest()
    options = options or {}
    references, dataset_cache = None, {}
    if stage == "paired":
        assert_workflows_finished(results_directory, ignore_stage="paired")
        references = master
        options = {
            "start": None,
            "days": 7,
            "segment_seconds": 3600,
            "frequency_min": 0.01,
            "frequency_max": 0.5,
        } | options
        # Both channels, epochs, snapshots and source metadata invalidate the cache.
        scope = {seed.rsplit(".", 2)[0] for seed in files}
        dependencies = []
        for _, candidate in master.iterrows():
            if candidate.Seed_id.rsplit(".", 2)[0] not in scope:
                continue
            source = candidate.get("ppsd_npz_file", "")
            dependencies.append(
                [
                    *[str(candidate[n]) for n in KEYS],
                    file_digest(candidate.inventory_file),
                    file_digest(Path(source).parent / "inventory.xml"),
                    str(candidate.get("ppsd_status", "")),
                    str(source),
                ]
            )
        options["paired_pool_signature"] = hashlib.sha256(
            json.dumps(dependencies, sort_keys=True).encode()
        ).hexdigest()
    if stage == "nearby":
        ensure_master_idle(results_directory)
        references = pressure_rows(master)
        options = {"criteria": asdict(Criteria())} | options
        options["reference_pool_fingerprint"] = pool_fingerprint(references)
    report["selected_epoch_keys"] = [epoch_key(row) for row in rows]
    pending, updates, completed = [], [], []

    def flush():
        if updates:
            UPDATE_FUNCTIONS[stage](master_directory, updates)
            updates.clear()
        # Record success after master updates; interrupted tasks safely retry.
        atomic_json(report_path, report)

    def record(row, stamp, status, patch, artifacts, error=None):
        entry = {
            "seed_id": row["Seed_id"],
            "epoch": {n: row[n] for n in KEYS},
            "fingerprint": stamp,
            "status": status,
            "updates": patch,
            "files": artifacts,
            "error": error,
        }
        report["entries"][epoch_key(row)] = entry
        completed.append(entry)
        updates.extend(patch)
        print(
            f"{stage}: {row['Seed_id']} -> {status}" + (f" ({error})" if error else ""),
            flush=True,
        )
        if len(completed) % 4 == 0:
            flush()

    for row in rows:
        paths = files[row["Seed_id"]]
        if stage == "paired":
            net, sta, _, _ = row["Seed_id"].split(".")
            paths = [
                p
                for seed, source_paths in files.items()
                if seed.startswith(f"{net}.{sta}.")
                for p in source_paths
            ]
        stamp = fingerprint(stage, row, paths, options, code)
        cached = report["entries"].get(epoch_key(row))
        if (
            not force
            and cached
            and cached["fingerprint"] == stamp
            and cached["status"] in SUCCESS
            and cached["files"]
            and all(Path(p).is_file() for p in cached["files"])
        ):
            completed.append(cached)
            updates.extend(cached["updates"])
            print(
                f"{stage}: {row['Seed_id']} already complete; reused saved results",
                flush=True,
            )
        else:
            pending.append((row, paths, stamp))
    # A previous completed run must not look idle while this batch is writing.
    report.pop("status_counts", None)
    flush()
    try:
        if stage == "ppsd":

            def accept(row, stamp, patch):
                artifacts = [
                    patch[name]
                    for name in ("ppsd_npz_file", "ppsd_provenance_file")
                    if patch.get(name)
                ]
                record(
                    row,
                    stamp,
                    patch["ppsd_status"],
                    [patch],
                    artifacts,
                    patch.get("ppsd_error"),
                )

            # Longest channels first, so one late large epoch cannot idle the
            # other workers. Scheduling only; each epoch is computed unchanged.
            pending.sort(key=lambda task: -sum(Path(p).stat().st_size for p in task[1]))
            if workers == 1:
                for row, paths, stamp in pending:
                    accept(row, stamp, ppsd_job(row, paths, results_directory, options))
            elif pending:
                for name in THREAD_VARIABLES:
                    os.environ.setdefault(name, "1")
                context = get_context("spawn")
                cancellation = context.Event()
                with ProcessPoolExecutor(
                    max_workers=min(workers, len(pending)),
                    mp_context=context,
                    initializer=initialize_worker,
                    initargs=(cancellation,),
                ) as pool:
                    futures = {
                        pool.submit(ppsd_job, row, paths, results_directory, options): (
                            row,
                            stamp,
                        )
                        for row, paths, stamp in pending
                    }
                    try:
                        for future in as_completed(futures):
                            row, stamp = futures[future]
                            try:
                                patch = future.result()
                            except Exception as exc:
                                patch = {n: row[n] for n in KEYS} | {
                                    "ppsd_status": "failed",
                                    "ppsd_error": str(exc),
                                }
                            accept(row, stamp, patch)
                    except KeyboardInterrupt:
                        cancellation.set()
                        for future in futures:
                            future.cancel()
                        print(
                            "Stopping after current chunks; rerun to resume", flush=True
                        )
                        raise
        else:
            for row, paths, stamp in pending:
                try:
                    status, patch, artifacts = diagnostic_job(
                        stage,
                        row,
                        paths,
                        master_directory,
                        options,
                        references,
                        dataset_cache,
                    )
                    for update in patch:
                        update[f"{PREFIXES[stage]}_error"] = None
                    record(row, stamp, status, patch, artifacts)
                except Exception as exc:
                    prefix = PREFIXES[stage]
                    message = f"{type(exc).__name__}: {exc}"
                    patch = {n: row[n] for n in KEYS} | {
                        f"{prefix}_status": "failed",
                        f"{prefix}_error": message,
                    }
                    record(row, stamp, "failed", [patch], [], message)
    finally:
        flush()
    report["status_counts"] = dict(Counter(entry["status"] for entry in completed))
    atomic_json(report_path, report)
    print(f"{stage}: {report['status_counts']}; details: {report_path}", flush=True)
    return report


def view_plots(master_directory, files, kind, *, show=True, workflow_path=None):
    """Open selected saved PNGs; do not recalculate or write any master fields."""
    prefixes = {
        "ppsd": "ppsd_",
        "pressure": "pressure_expected_",
        "infragravity": "infragravity_",
        "nearby": "nearby_",
        "paired": "paired_",
        "rayleigh": "rayleigh_",
    }
    rows = select_rows(load_master_dataset(master_directory), files)
    paths = []
    if workflow_path is not None:
        if kind != "rayleigh":
            raise ValueError("Explicit workflow viewing currently supports rayleigh")
        workflow = json.loads(Path(workflow_path).read_text())
        selected = {tuple(str(row[k]) for k in KEYS) for _, row in rows.iterrows()}
        for entry in workflow["entries"].values():
            if tuple(str(entry["epoch"][k]) for k in KEYS) not in selected:
                continue
            value = entry["update"].get("rayleigh_plot_file", "")
            path = Path(value)
            if value and path.is_file() and path not in paths:
                paths.append(path)
    for _, row in rows.iterrows():
        if workflow_path is not None:
            break
        for name, value in row.items():
            if not name.startswith(prefixes[kind]) or not str(value).endswith(".png"):
                continue
            if kind == "infragravity" and name.startswith("infragravity_pilot_"):
                continue
            path = Path(value)
            if path.is_file() and path not in paths:
                paths.append(path)
    if not paths:
        raise ValueError(f"No saved {kind} plots for the selected scope")
    for path in paths:
        figure, axis = plt.subplots(figsize=(12, 8))
        axis.imshow(plt.imread(path))
        axis.axis("off")
        figure.subplots_adjust(left=0, right=1, bottom=0, top=1)
        print(f"Opened {path}", flush=True)
    if show:
        plt.show()
    else:
        plt.close("all")
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=(*STAGES, "rayleigh", "all", "view"))
    parser.add_argument("--data-dir", type=Path, default=DATA_DIRECTORY)
    parser.add_argument("--inventory-dir", type=Path, default=XML_DIRECTORY)
    parser.add_argument("--master-dir", type=Path, default=OUTPUT_DIRECTORY)
    parser.add_argument("--results-dir", type=Path, default=RESULT_DIRECTORY)
    parser.add_argument(
        "--network", action="append", help="Omit to discover every network in Data"
    )
    parser.add_argument("--station", action="append")
    parser.add_argument("--seed-id", action="append")
    parser.add_argument(
        "--workers",
        type=worker_count,
        default=4,
        help="Parallel PPSD channel workers; `max` uses all CPU cores but one",
    )
    parser.add_argument(
        "--pressure-only",
        action="store_true",
        help="Restrict the ppsd and plots stages to pressure channels",
    )
    parser.add_argument(
        "--channels",
        action="append",
        help="Keep channel codes matching this shell pattern (repeatable), e.g. ?H?",
    )
    parser.add_argument(
        "--skip-uninventoried",
        action="store_true",
        help="Skip waveform channels without StationXML rows instead of stopping",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Retry saved stages; PPSD still resumes checkpoints",
    )
    parser.add_argument("--ppsd-start")
    parser.add_argument("--ppsd-end")
    parser.add_argument("--ppsd-length", type=float, default=10800)
    parser.add_argument("--minimum-frequency", type=float, default=0.001)
    parser.add_argument("--overlap", type=float, default=0.5)
    parser.add_argument("--chunk-seconds", type=float, default=86400)
    parser.add_argument("--rebuild-ppsd", action="store_true")
    parser.add_argument(
        "--view-kind",
        choices=(
            "ppsd",
            "pressure",
            "infragravity",
            "nearby",
            "paired",
            "rayleigh",
        ),
        default="ppsd",
    )
    parser.add_argument("--paired-start")
    parser.add_argument("--paired-days", type=int, default=7)
    parser.add_argument("--paired-segment-seconds", type=float, default=3600)
    parser.add_argument("--paired-frequency-min", type=float, default=0.01)
    parser.add_argument("--paired-frequency-max", type=float, default=0.5)
    parser.add_argument("--paired-vertical-id")
    parser.add_argument("--rayleigh-catalog", type=Path)
    model = parser.add_mutually_exclusive_group()
    model.add_argument("--rayleigh-dispersion-file", type=Path)
    model.add_argument(
        "--rayleigh-earth-model-file",
        type=Path,
        help="Layer CSV; generate phase/group dispersion with disba",
    )
    parser.add_argument("--rayleigh-frequency-min", type=float, default=0.03)
    parser.add_argument("--rayleigh-frequency-max", type=float, default=0.07)
    parser.add_argument("--rayleigh-native-rate-tolerance-ppm", type=float, default=0)
    parser.add_argument(
        "--rayleigh-event-resample-rate",
        type=float,
        help="Explicit working-copy anti-alias rate in Hz; requires no tilt",
    )
    parser.add_argument(
        "--rayleigh-pressure-sign", type=int, choices=(-1, 1), default=-1
    )
    parser.add_argument(
        "--rayleigh-min-events",
        type=int,
        default=3,
        help="Minimum accepted events for station estimate; not a cap in all mode",
    )
    parser.add_argument(
        "--rayleigh-max-events",
        type=int,
        default=8,
        help="Attempt cap per magnitude level, only with bounded event sampling",
    )
    parser.add_argument(
        "--rayleigh-event-sampling",
        choices=("all", "bounded"),
        default="all",
        help="all: exhaust eligible catalogue; bounded: historical pilot limits",
    )
    parser.add_argument(
        "--rayleigh-no-master-update",
        action="store_true",
        help="Keep Rayleigh comparison results separate; preserve master fields",
    )
    parser.add_argument(
        "--view-workflow",
        type=Path,
        help="View Rayleigh plots from a separate workflow rather than the master",
    )
    parser.add_argument(
        "--rayleigh-tilt-method",
        choices=("none", "tiskitpy", "local_eigen"),
        default="tiskitpy",
        help="none: pressure/BHZ only, without quiet training or tilt correction",
    )
    add_criteria_arguments(parser)
    args = parser.parse_args()
    if args.stage != "view":
        plt.switch_backend("Agg")
    stages = STAGES if args.stage == "all" else (args.stage,)
    failed = False
    files = None
    for stage in stages:
        if stage == "inventory":
            print(
                synchronize_inventory(args.inventory_dir, args.master_dir), flush=True
            )
            continue
        if files is None:
            files = apply_location_aliases(
                scan_files(args.data_dir, args.network, args.station),
                set(load_master_dataset(args.master_dir).Seed_id),
            )
            if args.seed_id:
                files = {s: p for s, p in files.items() if s in args.seed_id}
                if not files:
                    raise ValueError("No matching SDS waveform files for --seed-id")
            if args.skip_uninventoried and stage != "check":
                files = drop_uninventoried(files, args.master_dir, args.results_dir)
        if stage == "check":
            plan = build_plan(
                load_master_dataset(args.master_dir), files, args.data_dir
            )
            path = args.results_dir / "_workflow/plan.json"
            path.parent.mkdir(parents=True, exist_ok=True)
            atomic_json(path, plan)
            print(
                f"Networks: {plan['networks']}; {plan['station_count']} stations; "
                f"{plan['channel_count']} channels; {plan['waveform_files']} files; "
                f"{plan['waveform_bytes'] / 1e9:.2f} GB\nSaved {path}",
                flush=True,
            )
            if plan["unmatched_seed_ids"]:
                raise ValueError(
                    f"Add inventories for {plan['unmatched_seed_ids']} "
                    "and rerun inventory"
                )
        elif stage == "rayleigh":
            # Explicit stage only: network `all` must not start a calibration batch.
            from compy_pressure.rayleigh_network import run_network

            report = run_network(
                args.master_dir,
                args.results_dir,
                files,
                force=args.force,
                data_dir=args.data_dir,
                update_master=not args.rayleigh_no_master_update,
                options=dict(
                    catalog=args.rayleigh_catalog,
                    dispersion_file=args.rayleigh_dispersion_file,
                    earth_model_file=args.rayleigh_earth_model_file,
                    band=[args.rayleigh_frequency_min, args.rayleigh_frequency_max],
                    pressure_sign=args.rayleigh_pressure_sign,
                    min_events=args.rayleigh_min_events,
                    max_events=args.rayleigh_max_events,
                    event_sampling=args.rayleigh_event_sampling,
                    tilt_method=args.rayleigh_tilt_method,
                    event_resample_rate=args.rayleigh_event_resample_rate,
                    native_rate_tolerance_ppm=args.rayleigh_native_rate_tolerance_ppm,
                ),
            )
            failed |= "failed" in report["status_counts"]
        elif stage == "view":
            if not args.station and not args.seed_id:
                raise ValueError(
                    "Select --station or --seed-id when opening plot windows"
                )
            view_plots(
                args.master_dir,
                files,
                args.view_kind,
                workflow_path=args.view_workflow,
            )
        else:
            options = {}
            if stage == "ppsd":
                options = dict(
                    start=args.ppsd_start,
                    end=args.ppsd_end,
                    ppsd_length=args.ppsd_length,
                    overlap=args.overlap,
                    minimum_frequency=args.minimum_frequency,
                    chunk_seconds=args.chunk_seconds,
                    rebuild=args.rebuild_ppsd,
                )
            elif stage == "paired":
                options = dict(
                    start=args.paired_start,
                    days=args.paired_days,
                    segment_seconds=args.paired_segment_seconds,
                    frequency_min=args.paired_frequency_min,
                    frequency_max=args.paired_frequency_max,
                    vertical_id=args.paired_vertical_id,
                )
            elif stage == "nearby":
                options = {"criteria": asdict(criteria_from_args(args))}
            report = run_stage(
                stage,
                args.master_dir,
                args.results_dir,
                files,
                workers=args.workers,
                force=args.force,
                options=options,
                pressure_only=args.pressure_only,
                channels=args.channels,
            )
            failed |= any(status not in SUCCESS for status in report["status_counts"])
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
