"""Execute only a coverage plan's missing tasks, with locking and safe resume.

Without ``execute`` this is a dry run: it builds the plan, reports which tasks
are already complete in the campaign and which would run, and writes nothing
except the header cache. With ``execute`` it holds an exclusive campaign lock,
runs each task through the existing estimators (``network.ppsd_job`` and
``network.diagnostic_job``), records per-task state and logs, isolates failures,
and verifies that protected inputs are unchanged afterwards.

All writes stay inside the campaign folder. Reused saved products are never
opened for writing: checks on them run on verified clone copies (``views``),
because the estimators save outputs beside the PPSD file they read.
"""

import contextlib
import json
import os
import shutil
import socket
import subprocess
import sys
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from datetime import datetime, timezone
from multiprocessing import get_context
from pathlib import Path

import pandas as pd

from compy_pressure import network
from compy_pressure.coverage import (
    KEYS,
    build_plan,
    digest,
    fresh_name,
    json_hash,
    ppsd_row,
    public,
    summary_text,
    utc_now,
    write_plan,
)
from compy_pressure.inventory import RESULT_COLUMNS
from compy_pressure.nearby_test import Criteria
from compy_pressure.ppsd import atomic_json

try:
    import fcntl
except ImportError:  # Windows: the module imports, but campaigns do not run
    fcntl = None

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

DONE = "complete"
ORDER = ("ppsd", "ppsd_merge", "expected_pressure", "infragravity", "nearby")


class CampaignLocked(RuntimeError):
    """Another invocation is executing this campaign."""


@contextlib.contextmanager
def campaign_lock(campaign_dir):
    """Exclusive, non-waiting lock: a second invocation must not duplicate work."""
    if fcntl is None:
        raise OSError("Whole-drive campaigns need macOS or Linux (fcntl file locks)")
    campaign_dir.mkdir(parents=True, exist_ok=True)
    path = campaign_dir / ".campaign.lock"
    with path.open("a+") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            handle.seek(0)
            raise CampaignLocked(
                f"Campaign is locked by another invocation: {handle.read().strip()}"
            ) from exc
        handle.seek(0)
        handle.truncate()
        handle.write(
            json.dumps(
                {
                    "pid": os.getpid(),
                    "host": socket.gethostname(),
                    "started_utc": datetime.now(UTC).isoformat(),
                }
            )
        )
        handle.flush()
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


# --------------------------------------------------------------------------
# Task state


def state_path(campaign_dir, task):
    return campaign_dir / "tasks" / task["task_id"] / "state.json"


def log_path(campaign_dir, task):
    return campaign_dir / "tasks" / task["task_id"] / "log.txt"


def read_state(campaign_dir, task):
    path = state_path(campaign_dir, task)
    return json.loads(path.read_text()) if path.is_file() else None


def write_state(campaign_dir, task, **fields):
    path = state_path(campaign_dir, task)
    path.parent.mkdir(parents=True, exist_ok=True)
    previous = read_state(campaign_dir, task) or {
        "task_id": task["task_id"],
        "kind": task["kind"],
        "task": public(task),
        "attempts": 0,
    }
    state = previous | fields | {"updated_utc": datetime.now(UTC).isoformat()}
    atomic_json(path, state)
    return state


def outputs_valid(state):
    return bool(state.get("outputs")) and all(
        Path(item["path"]).is_file() and digest(item["path"]) == item["sha256"]
        for item in state["outputs"]
    )


def is_done(campaign_dir, task, pool=None):
    """Complete only with an unchanged recorded output set (and nearby pool)."""
    state = read_state(campaign_dir, task)
    if not state or state.get("status") != DONE:
        return False, state
    if state.get("result_status") in {"no_data", "insufficient_data"}:
        return True, state
    if pool is not None and state.get("pool_signature") != pool:
        return False, state
    return outputs_valid(state), state


def record_outputs(paths):
    return [
        {"path": str(Path(p).resolve()), "sha256": digest(p)}
        for p in dict.fromkeys(paths)
        if p and Path(p).is_file()
    ]


# --------------------------------------------------------------------------
# Protected inputs


def protected_snapshot(plan):
    """Hashes of small/curated inputs and stat records of every waveform file."""
    hashed = {plan["config"]["config_file"]}
    hashed.update(s["path"] for s in plan["inventory_sources"])
    for registry in plan["config"]["product_registries"]:
        hashed.add(registry)
        sibling = Path(registry).with_suffix(".xlsx")
        if sibling.is_file():
            hashed.add(str(sibling))
    hashed.update(plan["evidence_source_hashes"])
    for product in plan["products"]:
        for name in ("npz", "provenance"):
            if product.get(name) and Path(product[name]).is_file():
                hashed.add(product[name])
        snapshot = Path(product["npz"]).with_name("inventory.xml")
        if snapshot.is_file():
            hashed.add(str(snapshot))
    for channel in plan["channels"]:
        for entry in channel.get("checks", {}).values():
            for saved in (entry.get("result"), entry.get("historical")):
                if saved and saved.get("path") and Path(saved["path"]).is_file():
                    hashed.add(saved["path"])
    stats = {}
    for record in plan["_files"]:
        stat = Path(record["path"]).stat()
        stats[record["path"]] = [stat.st_size, stat.st_mtime_ns]
    return {
        "sha256": {p: digest(p) for p in sorted(hashed) if Path(p).is_file()},
        "stat": stats,
    }


def compare_snapshots(before, after):
    changed = [p for p, h in before["sha256"].items() if after["sha256"].get(p) != h]
    changed += [p for p, s in before["stat"].items() if after["stat"].get(p) != s]
    return changed


# --------------------------------------------------------------------------
# Inputs prepared inside the campaign


def day_token(record):
    """Year and day-of-year from the archive file name, else the header start."""
    parts = Path(record["path"]).name.split(".")
    if record["layout"] in {"sds", "sds_mseed"}:
        return parts[5], parts[6]
    if record["layout"] == "daily":
        day = datetime.strptime(parts[4], "%Y%m%d")
    else:
        day = datetime.fromisoformat(record["start"].replace("Z", "+00:00"))
    return f"{day:%Y}", f"{day.timetuple().tm_yday:03d}"


def symlink_view(campaign_dir, records):
    """Header-named read-only links, so existing readers see true SEED IDs."""
    paths = []
    for record in records:
        seed = record["seed_id"]
        net, sta, _, cha = seed.split(".")
        year, day = day_token(record)
        name = f"{seed}.D.{year}.{day}"
        source = Path(os.path.realpath(record["path"]))
        stat = source.stat()
        if [stat.st_size, stat.st_mtime_ns] != [record["size"], record["mtime_ns"]]:
            raise RuntimeError(f"Source changed after planning: {source}")
        for copy in range(1000):
            folder = (
                campaign_dir
                / "_inputs/waveforms"
                / net
                / sta
                / (f"{cha}.D" if not copy else f"{cha}.D.copy{copy}")
            )
            link = folder / name
            if link.is_symlink() and os.path.realpath(link) == str(source):
                break
            if not link.exists() and not link.is_symlink():
                folder.mkdir(parents=True, exist_ok=True)
                link.symlink_to(source)
                break
        else:
            raise RuntimeError(f"Too many same-day copies for {seed}")
        paths.append(link)
    return paths


def clone(source, destination):
    """Copy-on-write clone where supported (APFS), otherwise a plain copy."""
    if sys.platform == "darwin":
        result = subprocess.run(
            ["cp", "-c", str(source), str(destination)], capture_output=True
        )
        if result.returncode == 0:
            return
    shutil.copy2(source, destination)


def product_view(campaign_dir, npz, expected_sha256):
    """Verified copy of a saved product folder's PPSD, snapshot and provenance."""
    npz = Path(npz)
    if digest(npz) != expected_sha256:
        raise RuntimeError(f"Saved product changed since planning: {npz}")
    folder = campaign_dir / "views" / expected_sha256[:16] / npz.parent.name
    target = folder / "ppsd.npz"
    if target.is_file() and digest(target) == expected_sha256:
        return target
    folder.mkdir(parents=True, exist_ok=True)
    copied = {}
    for name in ("ppsd.npz", "inventory.xml", "processing.json"):
        source = npz.with_name(name)
        if source.is_file():
            temporary = folder / f".{name}.tmp"
            temporary.unlink(missing_ok=True)
            clone(source, temporary)
            if digest(temporary) != digest(source):
                raise RuntimeError(f"View copy differs from source: {source}")
            temporary.replace(folder / name)
            copied[name] = {"source": str(source), "sha256": digest(source)}
    atomic_json(folder / "view.json", {"sources": copied, "purpose": "read-only view"})
    return target


# --------------------------------------------------------------------------
# Workers


def ppsd_worker(row, paths, results_directory, options, log):
    """Run the existing PPSD job with this task's output captured in its log."""
    with open(log, "a") as stream:
        with contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
            print(f"--- {datetime.now(UTC).isoformat()} PPSD {row['Seed_id']}")
            return network.ppsd_job(row, paths, results_directory, options)


def run_ppsd_tasks(campaign_dir, pending, channels, workers):
    jobs = []
    for task in pending:
        channel = channels[tuple(task["key"][k] for k in KEYS)]
        try:
            paths = symlink_view(campaign_dir, channel["_files"])
        except Exception as exc:
            state = read_state(campaign_dir, task) or {}
            write_state(
                campaign_dir,
                task,
                status="failed",
                error=f"input view: {type(exc).__name__}: {exc}",
                attempts=state.get("attempts", 0) + 1,
            )
            continue
        options = {
            **task["parameters"],
            "start": task["interval"][0],
            "end": task["interval"][1],
            "rebuild": False,
        }
        output = Path(task["output_path"])
        output.mkdir(parents=True, exist_ok=True)
        jobs.append((task, ppsd_row(channel), paths, output, options))
        write_state(campaign_dir, task, status="running", started_utc=utc_now())

    def finish(task, patch):
        status = patch.get("ppsd_status")
        outputs = record_outputs(
            [patch.get("ppsd_npz_file"), patch.get("ppsd_provenance_file")]
        )
        final = (
            DONE
            if status in {"complete", "no_data", "insufficient_data"}
            else ("partial" if status == "partial" else "failed")
        )
        state = read_state(campaign_dir, task)
        write_state(
            campaign_dir,
            task,
            status=final,
            result_status=status,
            error=patch.get("ppsd_error"),
            outputs=outputs,
            updates=[patch],
            attempts=state.get("attempts", 0) + 1,
            finished_utc=utc_now(),
        )
        print(f"ppsd {task['seed_id']} {task['interval']}: {status}", flush=True)

    if workers == 1 or len(jobs) <= 1:
        for task, row, paths, output, options in jobs:
            try:
                patch = ppsd_worker(
                    row, paths, output, options, log_path(campaign_dir, task)
                )
            except KeyboardInterrupt:
                write_state(campaign_dir, task, status="interrupted")
                print("Interrupted: PPSD checkpoints are kept; rerun to resume")
                raise
            except Exception as exc:
                patch = {
                    "ppsd_status": "failed",
                    "ppsd_error": f"{type(exc).__name__}: {exc}",
                }
            finish(task, patch)
        return
    context = get_context("spawn")
    cancel = context.Event()
    with ProcessPoolExecutor(
        max_workers=min(workers, len(jobs)),
        mp_context=context,
        initializer=network.initialize_worker,
        initargs=(cancel,),
    ) as pool:
        futures = {
            pool.submit(
                ppsd_worker, row, paths, output, options, log_path(campaign_dir, task)
            ): task
            for task, row, paths, output, options in jobs
        }
        try:
            for future in as_completed(futures):
                task = futures[future]
                try:
                    patch = future.result()
                except Exception as exc:
                    patch = {"ppsd_status": "failed", "ppsd_error": repr(exc)}
                finish(task, patch)
        except KeyboardInterrupt:
            cancel.set()
            for future in futures:
                future.cancel()
            for future, task in futures.items():
                if not future.done():
                    write_state(campaign_dir, task, status="interrupted")
            print("Interrupted: PPSD checkpoints are kept; rerun to resume", flush=True)
            raise


def task_output_npz(campaign_dir, task_id):
    path = campaign_dir / "tasks" / task_id / "state.json"
    if not path.is_file():
        raise RuntimeError(f"Dependency {task_id} has not run")
    state = json.loads(path.read_text())
    npz = [o["path"] for o in state.get("outputs", []) if o["path"].endswith(".npz")]
    if state.get("status") != DONE or len(npz) != 1:
        raise RuntimeError(f"Dependency {task_id} produced no PPSD")
    return npz[0]


def run_merge(campaign_dir, task):
    """Base windows preserved; supplement windows added; no recomputation."""
    from obspy.signal import PPSD

    base = Path(task["base_product"])
    if digest(base) != task["base_product_sha256"]:
        raise RuntimeError("Base product changed since planning")
    supplements = [task_output_npz(campaign_dir, t) for t in task["dependencies"]]
    folder = Path(task["output_path"]) / base.parent.name
    folder.mkdir(parents=True, exist_ok=True)
    merged = PPSD.load_npz(str(base))
    counts = {"base_windows": len(merged.times_processed), "supplements": []}
    for path in supplements:
        before = len(merged.times_processed)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            merged.add_npz(path)
        counts["supplements"].append(
            {
                "path": path,
                "sha256": digest(path),
                "added_windows": len(merged.times_processed) - before,
                "warnings": [str(w.message) for w in caught],
            }
        )
    target = folder / "ppsd.npz"
    temporary = folder / "ppsd.tmp.npz"
    merged.save_npz(str(temporary))
    temporary.replace(target)
    clone(base.with_name("inventory.xml"), folder / "inventory.xml")
    manifest = folder / "merge.json"
    atomic_json(
        manifest,
        {
            "base": str(base),
            "base_sha256": task["base_product_sha256"],
            **counts,
            "total_windows": len(merged.times_processed),
            "note": "Derived by ObsPy PPSD.add_npz; overlapping windows are omitted.",
        },
    )
    return "complete", record_outputs([target, manifest, folder / "inventory.xml"])


def check_row(channel, npz, registry_scores):
    row = {k: v for k, v in channel.items() if not k.startswith("_")}
    row = {k: v for k, v in row.items() if not isinstance(v, (list, dict))}
    row.update(registry_scores)
    row.update(ppsd_status="complete", ppsd_npz_file=str(npz), ppsd_units="Pa^2/Hz")
    return row


def registry_scores(channel):
    """Existing manual score cells of the reused product's exact registry row."""
    for product in channel["products"]:
        if product["npz"] != channel.get("ppsd_product"):
            continue
        for entry in product["_product"].get("_registry_rows", []):
            if tuple(entry["row"].get(k, "") for k in KEYS) == tuple(
                channel[k] for k in KEYS
            ):
                return {c: entry["row"].get(c, "") for c in RESULT_COLUMNS}
    return {c: "" for c in RESULT_COLUMNS}


def resolve_target(campaign_dir, task):
    """A campaign product, or a verified view of a reused saved product."""
    source = task["ppsd_source"]
    if source.startswith(("ppsd-", "ppsd_merge-")):
        return Path(task_output_npz(campaign_dir, source))
    return product_view(campaign_dir, source, task["ppsd_source_sha256"])


def reference_pool(campaign_dir, plan, channels):
    """Available pressure products of this plan, read in place, plus a signature.

    Saved products are read where they are (never written); campaign products
    count only once complete. Members still missing are listed as exclusions.
    """
    rows, members, excluded = [], [], []
    hashes = {p["npz"]: p.get("npz_sha256", "") for p in plan["products"]}
    for channel in channels.values():
        if channel["channel_class"] != "pressure":
            continue
        target = channel.get("ppsd_target")
        if not target:
            continue
        try:
            if target.startswith(("ppsd-", "ppsd_merge-")):
                npz = task_output_npz(campaign_dir, target)
                sha = digest(npz)
            else:
                npz, sha = target, hashes.get(target, "")
        except RuntimeError:
            excluded.append(channel["Seed_id"])
            continue
        rows.append(check_row(channel, npz, registry_scores(channel)))
        members.append([*(channel[k] for k in KEYS), str(npz), sha])
    frame = pd.DataFrame(rows).astype(str) if rows else pd.DataFrame()
    signature = json_hash(sorted(members))
    return frame, signature, sorted(excluded)


def run_check(campaign_dir, task, channel, references, cache):
    npz = resolve_target(campaign_dir, task)
    row = check_row(channel, npz, registry_scores(channel))
    options = {"criteria": asdict(Criteria())} if task["stage"] == "nearby" else {}
    with open(log_path(campaign_dir, task), "a") as stream:
        with contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
            print(f"--- {datetime.now(UTC).isoformat()} {task['stage']} {npz}")
            status, updates, paths = network.diagnostic_job(
                task["stage"],
                row,
                [],
                campaign_dir / "_no_master",
                options,
                references if task["stage"] == "nearby" else None,
                cache,
            )
    return status, record_outputs(paths), updates


# --------------------------------------------------------------------------
# Orchestration


def selected(task, only_kinds=None, only_seeds=None):
    return (not only_kinds or task["kind"] in only_kinds) and (
        not only_seeds or task["seed_id"] in only_seeds
    )


def classify_tasks(campaign_dir, plan, only_kinds=None, only_seeds=None):
    done, pending, skipped = [], [], []
    channels = {tuple(c[k] for k in KEYS): c for c in plan["channels"]}
    pool = None
    if any(t["kind"] == "nearby" for t in plan["tasks"]):
        pool = reference_pool(campaign_dir, plan, channels)[1]
    for task in plan["tasks"]:
        if not selected(task, only_kinds, only_seeds):
            skipped.append(task)
            continue
        complete, _ = is_done(
            campaign_dir, task, pool if task["kind"] == "nearby" else None
        )
        (done if complete else pending).append(task)
    return done, pending, skipped


def run_campaign(
    config,
    campaign_dir,
    *,
    scope="pressure",
    full_headers=False,
    execute=False,
    workers=4,
    only_kinds=None,
    only_seeds=None,
):
    campaign_dir = Path(campaign_dir)
    if workers < 1:
        raise ValueError("Workers must be positive")
    if not execute:
        plan = build_plan(config, campaign_dir, scope=scope, full_headers=full_headers)
        done, pending, skipped = classify_tasks(
            campaign_dir, plan, only_kinds, only_seeds
        )
        print(summary_text(plan), flush=True)
        print(
            f"Dry run: {len(done)} task(s) already complete in this campaign; "
            f"{len(pending)} would run; {len(skipped)} excluded by --only-kind/seed. "
            "Nothing executed; add --execute to run.",
            flush=True,
        )
        for task in pending:
            print(
                f"  would run {task['task_id']} {task['seed_id']} "
                f"{task.get('interval', '')} ({task['reason']})"
            )
        return {"ok": True, "executed": False, "pending": len(pending)}
    with campaign_lock(campaign_dir):
        plan = build_plan(config, campaign_dir, scope=scope, full_headers=full_headers)
        run_dir = campaign_dir / "runs" / fresh_name(plan["plan_id"])
        write_plan(plan, run_dir / "plan")
        before = protected_snapshot(plan)
        report = {
            "plan_id": plan["plan_id"],
            "started_utc": datetime.now(UTC).isoformat(),
            "tasks": {},
        }
        try:
            execute_plan(campaign_dir, plan, workers, only_kinds, only_seeds, report)
        finally:
            after = protected_snapshot(plan)
            changed = compare_snapshots(before, after)
            report.update(
                finished_utc=datetime.now(UTC).isoformat(),
                protected_files_hashed=len(before["sha256"]),
                protected_files_stat_checked=len(before["stat"]),
                protected_inputs_unchanged=not changed,
                protected_changes=changed,
            )
            atomic_json(run_dir / "run.json", report)
        counts = {}
        for value in report["tasks"].values():
            counts[value] = counts.get(value, 0) + 1
        print(
            f"Run {run_dir}: {counts}; protected inputs unchanged: {not changed}",
            flush=True,
        )
        ok = not changed and all(
            v in {DONE, "already_complete", "excluded"}
            for v in report["tasks"].values()
        )
        return {"ok": ok, "executed": True, "run_dir": str(run_dir), "counts": counts}


def execute_plan(campaign_dir, plan, workers, only_kinds, only_seeds, report):
    channels = {tuple(c[k] for k in KEYS): c for c in plan["channels"]}
    by_id = {t["task_id"]: t for t in plan["tasks"]}
    outcome = report["tasks"]

    def dependencies_ok(task):
        for dependency in task["dependencies"]:
            complete, state = is_done(campaign_dir, by_id[dependency])
            if not complete or state.get("result_status") in {
                "no_data",
                "insufficient_data",
            }:
                return False
        return True

    for kind in ORDER:
        tasks = [t for t in plan["tasks"] if t["kind"] == kind]
        references = pool = excluded = None
        if kind == "nearby" and tasks:
            references, pool, excluded = reference_pool(campaign_dir, plan, channels)
        runnable = []
        for task in tasks:
            if not selected(task, only_kinds, only_seeds):
                outcome[task["task_id"]] = "excluded"
                continue
            complete, _ = is_done(campaign_dir, task, pool)
            if complete:
                outcome[task["task_id"]] = "already_complete"
            elif not dependencies_ok(task):
                outcome[task["task_id"]] = "blocked_by_dependency"
                write_state(campaign_dir, task, status="blocked_by_dependency")
            else:
                runnable.append(task)
        if not runnable:
            continue
        print(f"{kind}: {len(runnable)} task(s) to run", flush=True)
        if kind == "ppsd":
            run_ppsd_tasks(campaign_dir, runnable, channels, workers)
            for task in runnable:
                outcome[task["task_id"]] = read_state(campaign_dir, task)["status"]
            continue
        cache = {}
        for task in runnable:
            state = read_state(campaign_dir, task) or {}
            write_state(campaign_dir, task, status="running", started_utc=utc_now())
            try:
                if kind == "ppsd_merge":
                    status, outputs = run_merge(campaign_dir, task)
                    updates = []
                else:
                    channel = channels[tuple(task["key"][k] for k in KEYS)]
                    status, outputs, updates = run_check(
                        campaign_dir, task, channel, references, cache
                    )
                write_state(
                    campaign_dir,
                    task,
                    status=DONE,
                    result_status=status,
                    outputs=outputs,
                    updates=updates,
                    pool_signature=pool,
                    pool_excluded=excluded,
                    error=None,
                    attempts=state.get("attempts", 0) + 1,
                    finished_utc=utc_now(),
                )
                outcome[task["task_id"]] = DONE
                print(f"{kind} {task['seed_id']}: {status}", flush=True)
            except KeyboardInterrupt:
                write_state(campaign_dir, task, status="interrupted")
                raise
            except Exception as exc:
                with open(log_path(campaign_dir, task), "a") as stream:
                    traceback.print_exc(file=stream)
                write_state(
                    campaign_dir,
                    task,
                    status="failed",
                    error=f"{type(exc).__name__}: {exc}",
                    attempts=state.get("attempts", 0) + 1,
                    finished_utc=utc_now(),
                )
                outcome[task["task_id"]] = "failed"
                print(f"{kind} {task['seed_id']}: failed ({exc})", flush=True)
