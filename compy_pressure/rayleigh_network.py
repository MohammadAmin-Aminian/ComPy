"""Bounded, restartable network Rayleigh comparisons; no response or score edits."""

import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt

from compy_pressure.archive import index_waveform_files
from compy_pressure.inventory import (
    _lock_file,
    build_master_dataset,
    load_master_dataset,
    save_master_dataset,
    update_rayleigh_results,
)
from compy_pressure.paired_spectra import (
    KEYS,
    assert_workflows_finished,
    response_snapshot_path,
)
from compy_pressure.ppsd import atomic_json, is_pressure_unit
from compy_pressure.rayleigh_calibration import (
    require_clean_rotator,
    review_assumptions,
    run_calibration,
    select_calibration_channels,
    water_ratio,
)
from compy_pressure.rayleigh_dispersion import load_dispersion

FINISHED = {
    "candidate_comparison_under_assumed_sign",
    "unresolved",
    "unresolved_metadata",
    "reference_band_unusable",
}


def digest(path):
    path = Path(path)
    return (
        hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else "missing"
    )


def epoch_key(row):
    return hashlib.sha256(json.dumps([str(row[k]) for k in KEYS]).encode()).hexdigest()


def dependency_stamp(row, master, files, options, model_provenance=None):
    """Station epochs/responses, raw stats, catalogue and model control reuse."""
    prefix = row.Seed_id.rsplit(".", 2)[0]
    dependencies = []
    for _, partner in master[master.Seed_id.str.startswith(prefix + ".")].iterrows():
        saved = partner.get("ppsd_npz_file", "")
        dependencies.append(
            [
                *[str(partner[k]) for k in KEYS],
                saved,
                partner.get("ppsd_status", ""),
                digest(partner.inventory_file),
                digest(response_snapshot_path(partner)),
            ]
        )
    source = [
        [str(p), p.stat().st_size, p.stat().st_mtime_ns]
        for seed, paths in sorted(files.items())
        if seed.startswith(prefix + ".")
        for p in sorted(map(Path, paths))
    ]
    code = {p.name: digest(p) for p in Path(__file__).parent.glob("*.py")}
    value = dict(
        epoch=[str(row[k]) for k in KEYS],
        depth=row.Profondeur,
        coordinates=[str(row.Lat), str(row.Lon)],
        sensor=row.get("sensor_description", ""),
        dependencies=dependencies,
        waveforms=source,
        options=options,
        model_provenance=model_provenance,
        code=code,
        catalog_sha256=digest(options["catalog"]),
        dispersion_sha256=digest(options["dispersion_file"])
        if options.get("dispersion_file")
        else None,
        earth_model_sha256=digest(options["earth_model_file"])
        if options.get("earth_model_file")
        else None,
        window_sha256=digest(options["window_overrides"])
        if options.get("window_overrides")
        else None,
    )
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def patch_for(row, summary):
    """Only rayleigh_* fields; never a scientific score or instrument sensitivity."""
    gain = summary.get("median_pressure_multiplier")
    sign = summary.get("settings", {}).get("pressure_sign", -1)
    return {k: str(row[k]) for k in KEYS} | {
        "rayleigh_status": summary["status"],
        "rayleigh_result_file": summary["files"]["result"],
        "rayleigh_plot_file": summary["files"].get("paper_comparison", ""),
        "rayleigh_candidate_events": summary.get("candidate_events", 0),
        "rayleigh_pressure_multiplier": gain if gain is not None else "",
        "rayleigh_signed_pressure_multiplier": sign * gain if gain is not None else "",
        "rayleigh_sensitivity_multiplier": 1 / gain if gain is not None else "",
        "rayleigh_gain_of_median_ratio": summary.get("paper_comparison", {}).get(
            "candidate_pressure_multiplier"
        )
        or "",
        "rayleigh_pressure_sign": sign,
        "rayleigh_validation_status": "provisional_needs_review",
        "rayleigh_review_flags": ";".join(f["id"] for f in summary["review_flags"]),
        "rayleigh_error": summary.get("reason", ""),
    }


def calibration_master(inventory, master_dir, network, station, *, channels=None):
    """Master table for a calibration without PPSDs; return (master_dir, master).

    The calibration reads each channel's response from a saved snapshot,
    which the PPSD stage normally writes. Here each channel's own StationXML
    file is copied as its snapshot. ``inventory`` is a StationXML file or a
    folder of them. ``channels`` keeps only these channel codes of the
    station, e.g. the pressure and vertical codes without tilt cleaning, which
    also avoids pairing errors from wrong horizontal dips. A saved table is
    reused; a changed StationXML needs a new ``master_dir``.
    """
    inventory, master_dir = Path(inventory), Path(master_dir)
    sources = sorted(inventory.glob("*.xml")) if inventory.is_dir() else [inventory]
    copies = master_dir / "stationxml"
    if (master_dir / "master_inventory.csv").exists():
        if any((copies / s.name).read_bytes() != s.read_bytes() for s in sources):
            raise ValueError(f"{inventory} changed; use a new master folder")
        return master_dir, load_master_dataset(master_dir)
    copies.mkdir(parents=True, exist_ok=True)
    for source in sources:
        shutil.copy2(source, copies / source.name)
    master = build_master_dataset(copies).astype(object).fillna("")
    rows = master.Seed_id.str.startswith(f"{network}.{station}.")
    if not rows.any():
        raise ValueError(f"{network}.{station} is not in {inventory}")
    if channels is not None:
        code = master.Seed_id.str.rsplit(".", n=1).str[-1]
        master = master[~rows | code.isin(channels)].reset_index(drop=True)
        rows = master.Seed_id.str.startswith(f"{network}.{station}.")
    master["rayleigh_response_snapshot_file"] = ""
    for index in master.index[rows]:
        snapshot = master_dir / "responses" / master.at[index, "Seed_id"]
        snapshot.mkdir(parents=True, exist_ok=True)
        snapshot = snapshot / "inventory.xml"
        shutil.copy2(master.at[index, "inventory_file"], snapshot)
        master.at[index, "rayleigh_response_snapshot_file"] = str(snapshot)
    save_master_dataset(master, master_dir)
    return master_dir, master


def run_network(
    master_dir,
    results_dir,
    files,
    *,
    force=False,
    update_master=True,
    options=None,
    data_dir=None,
):
    """Serial station jobs; checkpoint each result and append master fields once."""
    results_dir = Path(results_dir).resolve()
    workflow = results_dir / "_workflow"
    workflow.mkdir(parents=True, exist_ok=True)
    # Prevent two Rayleigh writers; a leftover lock file is harmless after release.
    with (workflow / "rayleigh.lock").open("a+") as lock:
        if not _lock_file(lock, wait=False):
            raise RuntimeError("Another Rayleigh network run is active")
        return _run_network(
            master_dir,
            results_dir,
            files,
            force,
            update_master,
            options or {},
            data_dir,
        )


def _run_network(
    master_dir, results_dir, files, force, update_master, options, data_dir
):
    assert_workflows_finished(results_dir, ignore_stage="rayleigh")
    options = {
        "band": [0.03, 0.07],
        "pressure_sign": -1,
        "tilt_method": "tiskitpy",
        "min_events": 3,
        "max_events": 8,
        "event_sampling": "all",
    } | options
    options = {k: str(v) if isinstance(v, Path) else v for k, v in options.items()}
    if options["tilt_method"] not in {"none", "tiskitpy", "local_eigen"}:
        raise ValueError("Unknown tilt method")
    if options["event_sampling"] not in {"all", "bounded"}:
        raise ValueError("Unknown event sampling policy")
    for key in ("catalog",):
        if not options.get(key) or not Path(options[key]).is_file():
            raise ValueError(f"Supply an existing {key} file")
    resolved_dispersion = load_dispersion(
        dispersion_file=options.get("dispersion_file"),
        earth_model_file=options.get("earth_model_file"),
        band=options["band"],
    )
    # Require the real default package before starting a network-wide computation.
    if options["tilt_method"] == "tiskitpy":
        require_clean_rotator()
    master = load_master_dataset(master_dir)
    before = {
        str(Path(master_dir) / f"master_inventory.{ext}"): digest(
            Path(master_dir) / f"master_inventory.{ext}"
        )
        for ext in ("csv", "xlsx")
    }
    selected = master[
        master.Seed_id.isin(files) & master.input_units.map(is_pressure_unit)
    ]
    selected = selected.drop_duplicates(KEYS).sort_values(KEYS)
    report_path = results_dir / "_workflow/rayleigh.json"
    report = (
        json.loads(report_path.read_text())
        if report_path.exists()
        else {
            "stage": "rayleigh",
            "entries": {},
        }
    )
    report.update(
        execution_status="running",
        options=options,
        dispersion_model=resolved_dispersion[1],
        selected_epoch_keys=[epoch_key(r) for _, r in selected.iterrows()],
        implementation_status="working_method_with_provisional_assumptions",
    )
    report.pop("status_counts", None)
    atomic_json(report_path, report)
    completed, patches = [], []
    record_cache = {}
    for _, row in selected.iterrows():
        key = epoch_key(row)
        stamp = dependency_stamp(row, master, files, options, resolved_dispersion[1])
        cached = report["entries"].get(key)
        if (
            not force
            and cached
            and cached["fingerprint"] == stamp
            and cached["status"] in FINISHED
            and all(Path(p).is_file() for p in cached["files"])
        ):
            completed.append(cached)
            patches.append(cached["update"])
            print(f"rayleigh: {row.Seed_id} reused saved result", flush=True)
            continue
        net, sta, _, _ = row.Seed_id.split(".")
        summary = {
            "station": f"{net}.{sta}",
            "status": "unresolved_metadata",
            "epoch_keys": [{k: str(row[k]) for k in KEYS}],
            "implementation_status": "working_method_with_provisional_assumptions",
            "settings": options,
            "review_flags": review_assumptions(row, options["tilt_method"]),
            "candidate_events": 0,
            "events": [],
            "files": {},
        }
        try:
            rows, channels = select_calibration_channels(
                master,
                net,
                sta,
                row,
                options["tilt_method"],
                options.get("event_resample_rate"),
            )
            # This safeguard is not relaxed merely to produce a network gain.
            dispersion = resolved_dispersion[0]
            try:
                water_ratio(
                    [options["band"][1]],
                    dispersion[:, 1].max(),
                    float(row.Profondeur)
                    if options.get("depth") is None
                    else float(options["depth"]),
                    options.get("density", 1028),
                    options.get("sound_speed", 1500),
                )
            except ValueError:
                summary["status"] = "reference_band_unusable"
                raise
            ids = {r.Seed_id for r in rows}
            signature = tuple(sorted(ids))
            if signature not in record_cache:
                record_cache[signature] = index_waveform_files(
                    [p for seed in ids for p in files.get(seed, [])]
                )
            summary = run_calibration(
                master_dir=master_dir,
                data_dir=data_dir,
                network=net,
                station=sta,
                pressure_epoch=row,
                master_dataset=master,
                waveform_records=record_cache[signature],
                _resolved_dispersion=resolved_dispersion,
                **options,
            )
            if any("timing_review" in e for e in summary["events"]):
                summary["review_flags"].append(
                    {
                        "id": "clock_alignment_needs_review",
                        "status": "needs_review_or_testing",
                        "note": "Misaligned grids were excluded without interpolation",
                    }
                )
                atomic_json(Path(summary["files"]["result"]), summary)
        except (ValueError, FileNotFoundError) as exc:
            summary["reason"] = str(exc)
        except Exception as exc:
            summary.update(status="failed", reason=f"{type(exc).__name__}: {exc}")
        finally:
            plt.close("all")
        if "result" not in summary["files"]:
            output = results_dir / "_workflow/rayleigh_unresolved"
            output.mkdir(exist_ok=True)
            path = output / f"{key[:16]}_{uuid4().hex[:8]}.json"
            summary["files"]["result"] = str(path)
            atomic_json(path, summary)
        artifacts = [summary["files"]["result"]]
        if summary["files"].get("paper_comparison"):
            artifacts.append(summary["files"]["paper_comparison"])
        # Event data/plots and reviewed-window template are also required for reuse.
        artifacts.extend(
            p for e in summary["events"] for p in e.get("files", {}).values()
        )
        if summary["files"].get("dispersion"):
            artifacts.append(summary["files"]["dispersion"])
        if summary["files"].get("window_review_template"):
            artifacts.append(summary["files"]["window_review_template"])
        patch = patch_for(row, summary)
        entry = dict(
            seed_id=row.Seed_id,
            epoch={k: str(row[k]) for k in KEYS},
            fingerprint=stamp,
            status=summary["status"],
            update=patch,
            files=artifacts,
            error=summary.get("reason"),
        )
        report["entries"][key] = entry
        completed.append(entry)
        patches.append(patch)
        atomic_json(report_path, report)
        print(
            f"rayleigh: {row.Seed_id} -> {summary['status']}; "
            f"{summary.get('candidate_events', 0)} candidates",
            flush=True,
        )
    if any(digest(p) != h for p, h in before.items()):
        raise RuntimeError(
            "Master changed during the run; results saved, rerun to resume safely"
        )
    if update_master:
        update_rayleigh_results(master_dir, patches)
    report.update(
        execution_status="complete",
        master_updated=bool(update_master),
        status_counts=dict(Counter(e["status"] for e in completed)),
    )
    atomic_json(report_path, report)
    print(f"rayleigh: {report['status_counts']}; {report_path}", flush=True)
    return report
