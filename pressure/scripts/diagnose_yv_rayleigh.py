"""Read-only YV Rayleigh station diagnosis from saved products.

Reads the frozen no-tilt workflow, its saved comparison and event spectra, the
completed pole-only trial files, the exact response snapshots and, optionally,
a saved microseism phase run and raw day-file *names*. It fits no gain, starts
no calibration, opens no waveform and writes only into a new output folder.

Response hypotheses are fixed and parameter-free (no fitting):

* H0 nominal metadata;
* H1 positive-real pressure poles reflected (the completed trial factor);
* H2 H1 plus conjugate repair of unpaired complex vertical poles/zeros;
* H3 H2 with reversed relative P/Z polarity;
* H4 H2 with the positive pressure pole magnitude divided by ten.

Each hypothesis gives the phase that nominal-response processing would show if
it were true and the physical P/Z relation were real and positive (sign −1
applied, as in the saved spectra). Residuals are descriptive diagnostics.
"""

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HYPOTHESES = ("H0", "H1", "H2", "H3", "H4")
HYPOTHESIS_LABELS = {
    "H0": "nominal metadata",
    "H1": "positive-real pressure pole reflected",
    "H2": "H1 + vertical unpaired complex roots conjugate-repaired",
    "H3": "H2 + reversed relative P/Z polarity",
    "H4": "H2 + positive pressure pole magnitude / 10",
}
MICROSEISM_BAND_HZ = (0.05, 0.30)
MICROSEISM_COHERENCE = 0.5
PHASE_SCREEN_DEG = 20.0
EVENT_CATEGORIES = {
    "excluded_cluster": "selection_excluded_24h",
    "rejected_qc": "coherence_failed",
    "polarity_review": "phase_review",
    "arrival_overlap_review": "arrival_overlap_review",
    "candidate": "candidate",
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def record_hash(path, provenance, expected_sources=None, expected_sha=None):
    """Bind a read to the completed trial, and detect changes during this audit."""
    path = str(Path(path).resolve())
    actual = sha256(path)
    if expected_sources is not None:
        if path not in expected_sources:
            raise ValueError(f"Trial manifest does not bind input: {path}")
        if actual != expected_sources[path]:
            raise ValueError(f"Input changed since the completed trial: {path}")
    if expected_sha is not None and actual != expected_sha:
        raise ValueError(f"Response snapshot hash changed: {path}")
    if path in provenance and actual != provenance[path]:
        raise ValueError(f"Input changed during diagnosis: {path}")
    provenance[path] = actual
    return actual


def verify_inputs_unchanged(provenance):
    for path, expected in provenance.items():
        if sha256(path) != expected:
            raise ValueError(f"Input changed during diagnosis: {path}")


def wrap(degrees):
    return (np.asarray(degrees, dtype=float) + 180.0) % 360.0 - 180.0


def circular_mean(degrees):
    radians = np.radians(np.asarray(degrees, dtype=float))
    return float(np.degrees(np.arctan2(np.sin(radians).mean(), np.cos(radians).mean())))


def unpaired_roots(roots, rtol=1e-4):
    """Complex roots without a conjugate partner; a real response needs pairs."""
    roots = [complex(r) for r in roots]
    free = list(range(len(roots)))
    unpaired = []
    while free:
        i = free.pop(0)
        r = roots[i]
        if abs(r.imag) <= rtol * max(abs(r), 1e-30):
            continue
        partner = next(
            (j for j in free if np.isclose(roots[j], r.conjugate(), rtol=rtol)), None
        )
        if partner is None:
            unpaired.append(r)
        else:
            free.remove(partner)
    return unpaired


def conjugate_repair(roots, rtol=1e-4):
    """Replace one of each duplicated unpaired root by its conjugate.

    Returns the original list when nothing is unpaired. Raises when a lone
    unpaired root has no duplicate, because the repair would not be unique.
    """
    roots = [complex(r) for r in roots]
    bad = unpaired_roots(roots, rtol)
    if not bad:
        return roots
    repaired = list(roots)
    pending = list(bad)
    while pending:
        r = pending.pop(0)
        twin = next(
            (k for k, x in enumerate(pending) if np.isclose(x, r, rtol=rtol)), None
        )
        if twin is None:
            raise ValueError("Unpaired complex root without a duplicate")
        pending.pop(twin)
        index = max(k for k, x in enumerate(repaired) if np.isclose(x, r, rtol=rtol))
        repaired[index] = repaired[index].conjugate()
    if unpaired_roots(repaired, rtol):
        raise ValueError("Conjugate repair left unpaired roots")
    return repaired


def pz_value(frequency, poles, zeros):
    s = 2j * np.pi * np.asarray(frequency, dtype=float)
    value = np.ones(s.shape, dtype=complex)
    for z in zeros:
        value *= s - z
    for p in poles:
        value /= s - p
    return value


def stage_ratio(frequency, stages, alternative):
    """H_metadata/H_alternative for Laplace rad/s stages (A0 unchanged)."""
    ratio = np.ones(np.shape(frequency), dtype=complex)
    for (poles, zeros), (alt_poles, alt_zeros) in zip(stages, alternative, strict=True):
        ratio *= pz_value(frequency, poles, zeros) / pz_value(
            frequency, alt_poles, alt_zeros
        )
    return ratio


def reflect_positive(poles, divisor=1.0):
    if any(p.real > 0 and p.imag != 0 for p in poles):
        raise ValueError("Reflection hypothesis requires a real pressure pole")
    return [complex(-p.real / divisor, 0) if p.real > 0 else p for p in poles]


def alternative_stages(pressure, vertical):
    """Alternative (pressure, vertical, polarity_deg) stages per hypothesis."""
    repaired_vertical = [
        (conjugate_repair(p), conjugate_repair(z)) for p, z in vertical
    ]
    reflected = [(reflect_positive(p), z) for p, z in pressure]
    divided = [(reflect_positive(p, 10.0), z) for p, z in pressure]
    return {
        "H0": (pressure, vertical, 0.0),
        "H1": (reflected, vertical, 0.0),
        "H2": (reflected, repaired_vertical, 0.0),
        "H3": (reflected, repaired_vertical, 180.0),
        "H4": (divided, repaired_vertical, 0.0),
    }


def hypothesis_terms(frequency, pressure, vertical):
    """Rayleigh-convention model phase (deg) and amplitude factor per hypothesis.

    Pressure corrected with H_meta is P_true / R_P; vertical is Z_true / R_Z.
    The saved phase angle(P conj(A)) then carries -arg R_P + arg R_Z, and the
    saved |P|/|A| ratio carries |R_Z| / |R_P|.
    """
    terms = {}
    for name, (alt_p, alt_z, polarity) in alternative_stages(
        pressure, vertical
    ).items():
        r_p = stage_ratio(frequency, pressure, alt_p)
        r_z = stage_ratio(frequency, vertical, alt_z)
        phase = -np.angle(r_p, deg=True) + np.angle(r_z, deg=True) + polarity
        terms[name] = (wrap(phase), np.abs(r_z) / np.abs(r_p))
    return terms


def laplace_stages(channel):
    stages = []
    for stage in channel.response.response_stages:
        poles = list(getattr(stage, "poles", None) or [])
        zeros = list(getattr(stage, "zeros", None) or [])
        if not poles and not zeros:
            continue
        kind = getattr(stage, "pz_transfer_function_type", "")
        if kind != "LAPLACE (RADIANS/SECOND)":
            if (
                unpaired_roots(poles)
                or unpaired_roots(zeros)
                or any(p.real > 0 for p in poles)
            ):
                raise ValueError("Anomalous roots outside a Laplace rad/s stage")
            continue
        stages.append(([complex(p) for p in poles], [complex(z) for z in zeros]))
    return stages


def snapshot_channel(path, seed_id, epoch=None):
    """The single channel epoch stored in a per-channel response snapshot."""
    from obspy import UTCDateTime, read_inventory

    network, station, location, code = seed_id.split(".")
    matches = [
        c
        for n in read_inventory(str(path))
        if n.code == network
        for s in n
        if s.code == station
        for c in s
        if c.code == code and c.location_code == location
    ]
    if len(matches) != 1:
        raise ValueError(f"Snapshot does not contain exactly one {seed_id} epoch")
    channel = matches[0]
    if epoch is not None:
        for key, actual in (
            ("channel_start_time", channel.start_date),
            ("channel_end_time", channel.end_date),
        ):
            expected = UTCDateTime(epoch[key]) if epoch.get(key) else None
            if expected != actual:
                raise ValueError(f"Snapshot epoch mismatch for {seed_id}: {key}")
    return channel


def coherent_event(coherence, threshold):
    coherence = np.asarray(coherence, dtype=float)
    return bool(
        np.median(coherence) >= threshold
        and np.mean(coherence) >= threshold
        and (coherence >= threshold).sum() >= 3
    )


def event_category(event):
    status = event.get("status")
    if status in EVENT_CATEGORIES:
        return EVENT_CATEGORIES[status]
    reason = (event.get("reason") or "").lower()
    if "boundary" in reason:
        return "fit_bound_unresolved"
    if status in ("unresolved", "failed", "rejected_input"):
        return "input_refused"
    raise ValueError(f"Unknown saved event outcome: {status}")


def rayleigh_metrics(events, terms, threshold):
    """Phase residual and level-free amplitude slope over coherent events."""
    out = {}
    for name in HYPOTHESES:
        model, amplitude = terms[name]
        residuals, slopes, passing = [], [], 0
        for event in events:
            residual = wrap(event["phase_deg"] - model)
            mask = event["coherence_squared"] >= threshold
            residuals.append(residual[mask])
            passing += int((np.abs(residual[mask]) <= PHASE_SCREEN_DEG).sum() >= 3)
            ratio = event["psd_ratio_pa"] / event["expected_pa"] / amplitude
            slopes.append(float(np.log(ratio[-1] / ratio[0])))
        pooled = np.concatenate(residuals) if residuals else np.array([])
        out[name] = {
            "phase_rms_deg": float(np.sqrt(np.mean(pooled**2)))
            if pooled.size
            else None,
            "phase_low_bin_deg": circular_mean([r[0] for r in residuals if r.size])
            if pooled.size
            else None,
            "phase_high_bin_deg": circular_mean([r[-1] for r in residuals if r.size])
            if pooled.size
            else None,
            "fraction_bins_within_20deg": float(
                np.mean(np.abs(pooled) <= PHASE_SCREEN_DEG)
            )
            if pooled.size
            else None,
            "events_passing_phase_screen_only": passing,
            "median_log_ratio_slope": float(np.median(slopes)) if slopes else None,
        }
    return out


def microseism_metrics(bins, pressure, vertical):
    """Same hypotheses on saved microseism bins (convention D conj(P) = -Rayleigh)."""
    f = bins["frequency_hz"]
    m = (
        (bins["coherence_squared"] >= MICROSEISM_COHERENCE)
        & (f >= MICROSEISM_BAND_HZ[0])
        & (f <= MICROSEISM_BAND_HZ[1])
    )
    if m.sum() < 10:
        return {"coherent_bins": int(m.sum())}
    terms = hypothesis_terms(f[m], pressure, vertical)
    out = {
        "coherent_bins": int(m.sum()),
        "circular_mean_phase_deg": circular_mean(bins["phase_deg"][m]),
    }
    for name in HYPOTHESES:
        residual = wrap(bins["phase_deg"][m] + terms[name][0])
        out[f"{name}_phase_rms_deg"] = float(np.sqrt(np.mean(residual**2)))
    return out


def day_files(directory, channel):
    folder = Path(directory) / f"{channel}.D"
    days = set()
    if folder.is_dir():
        for path in folder.iterdir():
            parts = path.name.split(".")
            if len(parts) >= 7 and parts[-3] == "D":
                days.add((int(parts[-2]), int(parts[-1])))
    return days


def local_coverage(data_root, station, start, end):
    from obspy import UTCDateTime

    root = Path(data_root) / station
    pressure, vertical = day_files(root, "BDH"), day_files(root, "BHZ")
    first, last = UTCDateTime(start), UTCDateTime(end)
    epoch_days = set()
    day = UTCDateTime(first.date)
    while day <= last:
        epoch_days.add((day.year, day.julday))
        day += 86400
    both = pressure & vertical & epoch_days
    return {
        "pressure_day_files": len(pressure),
        "vertical_day_files": len(vertical),
        "epoch_days": len(epoch_days),
        "paired_epoch_day_files": len(both),
        "paired_epoch_fraction": round(len(both) / len(epoch_days), 4),
        "first_paired_day": "%04d.%03d" % min(both) if both else None,
        "last_paired_day": "%04d.%03d" % max(both) if both else None,
    }


def load_event(event):
    with np.load(event["files"]["spectra"], allow_pickle=False) as data:
        arrays = {
            key: np.asarray(data[key])
            for key in (
                "frequency_hz",
                "phase_deg",
                "coherence_squared",
                "psd_ratio_pa",
                "expected_pa",
            )
        }
    f = arrays["frequency_hz"]
    if (
        f.ndim != 1
        or len(f) < 3
        or np.any(f <= 0)
        or np.any(np.diff(f) <= 0)
        or any(a.shape != f.shape or not np.isfinite(a).all() for a in arrays.values())
        or np.any(arrays["psd_ratio_pa"] <= 0)
        or np.any(arrays["expected_pa"] <= 0)
        or np.any(arrays["coherence_squared"] < 0)
        or np.any(arrays["coherence_squared"] > 1 + 1e-10)
    ):
        raise ValueError("Invalid or mismatched saved spectral arrays")
    return arrays


def station_record(
    entry, trial_root, data_root, microseism, provenance, source_hashes=None
):
    seed_id = entry["seed_id"]
    station = seed_id.split(".")[1]
    update = entry["update"]
    comparison_path = Path(update["rayleigh_result_file"])
    record_hash(comparison_path, provenance, source_hashes)
    comparison = json.loads(comparison_path.read_text())
    trial_path = Path(trial_root) / f"YV_{station}" / "trial.json"
    record_hash(trial_path, provenance)
    trial = json.loads(trial_path.read_text())
    epoch = entry["epoch"]
    if (
        not seed_id.startswith("YV.")
        or not seed_id.endswith(".BDH")
        or epoch.get("Seed_id") != seed_id
        or epoch not in comparison.get("epoch_keys", [])
        or comparison.get("station") != f"YV.{station}"
    ):
        raise ValueError(f"Baseline pressure epoch binding mismatch: {seed_id}")
    if (
        trial.get("network") != "YV"
        or trial.get("station") != station
        or Path(trial.get("baseline_result", "")).resolve() != comparison_path.resolve()
    ):
        raise ValueError(f"Trial baseline binding mismatch: {seed_id}")
    settings = comparison["settings"]
    if (
        settings.get("pressure_sign") != -1
        or settings.get("tilt_method") != "none"
        or settings.get("band_hz", settings.get("band")) != [0.03, 0.07]
    ):
        raise ValueError("Diagnosis requires the frozen no-tilt, sign -1 baseline")
    pole = trial["summaries"].get("pole_absolute", {})
    record = {
        "station": station,
        "seed_id": seed_id,
        "epoch_start": entry["epoch"]["channel_start_time"],
        "epoch_end": entry["epoch"]["channel_end_time"],
        "inventory_file": epoch["inventory_file"],
        "baseline_status": entry["status"],
        "baseline_reason": comparison.get("reason"),
        "baseline_g": update.get("rayleigh_pressure_multiplier") or None,
        "baseline_candidates": comparison.get("candidate_events", 0),
        "trial_status": pole.get("status"),
        "trial_g": pole.get("median_pressure_multiplier"),
        "trial_candidates": pole.get("candidate_events"),
        "positive_pressure_poles_rad_s": trial.get(
            "positive_real_pressure_poles_rad_s"
        ),
        "depth_m": comparison.get("settings", {}).get("depth_m"),
    }
    events = comparison.get("events", [])
    counts = {}
    for event in events:
        key = event_category(event)
        counts[key] = counts.get(key, 0) + 1
    record["event_outcomes"] = counts
    record["events_attempted"] = sum(
        v for k, v in counts.items() if k != "selection_excluded_24h"
    )
    search = comparison.get("magnitude_search", [])
    record["magnitude_search"] = [
        [
            s["minimum_magnitude"],
            s["eligible_catalog_events"],
            s["attempted_events"],
            s["outcome"],
        ]
        for s in search
    ]
    record["eligible_lowest_level"] = (
        search[-1]["eligible_catalog_events"] if search else None
    )
    if data_root:
        record["local_coverage"] = local_coverage(
            data_root, station, record["epoch_start"], record["epoch_end"]
        )
    snapshots = comparison.get("response_snapshot_sha256", {})
    if not snapshots:
        record["response"] = None
        return record
    for path, expected in snapshots.items():
        record_hash(path, provenance, source_hashes, expected)
    channels = {}
    epochs = {e["Seed_id"]: e for e in comparison["epoch_keys"]}
    for path in snapshots:
        code = Path(path).parent.name.split("__")[0]
        if code not in epochs or code in channels:
            raise ValueError(f"Unbound or duplicate response snapshot: {code}")
        channels[code] = snapshot_channel(path, code, epochs[code])
    vertical_ids = {e["tilt"]["vertical_component"] for e in events if "tilt" in e}
    if len(vertical_ids) != 1:
        raise ValueError(f"{station}: vertical partner is not unique")
    vertical_id = vertical_ids.pop()
    if not vertical_id.startswith(f"YV.{station}.") or vertical_id not in epochs:
        raise ValueError(f"Unbound vertical partner: {vertical_id}")
    record["vertical_seed_id"] = vertical_id
    pressure_stages = laplace_stages(channels[seed_id])
    vertical_stages = laplace_stages(channels[vertical_id])
    record["response"] = {
        "pressure_sensor": getattr(channels[seed_id].sensor, "description", None),
        "vertical_sensor": getattr(channels[vertical_id].sensor, "description", None),
        "vertical_unpaired_roots": [
            [r.real, r.imag]
            for p, z in vertical_stages
            for r in unpaired_roots(p) + unpaired_roots(z)
        ],
        "pressure_unpaired_roots": [
            [r.real, r.imag]
            for p, z in pressure_stages
            for r in unpaired_roots(p) + unpaired_roots(z)
        ],
    }
    threshold = settings["coherence_squared_threshold"]
    if threshold != 0.95 or settings.get("phase_limit_deg") != PHASE_SCREEN_DEG:
        raise ValueError("Unexpected saved coherence or phase screen")
    trial_events = {e["id"]: e for e in trial.get("events", [])}
    loaded, reproduction, spectra_hashes = [], [], []
    for event in events:
        if "spectra" not in event.get("files", {}):
            continue
        path = event["files"]["spectra"]
        actual = record_hash(path, provenance, source_hashes)
        spectra = load_event(event)
        spectra_hashes.append(f"{actual}  {path}")
        if not coherent_event(spectra["coherence_squared"], threshold):
            continue
        loaded.append(spectra)
        saved = trial_events.get(event["id"], {}).get("reflected_phase_deg")
        if saved is not None:
            mine = wrap(
                spectra["phase_deg"]
                - hypothesis_terms(
                    spectra["frequency_hz"], pressure_stages, vertical_stages
                )["H1"][0]
            )
            reproduction.append(float(np.max(np.abs(wrap(mine - np.asarray(saved))))))
    record["event_spectra"] = {
        "files": len(spectra_hashes),
        "sha256_of_sorted_file_hash_lines": hashlib.sha256(
            "\n".join(sorted(spectra_hashes)).encode()
        ).hexdigest(),
    }
    record["coherent_events"] = len(loaded)
    record["trial_reflection_reproduction_max_abs_deg"] = (
        max(reproduction) if reproduction else None
    )
    if loaded:
        f = loaded[0]["frequency_hz"]
        if any(not np.array_equal(e["frequency_hz"], f) for e in loaded):
            raise ValueError("Event spectra use different frequency grids")
        terms = hypothesis_terms(f, pressure_stages, vertical_stages)
        record["rayleigh"] = rayleigh_metrics(loaded, terms, threshold)
        record["nominal_ratio_level_median"] = float(
            np.median([np.median(e["psd_ratio_pa"] / e["expected_pa"]) for e in loaded])
        )
        middle = hypothesis_terms(np.array([0.05]), pressure_stages, vertical_stages)
        record["amplitude_factor_at_0p05hz"] = {
            k: float(v[1][0]) for k, v in middle.items()
        }
    if microseism and station in microseism:
        paired = microseism[station]
        if (
            paired["pressure_seed_id"] != seed_id
            or paired["vertical_seed_id"] != vertical_id
            or paired["pressure_epoch"] != epoch
            or paired["vertical_epoch"] != epochs[vertical_id]
        ):
            raise ValueError(f"Microseism P/Z identity mismatch: {station}")
        bins_path = paired["bins"]
        record_hash(bins_path, provenance)
        with Path(bins_path).open() as stream:
            rows = list(csv.DictReader(stream))
        bins = {
            k: np.array([float(r[k]) for r in rows])
            for k in ("frequency_hz", "coherence_squared", "phase_deg")
        }
        record["microseism"] = microseism_metrics(
            bins, pressure_stages, vertical_stages
        )
    return record


def microseism_bins(run, provenance):
    run = Path(run)
    record_hash(run / "summary.json", provenance)
    summary = json.loads((run / "summary.json").read_text())
    paths = {}
    for item in summary["stations"]:
        if item.get("status") != "checked":
            continue
        record_hash(item["check_summary"], provenance)
        check = json.loads(Path(item["check_summary"]).read_text())
        pair_path = item["paired_result"]
        pair_hash = check.get("inputs", {}).get(pair_path, {}).get("sha256")
        if not pair_hash:
            raise ValueError("Microseism check does not bind its paired source")
        record_hash(pair_path, provenance, expected_sha=pair_hash)
        pair = json.loads(Path(pair_path).read_text())
        station = item["station"]
        if (
            station in paths
            or item.get("seed_id") != check.get("pressure_seed_id")
            or check.get("station") != f"YV.{station}"
            or pair["pressure_epoch"]["Seed_id"] != check.get("pressure_seed_id")
            or pair["vertical_epoch"]["Seed_id"] != check.get("vertical_seed_id")
        ):
            raise ValueError(f"Duplicate or mismatched microseism station: {station}")
        paths[station] = {
            "bins": Path(check["files"]["bins"]),
            "pressure_seed_id": check["pressure_seed_id"],
            "vertical_seed_id": check["vertical_seed_id"],
            "pressure_epoch": pair["pressure_epoch"],
            "vertical_epoch": pair["vertical_epoch"],
        }
    return paths


CSV_FIELDS = (
    "station",
    "seed_id",
    "epoch_start",
    "epoch_end",
    "inventory_file",
    "vertical_seed_id",
    "baseline_status",
    "baseline_g",
    "baseline_candidates",
    "trial_g",
    "trial_status",
    "trial_candidates",
    "positive_pressure_poles_rad_s",
    "vertical_unpaired_roots",
    "events_attempted",
    "selection_excluded_24h",
    "input_refused",
    "coherence_failed",
    "phase_review",
    "arrival_overlap_review",
    "fit_bound_unresolved",
    "candidate",
    "coherent_events",
    "eligible_lowest_level",
    "paired_epoch_fraction",
    *(f"rayleigh_{h}_phase_rms_deg" for h in HYPOTHESES),
    *(f"rayleigh_{h}_events_passing_phase_screen_only" for h in HYPOTHESES),
    *(f"rayleigh_{h}_median_log_ratio_slope" for h in HYPOTHESES),
    *(f"microseism_{h}_phase_rms_deg" for h in HYPOTHESES),
)


def csv_row(record):
    row = {k: record.get(k) for k in CSV_FIELDS if k in record}
    row.update(record.get("event_outcomes", {}))
    row["vertical_unpaired_roots"] = len(
        (record.get("response") or {}).get("vertical_unpaired_roots", [])
    )
    row["positive_pressure_poles_rad_s"] = ";".join(
        f"{p:g}" for p in record.get("positive_pressure_poles_rad_s") or []
    )
    row["paired_epoch_fraction"] = (record.get("local_coverage") or {}).get(
        "paired_epoch_fraction"
    )
    for h, values in (record.get("rayleigh") or {}).items():
        for key in (
            "phase_rms_deg",
            "events_passing_phase_screen_only",
            "median_log_ratio_slope",
        ):
            value = values[key]
            row[f"rayleigh_{h}_{key}"] = (
                round(value, 3) if isinstance(value, float) else value
            )
    for key, value in (record.get("microseism") or {}).items():
        if key.endswith("_phase_rms_deg") and key[:2] in HYPOTHESES:
            row[f"microseism_{key}"] = round(value, 3)
    return {k: ("" if row.get(k) is None else row.get(k)) for k in CSV_FIELDS}


def run(args):
    output = Path(args.output_dir)
    if output.exists():
        raise SystemExit(f"Refusing to reuse existing output folder: {output}")
    workflow_path = Path(args.workflow)
    trial_root = Path(args.trial_root)
    experiment = trial_root / "experiment.json"
    provenance = {
        str(workflow_path.resolve()): sha256(workflow_path),
        str(experiment.resolve()): sha256(experiment),
    }
    if (
        args.expect_trial_sha
        and provenance[str(experiment.resolve())] != args.expect_trial_sha
    ):
        raise SystemExit("Trial experiment.json hash differs from the expected value")
    trial_manifest = json.loads(experiment.read_text())
    source_hashes = trial_manifest.get("source_sha256")
    if (
        trial_manifest.get("status") != "complete"
        or trial_manifest.get("protected_inputs_unchanged") is not True
        or not isinstance(source_hashes, dict)
    ):
        raise ValueError("A completed, preserved trial source manifest is required")
    record_hash(workflow_path, provenance, source_hashes)
    workflow = json.loads(workflow_path.read_text())
    if workflow.get("execution_status") != "complete":
        raise ValueError("The frozen workflow must be complete")
    microseism = None
    if args.microseism_run:
        summary = Path(args.microseism_run) / "summary.json"
        record_hash(summary, provenance)
        microseism = microseism_bins(args.microseism_run, provenance)
    records = [
        station_record(
            entry, trial_root, args.data_root, microseism, provenance, source_hashes
        )
        for entry in sorted(workflow["entries"].values(), key=lambda e: e["seed_id"])
    ]
    import obspy

    if len({r["seed_id"] for r in records}) != len(records):
        raise ValueError("Duplicate pressure identities in this single-epoch YV audit")
    verify_inputs_unchanged(provenance)
    result = {
        "status": "diagnostic_only",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "hypotheses": HYPOTHESIS_LABELS,
        "workflow_options": workflow.get("options"),
        "phase_screen_deg": PHASE_SCREEN_DEG,
        "microseism_band_hz": MICROSEISM_BAND_HZ,
        "microseism_coherence_squared": MICROSEISM_COHERENCE,
        "gain_fitted": False,
        "gain_applied": False,
        "responses_or_masters_changed": False,
        "waveforms_read": False,
        "protected_inputs_unchanged": True,
        "code_sha256": sha256(__file__),
        "versions": {"numpy": np.__version__, "obspy": obspy.__version__},
        "stations": records,
        "input_sha256": dict(sorted(provenance.items())),
    }
    output.mkdir(parents=True)
    with (output / "evidence.json").open("x") as stream:
        json.dump(result, stream, indent=1, default=float, allow_nan=False)
        stream.write("\n")
    with (output / "stations.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(csv_row(r) for r in records)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--workflow", required=True)
    parser.add_argument("--trial-root", required=True)
    parser.add_argument("--expect-trial-sha")
    parser.add_argument("--microseism-run")
    parser.add_argument("--data-root", help="YV raw folder; only file names are listed")
    parser.add_argument("--output-dir", required=True)
    result = run(parser.parse_args())
    stations, inputs = len(result["stations"]), len(result["input_sha256"])
    print(f"{stations} stations; {inputs} hashed inputs")


if __name__ == "__main__":
    main()
