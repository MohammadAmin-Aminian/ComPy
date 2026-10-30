"""Export local-P campaign results keyed to master channel epochs, and compare.

Read-only for the campaign, the master and the saved Rayleigh/microseism
tables. Writes into a fresh results/LOCAL_P/ folder:

- local_p_master_export.csv: one row per master pressure epoch (Seed_id,
  channel_start_time, channel_end_time, inventory_file) with local_p_* fields.
  Numeric results, execution status and scientific outcome are separate
  columns. No score or check-status column is written: scoring is not agreed.
- local_p_comparison.csv: local-P band multipliers and polarity beside the
  saved Rayleigh multiplier (0.03-0.07 Hz) and microseism phase family. Bands
  never match and channels may differ; ratios are context, not validation.
"""

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

import numpy as np
from obspy import read_inventory

from compy_pressure.local_p_check import PROTOCOL_VERSION
from compy_pressure.local_p_io import local_p_output

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

ROOT = Path(__file__).resolve().parents[1]
MASTER = Path.home() / "Projects/bruit-fm/exports/master_inventory.csv"
RAYLEIGH = ROOT / "docs/calibration_network_summary.json"
MICROSEISM = ROOT / "docs/microseism_baseline_comparison.csv"
INVENTORY_MAP = ROOT / "docs/local_p_inventory_map.json"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def master_rows(path, seed_ids):
    with Path(path).open(newline="") as stream:
        return [r for r in csv.DictReader(stream) if r["Seed_id"] in seed_ids]


def rayleigh_table(path):
    data = json.loads(Path(path).read_text())
    table = {}
    for network in data["networks"]:
        for item in network["stations"]:
            table[network["network"], item["station"]] = {
                "rayleigh_seed_id": item["seed_id"],
                "rayleigh_status": item["status"],
                "rayleigh_pressure_multiplier": item["pressure_multiplier"],
                "rayleigh_band_hz": network["options"]["band"],
                "rayleigh_assumed_pressure_sign": network["options"]["pressure_sign"],
            }
    return table


def microseism_table(path):
    table = {}
    with Path(path).open(newline="") as stream:
        for row in csv.DictReader(stream):
            table[row["network"], row["station"].split(".")[-1]] = {
                "microseism_pressure_seed_id": row["pressure_seed_id"],
                "microseism_phase_family": row["phase_family"],
                "microseism_phase_deg": row["median_nominal_transfer_phase_deg"],
            }
    return table


def band_json(summary):
    return {
        f"{b['band_hz'][0]:g}-{b['band_hz'][1]:g}": {
            "events": b["events"],
            "gain_median": b.get("gain_median"),
            "gain_mad": b.get("gain_mad"),
            "noise_bracket": [
                b.get("gain_h2_lower_median"),
                b.get("gain_h1_upper_median"),
            ],
            "model_interval": [
                b.get("gain_model_low_median"),
                b.get("gain_model_high_median"),
            ],
            "regimes": b.get("regimes"),
        }
        for b in summary.get("bands", [])
        if b["events"]
    }


def result_status(summary, row):
    status = summary.get("outcome") or row.get("plan_status") or ""
    if row.get("execution_status") == "failed":
        return "failed"
    if (
        status == "conditional_gain_candidate"
        and summary.get("protocol_version") != PROTOCOL_VERSION
    ):
        return "legacy_protocol_unvalidated"
    return status


@lru_cache(maxsize=32)
def inventory_channels(path):
    return read_inventory(str(path))


def channel_identity(path, seed, epoch):
    net, sta, loc, code = seed.split(".")
    found = [
        c
        for n in inventory_channels(str(path)).select(
            network=net, station=sta, location=loc, channel=code
        )
        for s in n
        for c in s
        if [str(c.start_date or ""), str(c.end_date or "")] == list(epoch)
    ]
    if len(found) != 1:
        raise ValueError(f"No unique channel epoch in {path}: {seed} {epoch}")
    c = found[0]
    sensitivity = c.response.instrument_sensitivity if c.response else None
    return (
        float(c.sample_rate),
        float(c.latitude),
        float(c.longitude),
        float(c.elevation),
    ), (sensitivity.input_units.upper().replace(" ", "") if sensitivity else "")


def match_master_epochs(summary, row, candidates, mappings=()):
    """Resolve full keys, using only explicit hash-bound inventory aliases."""
    source = row.get("inventory") or next(
        (p for p in summary.get("inputs", {}) if Path(p).suffix.lower() == ".xml"), None
    )
    if not source:
        return [], "source_inventory_not_recorded"
    source = Path(source).resolve()
    recorded = summary.get("inputs", {}).get(str(source), {}).get("sha256")
    source_hash = sha256(source)
    if recorded and source_hash != recorded:
        raise ValueError(f"Source inventory changed: {source}")
    epochs = {tuple(e) for e in summary.get("pressure_epochs_in_coverage", [])}
    matched = []
    for m in candidates:
        epoch = (m["channel_start_time"], m["channel_end_time"])
        if epochs and epoch not in epochs:
            continue
        target = Path(m["inventory_file"]).resolve()
        direct = target == source
        aliases = [
            a
            for a in mappings
            if Path(a["source_inventory"]).resolve() == source
            and Path(a["master_inventory"]).resolve() == target
            and row["pressure_seed_id"] in a["seed_ids"]
        ]
        if not direct:
            if len(aliases) != 1:
                continue
            alias = aliases[0]
            if (
                source_hash != alias["source_sha256"]
                or sha256(target) != alias["master_sha256"]
            ):
                raise ValueError("Inventory alias hash mismatch")
        source_meta, source_units = channel_identity(
            source, row["pressure_seed_id"], epoch
        )
        target_meta, target_units = channel_identity(
            target, row["pressure_seed_id"], epoch
        )
        if (
            not np.allclose(source_meta, target_meta, rtol=0, atol=1e-6)
            or source_units != target_units
        ):
            raise ValueError("Mapped channel epoch metadata differ")
        matched.append(m)
    keys = [
        tuple(
            m[k]
            for k in (
                "Seed_id",
                "channel_start_time",
                "channel_end_time",
                "inventory_file",
            )
        )
        for m in matched
    ]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate full master epoch key")
    dates = [(m["channel_start_time"], m["channel_end_time"]) for m in matched]
    if len(dates) != len(set(dates)):
        raise ValueError("Ambiguous master inventory for the same channel epoch")
    return (
        matched,
        "verified_inventory_mapping" if matched else "no_verified_inventory_epoch",
    )


def epoch_summary(summary, epoch):
    """Never copy a pooled gain to an epoch without accepted-event evidence."""
    epoch = list(epoch)
    for part in summary.get("epoch_results", []):
        if part["pressure_epoch"] == epoch:
            return {
                **summary,
                **part,
                "event_reason_counts": {
                    "accepted": len(
                        {
                            event
                            for b in part.get("bands", [])
                            for event in b.get("event_ids", [])
                        }
                    )
                },
            }
    accepted = summary.get("pressure_epochs_accepted", [])
    if len(accepted) == 1 and list(accepted[0]) == epoch:
        return summary
    if not accepted and len(summary.get("pressure_epochs_in_coverage", [])) == 1:
        return summary  # no numerical results to misattribute
    return {
        **summary,
        "outcome": "epoch_attribution_unresolved",
        "bands": [],
        "frequency_dependence": {},
        "polarity_majority": None,
        "regime_counts": {},
        "event_reason_counts": {},
    }


def flags(summary, row):
    out = []
    if summary and summary.get("protocol_version") != PROTOCOL_VERSION:
        out.append("legacy_noise_guard_unverified")
        if summary.get("polarity_majority") is not None:
            out.append("legacy_polarity_summary_unverified")
    if summary.get("polarity_summary_basis"):
        out.append("polarity_dominant_classified_sign_not_confidence")
    regimes = summary.get("regime_counts") or {}
    if regimes.get("pre_multiple"):
        out.append("pre_multiple_diagnostic_only")
        out.append("pre_multiple_filter_leakage_not_field_validated")
    # Older campaigns recorded only lag_search_s; newer ones record timing_mode.
    mode = summary.get("timing_mode")
    if mode in {"search", "drift_model"} or (
        mode is None and float(row.get("lag_search_s") or 0) > 0
    ):
        out.append("timing_lag_estimated_phase_not_independent")
    if "timing_unverified" in (row.get("note") or ""):
        out.append("separate_logger_timing_unverified")
    if len(summary.get("pressure_epochs_accepted", [])) > 1:
        out.append("multiple_epochs_pooled")
    if summary.get("location_code_mismatch_in_filenames"):
        out.append("filename_location_differs_from_header")
    return out


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--master", type=Path, default=MASTER)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--inventory-map", type=Path, default=INVENTORY_MAP)
    args = parser.parse_args()
    out = local_p_output(args.output_dir)
    if out.exists() and any(out.iterdir()):
        raise SystemExit("Output must be a new or empty folder under results/LOCAL_P/")
    campaign = json.loads((args.campaign_dir / "campaign_summary.json").read_text())
    protected = [
        args.master,
        RAYLEIGH,
        MICROSEISM,
        args.campaign_dir / "campaign_summary.json",
    ]
    before = {str(p): sha256(p) for p in protected}
    mappings = (
        json.loads(args.inventory_map.read_text())
        if args.inventory_map.is_file()
        else []
    )
    if args.inventory_map.is_file():
        before[str(args.inventory_map)] = sha256(args.inventory_map)
    rows = campaign["rows"]
    master = master_rows(args.master, {r["pressure_seed_id"] for r in rows})
    rayleigh, microseism = rayleigh_table(RAYLEIGH), microseism_table(MICROSEISM)
    export, comparison = [], []
    for row in rows:
        summary = {}
        if row.get("result_dir"):
            path = Path(row["result_dir"]) / "summary.json"
            summary = json.loads(path.read_text()) if path.is_file() else {}
        status = result_status(summary, row)
        epochs = {tuple(e) for e in summary.get("pressure_epochs_in_coverage", [])}
        candidates = [m for m in master if m["Seed_id"] == row["pressure_seed_id"]]
        matched, match_note = match_master_epochs(summary, row, candidates, mappings)
        if summary:
            before[str(path)] = sha256(path)
            for p, state in summary.get("inputs", {}).items():
                if sha256(p) != state["sha256"]:
                    raise ValueError(f"Saved input changed: {p}")
                before[p] = state["sha256"]
        for m in matched:
            before[m["inventory_file"]] = sha256(m["inventory_file"])
        bands = band_json(summary)
        base = {
            "local_p_status": status,
            "local_p_protocol_version": summary.get(
                "protocol_version", "legacy_v1" if summary else ""
            ),
            "local_p_validation_status": "unvalidated",
            "local_p_inventory_resolution": match_note,
            "local_p_execution_status": row.get("execution_status") or "not_run",
            "local_p_vertical_seed_id": row.get("vertical_seed_id", ""),
            "local_p_regimes": json.dumps(summary.get("regime_counts", {})),
            "local_p_accepted_events": row.get("accepted_events", ""),
            "local_p_gain_by_band": json.dumps(bands) if bands else "",
            "local_p_frequency_assessment": (
                summary.get("frequency_dependence") or {}
            ).get("assessment", ""),
            "local_p_scalar_gain_if_consistent": (
                summary.get("frequency_dependence") or {}
            ).get("scalar_median_all_bands", ""),
            "local_p_polarity_majority": summary.get("polarity_majority", ""),
            "local_p_water_depth_m": summary.get("water_depth_m", ""),
            "local_p_result_file": str(Path(row["result_dir"]) / "summary.json")
            if row.get("result_dir")
            else "",
            "local_p_plot_file": summary.get("figure", ""),
            "local_p_review_flags": ";".join(flags(summary, row)),
            "local_p_error": summary.get("error", "") or row.get("error", ""),
            "local_p_note": row.get("note", ""),
            "local_p_score_note": "no score assigned; scoring policy not agreed",
        }
        if not matched:
            export.append(
                {
                    "Seed_id": row["pressure_seed_id"],
                    "channel_start_time": "",
                    "channel_end_time": "",
                    "inventory_file": "",
                    "master_epoch_match": "no_master_epoch",
                    **base,
                }
            )
        for m in matched:
            scoped = (
                epoch_summary(summary, (m["channel_start_time"], m["channel_end_time"]))
                if summary
                else {}
            )
            scoped_bands = band_json(scoped)
            scoped_base = {
                **base,
                "local_p_status": result_status(scoped, row),
                "local_p_gain_by_band": json.dumps(scoped_bands)
                if scoped_bands
                else "",
                "local_p_accepted_events": scoped.get("event_reason_counts", {}).get(
                    "accepted", 0
                )
                if summary
                else "",
                "local_p_regimes": json.dumps(scoped.get("regime_counts", {})),
                "local_p_frequency_assessment": scoped.get(
                    "frequency_dependence", {}
                ).get("assessment", ""),
                "local_p_scalar_gain_if_consistent": scoped.get(
                    "frequency_dependence", {}
                ).get("scalar_median_all_bands", ""),
                "local_p_polarity_majority": scoped.get("polarity_majority", ""),
            }
            export.append(
                {
                    "Seed_id": m["Seed_id"],
                    "channel_start_time": m["channel_start_time"],
                    "channel_end_time": m["channel_end_time"],
                    "inventory_file": m["inventory_file"],
                    "master_epoch_match": "exact" if epochs else "channel_not_run",
                    **scoped_base,
                }
            )
        if not bands:
            continue
        key = (row["network"], row["station"])
        ray = rayleigh.get(key, {})
        micro = microseism.get(key, {})
        item = {
            "network": row["network"],
            "station": row["station"],
            "local_p_pressure_seed_id": row["pressure_seed_id"],
            "local_p_status": status,
            "local_p_polarity_majority": summary.get("polarity_majority"),
            **{f"local_p_{k}_gain": v["gain_median"] for k, v in bands.items()},
            **{f"local_p_{k}_events": v["events"] for k, v in bands.items()},
            **ray,
            **micro,
            "same_pressure_channel_as_rayleigh": ray.get("rayleigh_seed_id")
            == row["pressure_seed_id"],
            "band_match": False,
            "comparison_note": (
                "different bands (and possibly channels); context only, not validation"
            ),
        }
        if ray.get("rayleigh_pressure_multiplier"):
            for k, v in bands.items():
                item[f"ratio_local_p_{k}_to_rayleigh"] = (
                    v["gain_median"] / ray["rayleigh_pressure_multiplier"]
                )
        polarity = summary.get("polarity_majority")
        if polarity and ray:
            item["polarity_agrees_with_rayleigh_assumed_sign"] = (
                polarity == ray["rayleigh_assumed_pressure_sign"]
            )
        comparison.append(item)
    out.mkdir(parents=True, exist_ok=True)
    for name, table in (
        ("local_p_master_export.csv", export),
        ("local_p_comparison.csv", comparison),
    ):
        fields = []
        for item in table:
            fields.extend(k for k in item if k not in fields)
        with (out / name).open("x", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=fields or ["empty"], lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(table)
    after = {p: sha256(p) for p in before}
    report = {
        "created_utc": datetime.now(UTC).isoformat(),
        "campaign": str(args.campaign_dir.resolve()),
        "master": str(args.master),
        "rows_exported": len(export),
        "rows_compared": len(comparison),
        "status_counts": {},
        "inputs_unchanged": before == after,
        "inputs_sha256": before,
    }
    for item in export:
        report["status_counts"][item["local_p_status"]] = (
            report["status_counts"].get(item["local_p_status"], 0) + 1
        )
    with (out / "export_report.json").open("x") as stream:
        stream.write(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps({k: v for k, v in report.items() if k != "inputs_sha256"}, indent=2)
    )
    if before != after:
        raise SystemExit("Protected inputs changed during export")


if __name__ == "__main__":
    main()
