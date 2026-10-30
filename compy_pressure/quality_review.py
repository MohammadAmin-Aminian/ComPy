"""Read-only, exact-epoch six-check evidence integration and approval planning.

No new scientific threshold, gain application or production master write.
Calculation outcomes are not pass/fail decisions. Approval requires a separately
adopted policy and explicit, fingerprint-bound reviewer decisions.
"""

import argparse
import copy
import csv
import hashlib
import json
import math
import re
from collections import Counter
from datetime import datetime
from pathlib import Path

from compy_pressure.inventory import RESULT_COLUMNS

ROOT = Path(__file__).resolve().parents[1] / "pressure"
KEYS = ("Seed_id", "channel_start_time", "channel_end_time", "inventory_file")
CHECKS = (
    "expected_pressure",
    "nearby",
    "infragravity",
    "rayleigh",
    "microseism",
    "local_p",
)
PREFIXES = ("pressure_expected", "nearby", "infragravity")
POLICY_VERSION = "six-check-review-v1"
VERDICTS = {"pass", "fail", "partial", "unresolved", "not_tested", "not_applicable"}


def digest(path):
    with Path(path).open("rb") as stream:
        # hashlib.file_digest needs Python 3.11; ComPy supports 3.10
        checksum = hashlib.sha256()
        for block in iter(lambda: stream.read(1 << 20), b""):
            checksum.update(block)
        return checksum.hexdigest()


def json_hash(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, allow_nan=False).encode()
    ).hexdigest()


def full_key(row):
    key = tuple(str(row.get(k, "")) for k in KEYS)
    parts = key[0].split(".")
    if (
        len(parts) != 4
        or not key[1]
        or not key[3]
        or any(not re.fullmatch(r"[A-Za-z0-9_-]+", parts[i]) for i in (0, 1, 3))
        or not re.fullmatch(r"[A-Za-z0-9_-]*", parts[2])
    ):
        raise ValueError("A complete SEED/epoch/source identity is required")
    return key


def resolve(path, root=ROOT):
    return (Path(root) / Path(path).expanduser()).resolve()


def fresh_output(path, root=ROOT):
    path = resolve(path, root)
    allowed = resolve("results/REVIEWS", root)
    if allowed not in path.parents:
        raise ValueError("Derived review outputs must be below results/REVIEWS/")
    if path.exists():
        raise FileExistsError(path)
    return path


class Inputs:
    def __init__(self, root=ROOT):
        self.root = Path(root)
        self.hashes = {}

    def protect(self, path, expected=None):
        path = resolve(path, self.root)
        current = digest(path)
        if (expected and current != expected) or (
            str(path) in self.hashes and self.hashes[str(path)] != current
        ):
            raise ValueError(f"Input hash changed: {path}")
        self.hashes[str(path)] = current
        return path

    def document(self, path, expected=None):
        return json.loads(self.protect(path, expected).read_text())

    def table(self, path):
        with self.protect(path).open(newline="") as stream:
            return list(csv.DictReader(stream))

    def verify(self):
        if any(digest(p) != value for p, value in self.hashes.items()):
            raise RuntimeError("Protected review input changed")


def evidence_state(status, flags=""):
    status = status or "not_tested"
    if status in {"failed", "processing_failed", "processing_error"}:
        return "processing_error"
    if (
        "legacy_noise_guard_unverified" in flags
        or status == "legacy_protocol_unvalidated"
    ):
        return "legacy_unvalidated"
    if status in {"conditional_repeatable_unvalidated", "conditional_gain_candidate"}:
        return "conditional"
    if status == "candidate_comparison_under_assumed_sign":
        return "provisional"
    if status in {"comparison_only", "exploratory_only", "diagnostic_only"}:
        return "diagnostic"
    if status in {
        "not_applicable",
        "not_applicable_sample_rate",
        "reference_band_unusable",
    }:
        return "not_applicable"
    if status in {"not_tested", "not_run", "no_catalog"}:
        return "not_tested"
    return "unresolved"


def evidence(
    status="not_tested", *, flags="", metrics=None, sources=None, binding="none"
):
    return {
        "state": evidence_state(status, flags),
        "source_status": status,
        "flags": flags,
        "metrics": metrics or {},
        "sources": sources or [],
        "binding": binding,
        "score": None,
        "verdict": "unresolved",
    }


def source_record(inputs, path, expected=None):
    path = inputs.protect(path, expected)
    return {"path": str(path), "sha256": inputs.hashes[str(path)]}


def pressure(row):
    return (
        str(row.get("input_units", "")).upper().replace(" ", "")
        in {"PA", "PASCAL", "PASCALS"}
        or row.get("ppsd_units") == "Pa^2/Hz"
    )


def load_masters(inputs, paths, networks=None):
    records = {}
    for name in paths:
        table = inputs.table(name)
        sibling = resolve(name, inputs.root).with_suffix(".xlsx")
        if sibling.is_file():
            inputs.protect(sibling)
        seen = set()
        for row in table:
            if not pressure(row) or (
                networks is not None and row["Seed_id"].split(".")[0] not in networks
            ):
                continue
            key = full_key(row)
            if key in seen:
                raise ValueError(f"Duplicate master full key: {key}")
            seen.add(key)
            item = records.setdefault(key, {"master_rows": [], "conflicts": [], **row})
            item["master_rows"].append(
                {"path": str(resolve(name, inputs.root)), "row": row}
            )
            for field in ("ppsd_npz_file", *[p + "_result_file" for p in PREFIXES]):
                if item.get(field) and row.get(field) and item[field] != row[field]:
                    item["conflicts"].append(field)
                elif row.get(field):
                    item[field] = row[field]
    return records


def base_review(key, record):
    metadata = {
        k: record.get(k, "")
        for k in ("Lat", "Lon", "Profondeur", "sample_rate_hz", "input_units")
    }
    return {
        **dict(zip(KEYS, key, strict=True)),
        "row_id": json_hash(key)[:24],
        "station": ".".join(key[0].split(".")[:2]),
        "metadata": metadata,
        "master_sources": [m["path"] for m in record.get("master_rows", [])],
        "checks": {name: evidence() for name in CHECKS},
        "ppsd": None,
        "issues": list(set(record.get("conflicts", []))),
        "existing_scores": [
            {
                "master": m["path"],
                "scores": {c: m["row"].get(c, "") for c in RESULT_COLUMNS},
            }
            for m in record.get("master_rows", [])
        ],
        "classification": "unresolved",
        "accepted_for_noise_model": False,
    }


def pressure_reference_warnings(summary):
    """Preserve old statistics, but never present their references as current."""
    references = summary.get("references", {})
    if not references:
        return []
    warnings = []
    if summary.get("reference_coverage", {}).get("webb_instructions", 0):
        if not references.get("webb", {}).get("verified_published_bounds", False):
            warnings.append("historical_unverified_handout_included")
    if (
        references.get("brown", {}).get("data_file")
        != "brown_hydroacoustic_fullband.csv"
    ):
        warnings.append("historical_brown_reference_support")
    return warnings


def add_master_evidence(inputs, item, record):
    for check, prefix in zip(CHECKS[:3], PREFIXES, strict=True):
        candidates = [
            (m, m["row"].get(prefix + "_result_file"))
            for m in record.get("master_rows", [])
        ]
        candidates = [(m, p) for m, p in candidates if p]
        if not candidates:
            continue
        if prefix + "_result_file" in item["issues"]:
            item["checks"][check] = evidence("ambiguous_master_evidence")
            continue
        master, path = candidates[0]
        row = master["row"]
        if not resolve(path, inputs.root).is_file():
            item["checks"][check] = evidence(
                "missing_saved_evidence", sources=[{"path": path}]
            )
            continue
        data = inputs.document(path)
        if data.get("seed_id") != item["Seed_id"]:
            raise ValueError(f"Scientific result SEED mismatch: {path}")
        if data.get("source_npz_sha256") and row.get("ppsd_npz_file"):
            inputs.protect(row["ppsd_npz_file"], data["source_npz_sha256"])
        status = data.get("status", row.get(prefix + "_status", "not_tested"))
        metrics = {
            k: v
            for k, v in row.items()
            if k.startswith(prefix + "_") and not k.endswith(("_file", "_status"))
        }
        # Original JSON keeps limits, reference caveats and compared frequency bands.
        metrics["saved_summary"] = {
            k: v for k, v in data.items() if k not in {"files", "events"}
        }
        item["checks"][check] = evidence(
            status,
            flags=";".join(pressure_reference_warnings(data))
            if check == "expected_pressure"
            else "",
            metrics=metrics,
            sources=[
                source_record(inputs, path),
                source_record(inputs, master["path"]),
            ],
            binding="exact_master_key_pointer",
        )
    p = record.get("ppsd_npz_file")
    if p and record.get("ppsd_status") in {"complete", "partial"}:
        if "ppsd_npz_file" in item["issues"]:
            return
        if not resolve(p, inputs.root).is_file():
            item["issues"].append("missing_ppsd_file")
            return
        item["ppsd"] = {
            **source_record(inputs, p),
            "status": record["ppsd_status"],
            "units": record.get("ppsd_units", ""),
            "gain_applied": False,
        }


def attach_staging(inputs, items, records, path, kind, networks, unjoined):
    seen = set()
    for row in inputs.table(path):
        if row["Seed_id"].split(".")[0] not in networks:
            continue
        if kind == "local_p" and row.get("master_epoch_match") != "exact":
            unjoined.append(
                {
                    "check": kind,
                    "seed_id": row["Seed_id"],
                    "reason": row.get("master_epoch_match", "missing_full_key"),
                }
            )
            continue
        key = full_key(row)
        if key in seen:
            raise ValueError(f"Duplicate {kind} staging full key")
        seen.add(key)
        item = items.setdefault(key, base_review(key, records.get(key, {})))
        if kind == "microseism":
            status = row[RESULT_COLUMNS[4].removesuffix("_score") + "_status"]
            flags = ";".join(
                str(row.get(k, ""))
                for k in (
                    "microseism_tilt_status",
                    "microseism_absolute_vertical_response_status",
                )
            )
        else:
            status = row["local_p_status"]
            flags = row.get("local_p_review_flags", "")
        item["checks"][kind] = evidence(
            status,
            flags=flags,
            metrics={k: v for k, v in row.items() if k not in KEYS},
            sources=[source_record(inputs, path)],
            binding="exact_staging_full_key",
        )


def attach_rayleigh(inputs, items, records, path, networks, unjoined):
    document = inputs.document(path)
    for network in document["networks"]:
        if network["network"] not in networks:
            continue
        for row in network["stations"]:
            result = resolve(row["result_file"], inputs.root)
            if not result.is_file():
                unjoined.append(
                    {
                        "check": "rayleigh",
                        "seed_id": row["seed_id"],
                        "reason": "missing_result_file",
                    }
                )
                continue
            data = inputs.document(result, row["result_sha256"])
            candidates = [
                k for k in data.get("epoch_keys", []) if k["Seed_id"] == row["seed_id"]
            ]
            if len(candidates) != 1:
                unjoined.append(
                    {
                        "check": "rayleigh",
                        "seed_id": row["seed_id"],
                        "reason": "no_unique_pressure_epoch_in_saved_result",
                    }
                )
                continue
            key = full_key(candidates[0])
            item = items.setdefault(key, base_review(key, records.get(key, {})))
            if item["checks"]["rayleigh"]["binding"] != "none":
                raise ValueError("Multiple baseline results for one full pressure key")
            item["checks"]["rayleigh"] = evidence(
                row["status"],
                flags=json.dumps(data.get("review_flags", []), sort_keys=True),
                metrics={
                    **row,
                    "band_hz": network["options"]["band"],
                    "limitations": data.get("limitations", []),
                    "assumed_pressure_sign": data.get("pressure_sign_assumption"),
                },
                sources=[source_record(inputs, path), source_record(inputs, result)],
                binding="exact_saved_result_epoch_key",
            )


def aware_time(value):
    dt = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if dt.tzinfo is None:
        raise ValueError("Reviewer timestamps must include a time zone")
    return dt


def validate_policy(policy):
    if (
        policy.get("version") != POLICY_VERSION
        or policy.get("classification_rule") != "all_reviewed_checks_v1"
        or policy.get("status") not in {"draft", "adopted"}
        or policy.get("production_master_updates_allowed") is not False
        or policy.get("gain_application_allowed") is not False
        or policy.get("noise_model_requires_explicit_channel_approval") is not True
        or policy.get("noise_model_requires_independent_reference_record") is not True
    ):
        raise ValueError("Unknown or unsafe review policy")
    if policy["status"] == "adopted":
        if not policy.get("approved_by") or not policy.get("rationale"):
            raise ValueError("An adopted policy requires recorded scientific approval")
        aware_time(policy["approved_utc"])


def classify(verdicts):
    performed = [v for v in verdicts if v not in {"not_tested", "not_applicable"}]
    if not performed:
        return "not_assessed"
    if any(v in {"partial", "unresolved"} for v in performed):
        return "unresolved"
    if all(v == "pass" for v in performed):
        return "valid_under_reviewed_checks"
    if len(performed) >= 2 and all(v == "fail" for v in performed):
        return "invalid_under_reviewed_checks"
    return "unresolved"


def apply_decisions(items, decisions, policy):
    """Apply decisions to the derived review only, never to original masters."""
    validate_policy(policy)
    if decisions.get("policy_version") != policy["version"]:
        raise ValueError("Decision policy version mismatch")
    entries = decisions.get("decisions", [])
    if entries and policy["status"] != "adopted":
        raise ValueError("Scientific policy is draft; no decisions may be adopted")
    by_id = {item["row_id"]: item for item in items}
    seen = set()
    for decision in entries:
        identifier = decision["row_id"]
        if identifier in seen or identifier not in by_id:
            raise ValueError("Duplicate or unknown decision row_id")
        seen.add(identifier)
        item = by_id[identifier]
        if decision.get("evidence_fingerprint") != item["evidence_fingerprint"]:
            raise ValueError("Stale scientific decision: evidence fingerprint mismatch")
        if not decision.get("reviewer") or not decision.get("rationale"):
            raise ValueError("Reviewer and rationale are required")
        aware_time(decision["reviewed_utc"])
        checks = decision["checks"]
        if set(checks) != set(CHECKS):
            raise ValueError("Explicit decisions for all six checks are required")
        for name, entry in checks.items():
            verdict, score = entry["verdict"], entry.get("score")
            if verdict not in VERDICTS or not entry.get("rationale"):
                raise ValueError("Unknown verdict or missing check rationale")
            numeric = isinstance(score, (int, float)) and not isinstance(score, bool)
            expected = {"pass": 1, "fail": 0}
            if verdict in expected:
                if not numeric or score != expected[verdict]:
                    raise ValueError("Pass/fail scores must be 1/0")
            elif verdict == "partial":
                if not numeric or not math.isfinite(score) or not 0 < score < 1:
                    raise ValueError(
                        "A partial score must lie strictly between 0 and 1"
                    )
            elif score is not None:
                raise ValueError(
                    "Unperformed/inapplicable/unresolved scores stay blank"
                )
            if verdict in {"pass", "fail", "partial"} and (
                not item["checks"][name]["sources"]
                or item["checks"][name]["state"] in {"not_tested", "processing_error"}
            ):
                raise ValueError(
                    "A processing failure or absent evidence cannot be scored"
                )
        item["classification"] = classify([checks[c]["verdict"] for c in CHECKS])
        item["review_decision"] = decision
        for name, entry in checks.items():
            item["checks"][name].update(
                verdict=entry["verdict"], score=entry.get("score")
            )
        if decision.get("noise_model_approved"):
            reference = decision.get("independent_reference", {})
            band = decision.get("approved_frequency_band_hz", [])
            period = decision.get("approved_time_range_utc", [])
            if (
                item["classification"] != "valid_under_reviewed_checks"
                or decision.get("nominal_response_accepted") is not True
                or not item["ppsd"]
                or item["issues"]
                or not reference.get("citation")
                or not reference.get("description")
                or len(band) != 2
                or any(
                    isinstance(x, bool)
                    or not isinstance(x, (int, float))
                    or not math.isfinite(x)
                    for x in band
                )
                or not 0 < band[0] < band[1]
                or len(period) != 2
                or aware_time(period[0]) >= aware_time(period[1])
            ):
                raise ValueError(
                    "Noise approval needs valid checks, PPSD, reference, band and times"
                )
            item["accepted_for_noise_model"] = True


def write_csv(path, rows, fields):
    with Path(path).open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def station_summary(rows):
    """Roll up availability, not channel gains or station validity."""
    result = []
    for station in sorted({r["station"] for r in rows}):
        members = [r for r in rows if r["station"] == station]
        record = {
            "station": station,
            "source_epochs": len(members),
            "saved_ppsd_epochs": sum(r["ppsd"] is not None for r in members),
            "approved_noise_epochs": sum(
                r["accepted_for_noise_model"] for r in members
            ),
        }
        for check in CHECKS:
            record[check + "_states"] = ";".join(
                sorted({r["checks"][check]["state"] for r in members})
            )
        record["rayleigh_number_available"] = any(
            r["checks"]["rayleigh"]["metrics"].get("pressure_multiplier") is not None
            for r in members
            if r["checks"]["rayleigh"]["binding"] != "none"
        )
        record["microseism_number_available"] = any(
            bool(r["checks"]["microseism"]["metrics"].get("microseism_conditional_g"))
            for r in members
        )
        record["local_p_number_available"] = any(
            bool(
                json.loads(
                    r["checks"]["local_p"]["metrics"].get("local_p_gain_by_band", "{}")
                    or "{}"
                )
            )
            for r in members
        )
        record["calibration_methods_with_numbers"] = sum(
            record[c + "_number_available"] for c in CHECKS[3:]
        )
        result.append(record)
    return result


def build_review(config, output, *, root=ROOT, decisions=None):
    inputs = Inputs(root)
    for code in (Path(__file__), ROOT / "../compy_pressure/inventory.py"):
        inputs.protect(code)
    config_doc = inputs.document(config)
    policy = inputs.document(config_doc["policy"])
    validate_policy(policy)
    networks = set(config_doc["networks"])
    records = load_masters(inputs, config_doc["masters"], networks)
    items = {
        key: base_review(key, row)
        for key, row in records.items()
        if key[0].split(".")[0] in networks and row.get("ppsd_npz_file")
    }
    unjoined = []
    for kind in ("microseism", "local_p"):
        attach_staging(
            inputs, items, records, config_doc[kind], kind, networks, unjoined
        )
    attach_rayleigh(inputs, items, records, config_doc["rayleigh"], networks, unjoined)
    for key, item in items.items():
        if resolve(key[3], root).is_file():
            item["inventory"] = source_record(inputs, key[3])
        else:
            item["issues"].append("missing_inventory")
        record = records.get(key, {})
        if not record:
            item["issues"].append("no_exact_master_metadata")
        else:
            add_master_evidence(inputs, item, record)
        item["evidence_fingerprint"] = json_hash(item)
    ordered = [items[k] for k in sorted(items)]
    if len({r["row_id"] for r in ordered}) != len(ordered):
        raise ValueError("Review identifier collision")
    if decisions:
        apply_decisions(ordered, inputs.document(decisions), policy)
    fields = list(KEYS) + [
        "row_id",
        "station",
        "evidence_fingerprint",
        "classification",
        "accepted_for_noise_model",
        "ppsd_file",
        "issues",
    ]
    for name in CHECKS:
        fields += [
            name + suffix
            for suffix in ("_state", "_source_status", "_verdict", "_score")
        ]
    flat = []
    for item in ordered:
        row = {k: item[k] for k in KEYS}
        row.update(
            {
                k: item[k]
                for k in (
                    "row_id",
                    "station",
                    "evidence_fingerprint",
                    "classification",
                    "accepted_for_noise_model",
                )
            }
        )
        row.update(
            ppsd_file=(item["ppsd"] or {}).get("path", ""),
            issues=";".join(item["issues"]),
        )
        for name, entry in item["checks"].items():
            row.update(
                {
                    name + "_" + k: entry[k]
                    for k in ("state", "source_status", "verdict", "score")
                }
            )
        flat.append(row)
    report = {
        "schema_version": 1,
        "policy": policy,
        "rows": ordered,
        "unjoined_evidence": unjoined,
        "protected_sha256": inputs.hashes,
        "production_master_updated": False,
        "gain_applied": False,
        "protected_inputs_unchanged": True,
        "summary": {
            "stations": len({r["station"] for r in ordered}),
            "full_source_epochs": len(ordered),
            "accepted_for_noise_model": sum(
                r["accepted_for_noise_model"] for r in ordered
            ),
            "with_saved_ppsd": sum(r["ppsd"] is not None for r in ordered),
            "unjoined_evidence_rows": len(unjoined),
            "check_states": {
                name: dict(Counter(r["checks"][name]["state"] for r in ordered))
                for name in CHECKS
            },
        },
    }
    inputs.verify()
    out = fresh_output(output, root)
    out.mkdir(parents=True)
    (out / "review.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    write_csv(out / "station_channel_review.csv", flat, fields)
    write_csv(
        out / "accepted_channels.csv",
        [r for r in flat if r["accepted_for_noise_model"]],
        fields,
    )
    stations = station_summary(ordered)
    if stations:
        write_csv(out / "station_summary.csv", stations, list(stations[0]))
    (out / "decisions_template.json").write_text(
        json.dumps(
            {
                "policy_version": policy["version"],
                "decisions": [],
                "pending_review": [
                    {
                        "row_id": item["row_id"],
                        "evidence_fingerprint": item["evidence_fingerprint"],
                        "reviewer": "",
                        "reviewed_utc": "",
                        "rationale": "",
                        "checks": {
                            name: {
                                "verdict": "unresolved",
                                "score": None,
                                "rationale": "",
                            }
                            for name in CHECKS
                        },
                        "noise_model_approved": False,
                        "nominal_response_accepted": False,
                        "independent_reference": {"citation": "", "description": ""},
                        "approved_frequency_band_hz": [],
                        "approved_time_range_utc": [],
                    }
                    for item in ordered
                ],
                "instructions": "Read docs/six_check_review.md before approval.",
            },
            indent=2,
        )
        + "\n"
    )
    inputs.verify()
    print(json.dumps(report["summary"], indent=2))
    print(f"Review (no production writes): {out}")
    return report


def validate_review(report):
    """Reject hand-edited acceptance flags or stale evidence in downstream tools."""
    if report.get("schema_version") != 1:
        raise ValueError("Unknown review schema")
    base = copy.deepcopy(report["rows"])
    decisions = []
    for item in base:
        decision = item.pop("review_decision", None)
        if decision:
            decisions.append(decision)
        fingerprint = item.pop("evidence_fingerprint")
        item["classification"] = "unresolved"
        item["accepted_for_noise_model"] = False
        for check in item["checks"].values():
            check.update(score=None, verdict="unresolved")
        if json_hash(item) != fingerprint:
            raise ValueError("Review row evidence fingerprint does not match")
        item["evidence_fingerprint"] = fingerprint
    apply_decisions(
        base,
        {"policy_version": report["policy"]["version"], "decisions": decisions},
        report["policy"],
    )
    if base != report["rows"]:
        raise ValueError(
            "Review decisions/classifications were changed without approval"
        )
    return base


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sources", type=Path, default=ROOT / "docs/quality_review_sources.json"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--decisions", type=Path)
    args = parser.parse_args()
    build_review(args.sources, args.output_dir, decisions=args.decisions)


if __name__ == "__main__":
    main()
