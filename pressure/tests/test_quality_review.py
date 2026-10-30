"""Six-check integration never equates calculation success with acceptance."""

import copy
import csv
import json

import pytest

from compy_pressure import quality_review as review


def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return path


def save_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def create_review_case(tmp_path):
    inventory = tmp_path / "inventory.xml"
    inventory.write_text("test-only source")
    ppsd = tmp_path / "ppsd.npz"
    ppsd.write_bytes(b"test-only saved spectra")
    row = dict(
        zip(
            review.KEYS,
            ["YV.RR52.00.BDH", "2020-01-01", "2021-01-01", str(inventory)],
            strict=True,
        )
    )
    row.update(
        input_units="PA",
        ppsd_units="Pa^2/Hz",
        ppsd_status="complete",
        ppsd_npz_file=str(ppsd),
        sample_rate_hz="4",
        Lat="1",
        Lon="2",
        Profondeur="1000",
    )
    for score in review.RESULT_COLUMNS:
        row[score] = ""
    for prefix in review.PREFIXES:
        result = save_json(
            tmp_path / (prefix + ".json"),
            {
                "seed_id": row["Seed_id"],
                "status": "comparison_only",
                "source_npz_sha256": review.digest(ppsd),
            },
        )
        row[prefix + "_result_file"] = str(result)
    save_csv(tmp_path / "master.csv", [row])
    policy = json.loads((review.ROOT / "docs/six_check_policy_v1.json").read_text())
    save_json(tmp_path / "policy.json", policy)
    micro = {k: row[k] for k in review.KEYS}
    micro.update(
        {
            review.RESULT_COLUMNS[4].removesuffix("_score")
            + "_status": "conditional_repeatable_unvalidated",
            "microseism_conditional_g": ".68",
        }
    )
    save_csv(tmp_path / "micro.csv", [micro])
    local = {k: row[k] for k in review.KEYS}
    local.update(
        master_epoch_match="exact",
        local_p_status="diagnostic_only",
        local_p_review_flags="legacy_noise_guard_unverified",
        local_p_gain_by_band='{"4-8":{"gain_median":0.55}}',
    )
    save_csv(tmp_path / "local.csv", [local])
    result = save_json(
        tmp_path / "ray.json",
        {
            "epoch_keys": [{k: row[k] for k in review.KEYS}],
            "review_flags": [{"independent_vertical": "not_verified"}],
        },
    )
    save_json(
        tmp_path / "rayleigh.json",
        {
            "networks": [
                {
                    "network": "YV",
                    "options": {"band": [0.03, 0.07]},
                    "stations": [
                        {
                            "station": "RR52",
                            "seed_id": row["Seed_id"],
                            "result_file": str(result),
                            "result_sha256": review.digest(result),
                            "pressure_multiplier": 0.67,
                            "status": "candidate_comparison_under_assumed_sign",
                        }
                    ],
                }
            ]
        },
    )
    config = {
        "networks": ["YV"],
        "masters": ["master.csv"],
        "policy": "policy.json",
        "microseism": "micro.csv",
        "local_p": "local.csv",
        "rayleigh": "rayleigh.json",
    }
    save_json(tmp_path / "sources.json", config)
    return tmp_path, row, policy


@pytest.fixture
def review_case(tmp_path):
    return create_review_case(tmp_path)


def run_review(case, name="review", decisions=None):
    root, _, _ = case
    return review.build_review(
        root / "sources.json",
        root / "results/REVIEWS" / name,
        root=root,
        decisions=decisions,
    )


def approval(item):
    checks = {
        name: {
            "verdict": "not_applicable",
            "score": None,
            "rationale": "synthetic test only",
        }
        for name in review.CHECKS
    }
    checks["expected_pressure"] = {
        "verdict": "pass",
        "score": 1,
        "rationale": "synthetic test only",
    }
    return {
        "row_id": item["row_id"],
        "evidence_fingerprint": item["evidence_fingerprint"],
        "reviewer": "TEST ONLY",
        "reviewed_utc": "2026-10-19T12:00:00Z",
        "rationale": "synthetic approval, not a real instrument",
        "checks": checks,
        "noise_model_approved": True,
        "nominal_response_accepted": True,
        "approved_frequency_band_hz": [0.05, 10],
        "approved_time_range_utc": ["2020-01-01T00:00:00Z", "2020-01-01T00:01:00Z"],
        "independent_reference": {
            "citation": "test-only reference",
            "description": "synthetic fixture",
        },
    }


def adopted_review(case):
    root, _, policy = case
    policy = {
        **policy,
        "status": "adopted",
        "approved_by": "TEST ONLY",
        "approved_utc": "2026-10-19T12:00:00Z",
        "rationale": "synthetic only",
    }
    save_json(root / "policy.json", policy)
    base = run_review(case, "adopted_base")
    decision = {
        "policy_version": policy["version"],
        "decisions": [approval(base["rows"][0])],
    }
    path = save_json(root / "decisions.json", decision)
    return run_review(case, "adopted", path)


def test_dry_run_preserves_inputs_and_all_six_scientific_limits(review_case):
    root, row, _ = review_case
    original = (root / "master.csv").read_bytes()
    report = run_review(review_case)
    assert (root / "master.csv").read_bytes() == original
    assert report["summary"]["stations"] == 1
    assert report["summary"]["accepted_for_noise_model"] == 0
    item = report["rows"][0]
    assert [item["checks"][c]["state"] for c in review.CHECKS] == [
        "diagnostic",
        "diagnostic",
        "diagnostic",
        "provisional",
        "conditional",
        "legacy_unvalidated",
    ]
    assert all(item["checks"][c]["score"] is None for c in review.CHECKS)
    assert review.validate_review(report) == report["rows"]
    assert (
        review.station_summary(report["rows"])[0]["calibration_methods_with_numbers"]
        == 3
    )
    assert not report["production_master_updated"]
    assert set(row[k] for k in review.KEYS)  # full identities survive.
    with pytest.raises(FileExistsError):
        run_review(review_case)


def test_explicit_synthetic_approval_can_flow_but_never_writes_master(review_case):
    root, _, _ = review_case
    before = (root / "master.csv").read_bytes()
    report = adopted_review(review_case)
    assert report["rows"][0]["accepted_for_noise_model"]
    assert review.validate_review(report) == report["rows"]
    assert (root / "master.csv").read_bytes() == before


@pytest.mark.parametrize(
    "status,expected",
    [
        ("no_passing_weeks", "unresolved"),
        ("no_suitable_events", "unresolved"),
        ("no_suitable_neighbors", "unresolved"),
        ("failed", "processing_error"),
        ("no_full_pressure_response", "unresolved"),
        ("not_applicable_sample_rate", "not_applicable"),
        ("new_unknown_status", "unresolved"),
    ],
)
def test_missing_or_failed_calculation_is_not_instrument_failure(status, expected):
    assert review.evidence_state(status) == expected


@pytest.mark.parametrize(
    "verdicts,expected",
    [
        (["not_tested", "not_applicable"], "not_assessed"),
        (["pass", "not_applicable"], "valid_under_reviewed_checks"),
        (["fail"], "unresolved"),
        (["fail", "fail"], "invalid_under_reviewed_checks"),
        (["pass", "fail"], "unresolved"),
        (["partial", "fail"], "unresolved"),
        (["pass", "unresolved"], "unresolved"),
    ],
)
def test_classification_is_categorical_not_an_average(verdicts, expected):
    assert review.classify(verdicts) == expected


def test_draft_policy_refuses_decisions(review_case):
    report = run_review(review_case)
    with pytest.raises(ValueError, match="draft"):
        review.apply_decisions(
            report["rows"],
            {
                "policy_version": report["policy"]["version"],
                "decisions": [approval(report["rows"][0])],
            },
            report["policy"],
        )


@pytest.mark.parametrize(
    "mutation", ["stale", "score", "reference", "nominal", "timezone", "unknown"]
)
def test_unsafe_approval_refused(review_case, mutation):
    report = adopted_review(review_case)
    item = report["rows"][0]
    decision = copy.deepcopy(item["review_decision"])
    if mutation == "stale":
        decision["evidence_fingerprint"] = "old"
    elif mutation == "score":
        decision["checks"]["expected_pressure"]["score"] = True
    elif mutation == "reference":
        decision["independent_reference"] = {}
    elif mutation == "nominal":
        decision["nominal_response_accepted"] = False
    elif mutation == "timezone":
        decision["reviewed_utc"] = "2026-10-19T12:00:00"
    else:
        decision["row_id"] = "unknown"
    with pytest.raises(ValueError):
        review.apply_decisions(
            [item],
            {"policy_version": report["policy"]["version"], "decisions": [decision]},
            report["policy"],
        )


def test_tampered_acceptance_or_evidence_is_refused(review_case):
    report = run_review(review_case)
    report["rows"][0]["accepted_for_noise_model"] = True
    with pytest.raises(ValueError, match="changed without"):
        review.validate_review(report)
    report["rows"][0]["checks"]["rayleigh"]["metrics"]["pressure_multiplier"] = 1
    with pytest.raises(ValueError, match="fingerprint"):
        review.validate_review(report)


def test_duplicate_source_epochs_and_stale_result_hashes_refused(review_case):
    root, row, _ = review_case
    save_csv(root / "master.csv", [row, row])
    with pytest.raises(ValueError, match="Duplicate master"):
        run_review(review_case)
    save_csv(root / "master.csv", [row])
    (root / "ray.json").write_text("{}")
    with pytest.raises(ValueError, match="hash changed"):
        run_review(review_case)


def test_join_never_falls_back_to_station_or_channel_only(review_case):
    root, _, _ = review_case
    local = list(csv.DictReader((root / "local.csv").open()))[0]
    local["master_epoch_match"] = "ambiguous"
    local["channel_start_time"] = ""
    save_csv(root / "local.csv", [local])
    report = run_review(review_case)
    assert len(report["unjoined_evidence"]) == 1
    assert report["rows"][0]["checks"]["local_p"]["state"] == "not_tested"


def test_output_and_identity_safety(tmp_path):
    with pytest.raises(ValueError, match="below results"):
        review.fresh_output(tmp_path / "exports", tmp_path)
    with pytest.raises(ValueError, match="identity"):
        review.full_key(
            dict(zip(review.KEYS, ["../A.STA..BDH", "x", "", "a"], strict=True))
        )
