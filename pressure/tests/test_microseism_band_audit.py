import hashlib
import json

import numpy as np
import pytest
from obspy import UTCDateTime
from obspy.core.inventory import Channel, Inventory, Network, Station
from obspy.core.inventory.response import Response

from compy_pressure.microseism_band_audit import (
    hypotheses,
    model_band_context,
    phase_band,
)
from scripts.audit_microseism_gain_bands import report, response_audit, table_export
from scripts.combine_microseism_weeks import InputGuard


def test_half_open_band_and_prescribed_hypothesis_is_not_reselected():
    f = np.linspace(0.03, 0.06, 31)
    phase = np.full(f.shape, 180.0)
    result = phase_band(f, phase, np.ones_like(f), [0.03, 0.05], "nominal")
    assert result["bins"] == 20
    assert result["status"] == "inconsistent_with_prescribed"
    assert result["prescribed_hypothesis"] == "nominal"
    assert result["prescribed_rms_deg"] == 180
    assert result["best_alternative_diagnostic"] == "reversed_nominal"


def test_phase_wrapping_and_coherence_selection_are_explicit():
    f = np.linspace(0.03, 0.05, 20, endpoint=False)
    phase = np.tile([-179.0, 179.0], 10)
    result = phase_band(f, phase, np.ones_like(f), [0.03, 0.05], "reversed_nominal")
    assert result["prescribed_rms_deg"] == pytest.approx(1)
    assert result["status"] == "consistent_with_prescribed"
    co = np.r_[np.ones(9), np.full(11, 0.49)]
    result = phase_band(f, phase, co, [0.03, 0.05], "reversed_nominal")
    assert result["coherent_bins"] == 9
    assert result["status"] == "insufficient_coherent_bins"


def test_phase_cutoffs_are_descriptive_and_empty_band_is_not_a_pass():
    f = np.linspace(0.03, 0.05, 20, endpoint=False)
    for degrees, expected in (
        (10, "consistent_with_prescribed"),
        (10.1, "approximate_prescribed"),
        (20.1, "inconsistent_with_prescribed"),
    ):
        r = phase_band(
            f, np.full_like(f, degrees), np.ones_like(f), [0.03, 0.05], "nominal"
        )
        assert r["status"] == expected
    result = phase_band(f, f, np.ones_like(f), [0.1, 0.2], "nominal")
    assert result["status"] == "insufficient_coherent_bins"
    assert result["prescribed_rms_deg"] is None


def test_reflected_phase_matches_complex_all_pass_ratio_with_vertical_sign():
    f = np.linspace(0.03, 0.05, 20, endpoint=False)
    omega = 2 * np.pi * f
    # Each nominal physical trace = true trace * R_true/R_nominal.
    # H_nominal/H_true = vertical bias / pressure bias because H = D/P.
    p, z = 0.1, 0.04
    pressure = (1j * omega - p) / (1j * omega + p)
    vertical = (1j * omega - z) / (1j * omega + z)
    expected = np.angle(vertical / pressure, deg=True)
    actual = hypotheses(f, [p], [z])["reflected_pole"]
    np.testing.assert_allclose(actual, expected, atol=1e-12)
    r = phase_band(
        f, expected, np.ones_like(f), [0.03, 0.05], "reflected_pole", [p], [z]
    )
    assert r["prescribed_rms_deg"] < 1e-10


@pytest.mark.parametrize("bad", [np.nan, -0.1, 1.1])
def test_bad_coherence_is_refused(bad):
    with pytest.raises(ValueError, match="invalid coherence"):
        phase_band(np.array([0.04]), [0], [bad], [0.03, 0.05], "nominal")


def test_missing_pole_invalid_grid_and_approximate_label_are_refused():
    for f, label, match in (
        ([0.04], "reflected_pole", "requires recorded"),
        ([0.04, 0.04], "nominal", "increase strictly"),
        ([0.04], "approximate_nominal", "exact frozen"),
    ):
        with pytest.raises(ValueError, match=match):
            phase_band(f, np.zeros(len(f)), np.ones(len(f)), [0.03, 0.05], label)


def test_real_coherent_contamination_can_pass_phase_while_biasing_amplitude():
    # Z = physical Z + 0.1*physical Z: perfectly coherent real contamination.
    # A phase test cannot establish that the amplitude is the desired signal.
    f = np.linspace(0.03, 0.05, 20, endpoint=False)
    h_true = np.ones_like(f, dtype=complex)
    contaminated = 1.1 * h_true
    result = phase_band(
        f, np.angle(contaminated, deg=True), np.ones_like(f), [0.03, 0.05], "nominal"
    )
    assert result["status"] == "consistent_with_prescribed"
    assert abs(contaminated[0] / h_true[0]) == pytest.approx(1.1)


def test_cached_model_context_is_same_bin_sensitivity_not_new_gain_fit():
    bins = {
        "g_mode0_water_h1": [0.7, 0.8],
        "g_mode0_water_h2": [0.8, 0.9],
        "g_mode0_no_water_h1": [0.707, 0.808],
        "g_mode0_no_water_h2": [0.808, 0.909],
        "upper_bound_h2": [0.9, 1.0],
    }
    r = model_band_context(bins, [True, False], reference=0.77)
    assert r["dry_vs_water_h1_pct"] == pytest.approx(1)
    assert r["reference_inside_model_h1_h2"]
    assert r["required_height_over_model_h1"] == pytest.approx(1.1)
    assert r["required_height_over_model_h2"] == pytest.approx(0.9625)
    with pytest.raises(ValueError, match="positive and finite"):
        model_band_context(bins, [False, False])


def test_response_reproduction_does_not_claim_absolute_calibration(tmp_path):
    start, end = UTCDateTime("2010-01-01"), UTCDateTime("2011-01-01")
    channels = []
    for code, units, dip in (("BDH", "PA", 0), ("BHZ", "M", -90)):
        channels.append(
            Channel(
                code,
                "00",
                1,
                2,
                -3000,
                0,
                dip=dip,
                sample_rate=1,
                start_date=start,
                end_date=end,
                response=Response.from_paz(
                    [], [], stage_gain=4, input_units=units, output_units="COUNTS"
                ),
            )
        )
    inv = Inventory(
        [Network("XX", stations=[Station("SYN", 1, 2, -3000, channels=channels)])],
        "test",
    )
    path = tmp_path / "response.xml"
    inv.write(str(path), format="STATIONXML")
    paired = {
        "response_snapshots": {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest()
        },
        "start_utc": "2010-06-01",
        "end_utc": "2010-06-08",
    }
    for name, code in (("pressure", "BDH"), ("vertical", "BHZ")):
        paired[f"{name}_epoch"] = {
            "Seed_id": f"XX.SYN.00.{code}",
            "channel_start_time": str(start),
            "channel_end_time": str(end),
        }
    data = {"frequency_hz": np.array([0.03, 0.04])}
    data["pressure_response"] = np.full(2, 4 + 0j)
    data["displacement_response"] = np.full(2, 4 + 0j)
    result = response_audit(paired, data, InputGuard())
    assert result["vertical"]["maximum_complex_relative_error"] < 1e-10
    assert result["vertical"]["absolute_response_independently_verified"] is False
    data["displacement_response"] *= 2
    with pytest.raises(ValueError, match="Saved complex response differs"):
        response_audit(paired, data, InputGuard())


def test_review_hash_drift_is_refused_before_output(tmp_path):
    protected = tmp_path / "saved_input"
    protected.write_text("changed")
    review = tmp_path / "review.json"
    review.write_text(
        json.dumps(
            {
                "frequency_rule": {"id": "claude-frequency-v1-20261015"},
                "protected_inputs_unchanged": True,
                "input_sha256": {str(protected): "wrong"},
            }
        )
    )
    output = tmp_path / "outputs"
    with pytest.raises(ValueError, match="Recorded input hash mismatch"):
        report(review, output)
    assert not output.exists()


def staging_fixture(tmp_path, *, checked=3, reportable=True):
    epoch = {
        "Seed_id": "XX.SYN.00.BDH",
        "channel_start_time": "2010-01-01",
        "channel_end_time": "2011-01-01",
        "inventory_file": "original.xml",
    }
    pair = tmp_path / "pair.json"
    pair.write_text(json.dumps({"pressure_epoch": epoch}))
    record = {
        "network": "XX",
        "station": "SYN",
        "seed_id": epoch["Seed_id"],
        "paired_result": str(pair),
        "status": "checked" if checked else "unresolved_pairing",
        "paired_reason": None if checked else "no common data",
    }
    (tmp_path / "summary.json").write_text(json.dumps({"stations": [record]}))
    row = {
        "station": "XX.SYN",
        "weeks_checked": checked,
        "weeks": [{"shift": 0, "run": str(tmp_path)}],
        "frequency_consistency": {
            "conditional_reportable": reportable,
            "median_g": 0.8,
            "passing_weeks": checked,
            "status": "stable" if reportable else "no_passing_weeks",
        },
    }
    return {"rows": [row]}, epoch


def test_table_staging_keeps_full_epoch_scores_blank_and_conditional_status(tmp_path):
    review, epoch = staging_fixture(tmp_path)
    row = table_export(review, [], InputGuard())[0]
    assert all(row[k] == v for k, v in epoch.items())
    assert row["microseism_conditional_g"] == 0.8
    assert row["calibration_check_pressure_vertical_microseisms_score"] is None
    assert (
        row["calibration_check_pressure_vertical_microseisms_status"]
        == "conditional_repeatable_unvalidated"
    )
    assert row["microseism_phase_audit_status"] == "not_audited"
    assert row["microseism_tilt_status"] == "not_assessed"
    assert json.loads((tmp_path / "pair.json").read_text()) == {"pressure_epoch": epoch}


def test_table_staging_missing_data_is_not_failure_or_a_reportable_gain(tmp_path):
    review, _ = staging_fixture(tmp_path, checked=0, reportable=False)
    row = table_export(review, [], InputGuard())[0]
    assert (
        row["calibration_check_pressure_vertical_microseisms_status"]
        == "unresolved_pairing"
    )
    assert row["microseism_conditional_g"] is None
    assert row["microseism_base_pair_reason"] == "no common data"
    assert row["calibration_check_pressure_vertical_microseisms_score"] is None


def test_table_staging_rejects_duplicate_epoch_and_misbound_station(tmp_path):
    review, _ = staging_fixture(tmp_path)
    review["rows"] *= 2
    with pytest.raises(ValueError, match="Duplicate full epoch"):
        table_export(review, [], InputGuard())
    review["rows"] = review["rows"][:1]
    pair = tmp_path / "pair.json"
    saved = json.loads(pair.read_text())
    saved["pressure_epoch"]["Seed_id"] = "XX.OTHER.00.BDH"
    pair.write_text(json.dumps(saved))
    with pytest.raises(ValueError, match="differs from source"):
        table_export(review, [], InputGuard())
