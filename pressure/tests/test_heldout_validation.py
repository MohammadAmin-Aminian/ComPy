"""Leakage prevention and frozen predictions on independent analytic cases."""

import numpy as np
import pytest
from obspy import UTCDateTime

from compy_pressure.heldout_validation import (
    frozen_gain_check,
    select_unused_events,
    summarize_errors,
)


def event(eid, hours, magnitude=6):
    return {
        "id": eid,
        "origin_utc": str(UTCDateTime("2020-01-01") + hours * 3600),
        "magnitude": magnitude,
        "depth_km": 20,
        "magnitude_type": "mww",
    }


def test_excludes_recorded_failures_and_origin_aliases_not_only_training_fits():
    reviewed = [dict(event("old", 0), status="rejected_qc")]
    selected, excluded = select_unused_events(
        [event("old", 100, 7), event("alias", 1, 7), event("new", 48)], reviewed
    )
    assert [e["id"] for e in selected] == ["new"]
    assert {e["reason"] for e in excluded} == {
        "previously_recorded_in_YV_trial",
        "within_24h_of_reviewed_event",
    }


def test_selects_all_mainshocks_without_three_event_cap_and_keeps_strongest_cluster():
    catalogue = [event(str(i), 48 * i) for i in range(8)] + [event("strong", 1, 7)]
    selected, _ = select_unused_events(catalogue, [])
    assert len(selected) == 8
    assert "strong" in {e["id"] for e in selected}
    assert "0" not in {e["id"] for e in selected}
    assert select_unused_events([event("boundary", 24)], [event("old", 0)])[0]


def spectra():
    return {
        "coherence_squared": np.full(9, 0.99),
        "cross_pa": np.ones(9, dtype=complex),
        "psd_ratio_pa": np.ones(9),
        "expected_pa": np.ones(9),
    }


SETTINGS = {"coherence_squared_threshold": 0.95, "phase_limit_deg": 20}


def test_wrong_frozen_gain_is_reported_without_refitting_or_amplitude_selection():
    data = spectra()
    before = data["psd_ratio_pa"].copy()
    result = frozen_gain_check(data, SETTINGS, 2.0, np.ones(9))
    assert result["status"] == "qc_usable"
    assert result["median_error_percent"] == pytest.approx(100)
    assert result["rms_error_percent"] == pytest.approx(100)
    assert result["frozen_gain"] == 2
    assert result["gain_refitted"] is False
    np.testing.assert_array_equal(before, data["psd_ratio_pa"])


def test_phase_and_coherence_failures_keep_prediction_diagnostics():
    data = spectra()
    data["cross_pa"] *= np.exp(1j * np.radians(45))
    result = frozen_gain_check(data, SETTINGS, 1.0, np.ones(9))
    assert result["status"] == "phase_review"
    assert result["all_band_median_error_percent"] == 0
    data["coherence_squared"][:] = 0.5
    assert (
        frozen_gain_check(data, SETTINGS, 1.0, np.ones(9))["status"]
        == "coherence_review"
    )


def test_primary_summary_excludes_arrival_contamination_and_reports_insufficient_data():
    base = {"reflected": {"status": "qc_usable", "median_error_percent": 10.0}}
    events = [
        base,
        dict(base, body_arrivals_in_window=["P"]),
        dict(base, other_packet_ids=["other"]),
    ]
    primary = summarize_errors(events)
    assert primary["usable_events"] == 1
    assert primary["status"] == "insufficient_heldout_events"
    assert summarize_errors(events, require_isolated=False)["usable_events"] == 2


def test_conditional_bootstrap_does_not_claim_absolute_calibration_confidence():
    events = [
        {"reflected": {"status": "qc_usable", "median_error_percent": 5.0}}
        for _ in range(6)
    ]
    result = summarize_errors(events)
    assert result["event_resampling_interval_of_median_error_percent"] == [5.0, 5.0]
    assert "not calibration confidence" in result["interval_note"]
