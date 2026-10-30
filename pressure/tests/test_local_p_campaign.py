import json

import pytest

from compy_pressure.local_p_check import PROTOCOL_VERSION
from scripts.export_local_p_results import band_json, flags
from scripts.run_local_p_campaign import archive_channels, summary_row, timing_options


def test_archive_channels_reads_names_from_both_layouts(tmp_path):
    (tmp_path / "BDH.D").mkdir()
    (tmp_path / "BDH.D" / "YV.RR52.00.BDH.D.2013.100").write_text("")
    (tmp_path / "waveforms").mkdir()
    (tmp_path / "waveforms" / "XO.WD48.--.LDH.20180819.mseed").write_text("")
    (tmp_path / "inventory.xml").write_text("")
    assert archive_channels(tmp_path) == {"BDH", "LDH"}


def test_summary_row_keeps_status_outcome_and_bands_separate():
    plan = {
        "network": "XO",
        "station": "LT02",
        "pressure_seed_id": "XO.LT02..HDH",
        "vertical_seed_id": "XO.LT02..HHZ",
        "plan_status": "eligible",
        "timing_mode": "drift_model",
    }
    summary = {
        "execution_status": "completed",
        "outcome": "conditional_gain_candidate",
        "event_reason_counts": {"accepted": 7, "lag_refinement_failed": 63},
        "regime_counts": {"steady_state": 7},
        "frequency_dependence": {"assessment": "consistent_within_spread"},
        "polarity_majority": 1,
        "bands": [
            {"band_hz": [1.0, 2.0], "events": 7, "gain_median": 1.02, "gain_mad": 0.04},
            {"band_hz": [8.0, 16.0], "events": 0},
        ],
    }
    row = summary_row(plan, summary, "dir", "log", "ran")
    assert row["execution_status"] == "completed"
    assert row["outcome"] == "conditional_gain_candidate"
    assert row["accepted_events"] == 7
    assert row["1_2hz_gain_median"] == 1.02 and row["8_16hz_events"] == 0
    assert "score" not in json.dumps(row)
    assert summary_row({**plan, "plan_status": "not_applicable_sample_rate"}) == {
        **{k: plan.get(k, "") for k in row if k in plan or k == "note"},
        "plan_status": "not_applicable_sample_rate",
        "pressure_rate_hz": "",
        "vertical_rate_hz": "",
        "inventory": "",
    }


def test_export_flags_and_bands():
    summary = {
        "protocol_version": PROTOCOL_VERSION,
        "regime_counts": {"pre_multiple": 3},
        "pressure_epochs_accepted": [["a", "b"]],
        "bands": [
            {
                "band_hz": [2.0, 4.0],
                "events": 3,
                "gain_median": 0.78,
                "regimes": ["pre_multiple"],
            },
            {"band_hz": [4.0, 8.0], "events": 0},
        ],
    }
    expected = [
        "pre_multiple_diagnostic_only",
        "pre_multiple_filter_leakage_not_field_validated",
    ]
    assert flags(summary, {"lag_search_s": 0}) == expected
    assert "timing_lag_estimated_phase_not_independent" in flags(
        summary, {"lag_search_s": 150.0}
    )
    drift = {**summary, "timing_mode": "drift_model"}
    assert "timing_lag_estimated_phase_not_independent" in flags(drift, {})
    assert (
        flags({**summary, "timing_mode": "none"}, {"lag_search_s": 150.0}) == expected
    )
    assert list(band_json(summary)) == ["2-4"]


def test_timing_options_follow_recorded_decisions():
    decisions = {
        "XO.LT02..HDH": {"mode": "drift_model", "drift": {"slope_s_per_day": -0.48}},
        "XO.LA33..HDH": {"mode": "none"},
    }
    options, mode, note = timing_options("XO.LT02..HDH", "APG_logger", decisions)
    assert mode == "drift_model" and options["timing_drift"]["slope_s_per_day"] == -0.48
    options, mode, note = timing_options("XO.LA33..HDH", "APG_logger", decisions)
    assert mode == "none" and options["timing_lag_search_s"] == 0 and note == ""
    _, mode, note = timing_options("XO.LA99..HDH", "Parosci APG/LDEO_APG_logger", {})
    assert mode == "none" and note == "separate_logger_timing_unverified"
    with pytest.raises(ValueError):
        timing_options("X", "", {"X": {"mode": "guess"}})


def test_drift_fit_recovers_line_and_rejects_outlier():
    import numpy as np
    from obspy import UTCDateTime

    from scripts.fit_local_p_clock_drift import fit

    reference = UTCDateTime("2018-06-01")
    days = np.linspace(0, 200, 12)
    rng = np.random.default_rng(7)
    lags = -4.6 - 0.4834 * days + 0.05 * rng.standard_normal(days.size)
    lags[5] += 1.5  # one bad refinement
    times = [reference + d * 86400 for d in days]
    slope, intercept, _, inlier, residual = fit(times, lags, reference, 1.0)
    assert slope == pytest.approx(-0.4834, abs=0.002)
    assert intercept == pytest.approx(-4.6, abs=0.1)
    assert not inlier[5] and inlier.sum() == 11
    assert np.sqrt(np.mean(residual**2)) < 0.1


def test_summary_tables_read_saved_microseism_and_rayleigh_values():
    from scripts.summarize_local_p_campaign import (
        microseism_values,
        rayleigh_values,
    )

    micro = microseism_values()
    assert micro["YV.RR40"][0] == pytest.approx(0.748529445)
    assert micro["YV.RR40"][1] == "frozen rule (conditional_unvalidated)"
    assert "YV.RR36" not in micro
    assert len(micro) == 7
    rayleigh = rayleigh_values()
    assert rayleigh["XO.WD66"][1] == "XO.WD66..LDH"
    assert rayleigh["YV.RR52"][0] == pytest.approx(0.6787, abs=1e-4)
