"""Tests for completed-campaign microseism comparisons."""

import json

import pytest

from compy_pressure import microseism_compare as comparison


def row(station, coherence, phase, nominal, fixed, selected, available, ruan):
    return {
        "network": station.split(".")[0],
        "station": station,
        "retained_segments": 144,
        "selected_days": selected,
        "available_days": available,
        "median_coherence_squared": coherence,
        "ruan_high_coherence_bins": int(ruan > 0),
        "ruan_high_coherence_fraction": ruan,
        "median_nominal_transfer_amplitude_m_pa": nominal,
        "median_fixed_reference_transfer_amplitude_m_pa": fixed,
        "median_nominal_transfer_phase_deg": phase,
        "phase_family": comparison.phase_family(phase)[0],
        "published_coefficient": None,
        "first_method_publication_difference_percent": None,
        "publication_comparison_status": "none",
    }


def test_trailing_json_parser_requires_one_final_object(tmp_path):
    path = tmp_path / "station.log"
    path.write_text('warning before result\n{"station": "XX.ONE"}\n')
    assert comparison.parse_trailing_json(path) == {"station": "XX.ONE"}
    path.write_text("warning only\n")
    with pytest.raises(ValueError, match="trailing JSON"):
        comparison.parse_trailing_json(path)


def test_phase_families_and_circular_mean():
    assert comparison.phase_family(-7) == ("near_0_deg", 7)
    family, distance = comparison.phase_family(-175)
    assert family == "near_180_deg"
    assert distance == 5
    assert abs(abs(comparison.circular_mean_degrees([179, -179])) - 180) < 1e-10


def test_group_summary_keeps_references_descriptive():
    rows = [
        row("XX.ONE", 0.9, 4, 1, 2, 2, 2, 0.2),
        row("XX.TWO", 0.7, 176, 4, 4, 0, 2, 0),
    ]
    result = comparison.summarize_group(rows)
    assert result["stations_with_weekly_median_coherence_at_least_0_8"] == 1
    assert result["stations_with_any_ruan_gt_0_95_bins"] == 1
    assert result["saikia_daily_reference"]["selected_days"] == 2
    assert result["saikia_daily_reference"]["stations_with_none"] == 1
    assert result["phase"]["family_counts"] == {
        "near_0_deg": 1,
        "near_180_deg": 1,
    }
    assert (
        result["amplitude_dispersion"]["fixed_reference_geometric_spread_factor"]
        < result["amplitude_dispersion"]["nominal_geometric_spread_factor"]
    )


def test_publication_statistics_preserve_direction_and_c07_exception():
    rows = [
        row("YV.RR52", 0.8, 0, 1, 1, 1, 1, 0),
        row("YL.A02W", 0.8, 180, 1, 1, 1, 1, 0),
        row("YL.C07W", 0.8, 180, 1, 1, 1, 1, 0),
    ]
    values = [(0.67, 1.3), (1.5, -4.2), (1.99, 86.8)]
    for item, (published, difference) in zip(rows, values, strict=True):
        item["published_coefficient"] = published
        item["first_method_publication_difference_percent"] = difference
        item["publication_comparison_status"] = "matched"
    stats = comparison.publication_statistics(rows)
    assert stats["YV_Aminian2025_same_g_direction"]["stations"] == 1
    peers = stats["YL_ZhaWebb2016_reciprocal_direction_excluding_C07W"]
    assert peers["stations"] == 1
    assert peers["median_difference_percent"] == -4.2


def test_json_fixture_is_serializable():
    payload = comparison.summarize_group([row("XX.ONE", 0.9, 4, 1, 2, 2, 2, 0.2)])
    json.dumps(payload)
