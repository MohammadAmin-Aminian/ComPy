"""All-network trial routing and preservation of baseline exclusions."""

import json

import pytest

from scripts.run_calibration_experiment import (
    DEFAULTS,
    Guard,
    comparison_rows,
    selected_stations,
    station_trial,
)


def test_all_station_selection_includes_workflow_stations_outside_small_defaults():
    workflow = {
        "entries": {
            str(i): {"seed_id": f"YV.{station}.00.BDH"}
            for i, station in enumerate(["RR03", "RR41", "RR57"])
        }
    }
    assert selected_stations(workflow, "YV") == DEFAULTS["YV"]
    assert selected_stations(workflow, "YV", all_stations=True) == [
        "RR03",
        "RR41",
        "RR57",
    ]
    assert selected_stations(workflow, "YV", ["RR57"]) == ["RR57"]
    assert selected_stations(workflow, "7D") == ["RR03", "RR41", "RR57"]


def test_reference_band_exclusion_retained_without_loading_master(tmp_path):
    source = tmp_path / "original.json"
    source.write_text(
        json.dumps(
            {
                "status": "reference_band_unusable",
                "settings": {},
                "reason": "Model pole preflight",
                "events": [],
            }
        )
    )
    before = source.read_bytes()
    guard = Guard()
    report = station_trial(
        "YV",
        "RR57",
        {"update": {"rayleigh_result_file": str(source)}},
        {},
        tmp_path,
        tmp_path,
        tmp_path / "new_trial",
        guard,
    )
    assert report["baseline_reason"] == "Model pole preflight"
    assert all(
        v["status"] == "reference_band_unusable"
        and v["median_pressure_multiplier"] is None
        for v in report["summaries"].values()
    )
    guard.verify()
    assert source.read_bytes() == before


@pytest.mark.parametrize("network,coefficient", [("7D", 2.0), ("XO", 2.0), ("YL", 0.5)])
def test_network_comparison_direction_and_no_unverified_7d_xo_benchmark(
    network, coefficient
):
    rows = comparison_rows(
        [
            {
                "network": network,
                "station": "UNPUBLISHED",
                "saved_station_gain": 2.0,
                "summaries": {},
            }
        ]
    )
    assert rows[0]["comparison_coefficient"] == coefficient
    assert rows[0]["published_coefficient"] is None
