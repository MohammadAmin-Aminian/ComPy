"""The explicit local-P updater must never alter other checks or scores."""

import pandas as pd
import pytest

from compy_pressure.inventory import (
    load_master_dataset,
    save_master_dataset,
    update_local_p_results,
)


def master(tmp_path, duplicate=False):
    row = {
        "Seed_id": "XX.A..HDH",
        "channel_start_time": "2020-01-01",
        "channel_end_time": "2020-02-01",
        "inventory_file": "one.xml",
        "sample_rate_hz": "100",
        "calibration_check_pressure_vertical_local_p_waves_score": "",
        "calibration_check_pressure_vertical_local_p_waves_status": "not_tested",
        "rayleigh_pressure_multiplier": ".7",
        "notes": "preserve",
    }
    rows = [
        row,
        dict(row) if duplicate else {**row, "channel_start_time": "2021-01-01"},
    ]
    save_master_dataset(pd.DataFrame(rows), tmp_path)
    return {
        k: row[k]
        for k in ["Seed_id", "channel_start_time", "channel_end_time", "inventory_file"]
    }


def test_local_p_updater_preserves_scores_statuses_other_epochs_and_fields(tmp_path):
    key = master(tmp_path)
    before = load_master_dataset(tmp_path)
    update_local_p_results(
        tmp_path,
        [
            {
                **key,
                "local_p_status": "diagnostic_only",
                "local_p_gain_by_band": '{"1-2":{"gain_median":0.8}}',
            }
        ],
    )
    after = load_master_dataset(tmp_path)
    assert after[before.columns].equals(before)
    assert after.local_p_status.tolist() == ["diagnostic_only", ""]
    assert list((tmp_path / "backups").glob("*/master_inventory.csv"))


@pytest.mark.parametrize(
    "field",
    [
        "rayleigh_status",
        "calibration_check_pressure_vertical_local_p_waves_score",
        "calibration_check_pressure_vertical_local_p_waves_status",
        "local_p_fake_score",
    ],
)
def test_local_p_updater_rejects_other_results_and_scores(tmp_path, field):
    key = master(tmp_path)
    before = (tmp_path / "master_inventory.csv").read_bytes()
    with pytest.raises(ValueError, match="only local_p"):
        update_local_p_results(tmp_path, [{**key, field: 1}])
    assert (tmp_path / "master_inventory.csv").read_bytes() == before


def test_local_p_updater_refuses_duplicate_or_unmatched_keys(tmp_path):
    key = master(tmp_path)
    update = {**key, "local_p_status": "unresolved"}
    with pytest.raises(ValueError, match="duplicate"):
        update_local_p_results(tmp_path, [update, update])
    with pytest.raises(ValueError, match="exactly one"):
        update_local_p_results(tmp_path, [{**update, "inventory_file": "other.xml"}])
    with pytest.raises(ValueError, match="four"):
        update_local_p_results(
            tmp_path, [{"Seed_id": "XX.A..HDH", "local_p_status": "unresolved"}]
        )


def test_local_p_updater_refuses_duplicate_master_rows(tmp_path):
    key = master(tmp_path, duplicate=True)
    with pytest.raises(ValueError, match="exactly one"):
        update_local_p_results(tmp_path, [{**key, "local_p_status": "unresolved"}])
