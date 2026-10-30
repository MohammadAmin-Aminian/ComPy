"""Regression contracts from the independent local-P review (no field retuning)."""

import json
from pathlib import Path

import pytest

from compy_pressure import local_p_check as lc
from compy_pressure.local_p_io import fresh_file, local_p_output
from scripts import export_local_p_results as export
from scripts import run_local_p_campaign as campaign


@pytest.mark.parametrize("earliest,expected", [(1000.7, 800), (1500.1, 1200)])
def test_noise_window_precedes_uncertain_p_and_signal_gap(earliest, expected):
    # Signal starts at 1500, gap 100 samples, window length 200.
    start = lc.noise_window_start(1500, 200, 100, earliest, 1.0)
    assert start == expected
    assert start + 200 <= earliest
    assert start + 200 <= 1400


def test_uncertainty_moves_noise_not_signal():
    start = lc.noise_window_start(1500, 200, 100, 1000.7, 1.0)
    assert start + 200 == 1000
    assert lc.noise_window_start(1500, 200, 100, 2000, 1.0) == 1200


def test_epoch_summary_never_copies_gain_to_unaccepted_epoch():
    s = {
        "outcome": "conditional_gain_candidate",
        "pressure_epochs_in_coverage": [["a", "b"], ["c", "d"]],
        "pressure_epochs_accepted": [["a", "b"]],
        "bands": [{"events": 3}],
        "frequency_dependence": {"scalar_median_all_bands": 1.2},
    }
    assert export.epoch_summary(s, ["a", "b"])["bands"]
    other = export.epoch_summary(s, ["c", "d"])
    assert other["bands"] == []
    assert other["frequency_dependence"] == {}
    assert other["outcome"] == "epoch_attribution_unresolved"
    assert (
        export.epoch_summary(
            {**s, "pressure_epochs_accepted": [["a", "b"], ["c", "d"]]}, ["a", "b"]
        )["bands"]
        == []
    )


def test_inventory_join_requires_source_not_copied_master_path(tmp_path, monkeypatch):
    source, unrelated = tmp_path / "source.xml", tmp_path / "other.xml"
    source.write_text("source")
    unrelated.write_text("different")
    monkeypatch.setattr(
        export, "channel_identity", lambda *a: ((1.0, 2.0, 3.0, 4.0), "PA")
    )
    summary = {
        "inputs": {str(source): {"sha256": export.sha256(source)}},
        "pressure_epochs_in_coverage": [["a", "b"]],
    }
    row = {"pressure_seed_id": "XX.A..HDH"}
    key = {
        "Seed_id": row["pressure_seed_id"],
        "channel_start_time": "a",
        "channel_end_time": "b",
    }
    candidates = [
        {**key, "inventory_file": str(source)},
        {**key, "inventory_file": str(unrelated)},
    ]
    matched, _ = export.match_master_epochs(summary, row, candidates)
    assert matched == candidates[:1]
    with pytest.raises(ValueError, match="Duplicate"):
        export.match_master_epochs(summary, row, candidates[:1] * 2)
    mapping = {
        "source_inventory": str(source),
        "master_inventory": str(unrelated),
        "source_sha256": export.sha256(source),
        "master_sha256": "bad",
        "seed_ids": [row["pressure_seed_id"]],
    }
    with pytest.raises(ValueError, match="hash mismatch"):
        export.match_master_epochs(summary, row, candidates[1:], [mapping])


def test_polarity_tie_is_not_positive():
    band = {
        "band_hz": [1.0, 2.0],
        "gain_noise_corrected": 1.0,
        "gain_h1_upper": 1.0,
        "gain_h2_lower": 1.0,
        "gain_model_low": 0.9,
        "gain_model_high": 1.1,
        "phase_deg": 0.0,
    }
    events = [
        {
            "event_id": str(i),
            "status": "accepted",
            "regime": "steady_state",
            "bands": [{**band, "polarity": p}],
        }
        for i, p in enumerate([1, -1, 0])
    ]
    result = lc.station_summary(events, lc.DEFAULTS, 3)
    assert result["polarity_majority"] is None
    assert result["polarity_classified_fraction"] == pytest.approx(2 / 3)
    assert "not confidence" in result["polarity_summary_basis"]


def test_fingerprint_tracks_root_added_and_modified_waveforms(tmp_path, monkeypatch):
    root = tmp_path / "waves"
    root.mkdir()
    wave = root / "XX.A..HDH.D.2020.001"
    wave.write_bytes(b"one")
    other = root / "XX.A..HHZ.D.2020.001"
    other.write_bytes(b"one")
    inv, cat = tmp_path / "inv.xml", tmp_path / "cat.csv"
    inv.write_text("inventory")
    cat.write_text("catalogue")
    monkeypatch.setattr(campaign, "CODE_FILES", [])
    row = {
        "pressure_seed_id": "XX.A..HDH",
        "vertical_seed_id": "XX.A..HHZ",
        "inventory": str(inv),
        "catalog": str(cat),
        "waveform_dir": str(root),
    }
    initial = campaign.fingerprint(row, {})
    wave.write_bytes(b"changed length")
    changed = campaign.fingerprint(row, {})
    assert changed != initial
    (root / "XX.A..HDH.D.2020.002").write_bytes(b"added")
    assert campaign.fingerprint(row, {}) != changed
    relocated = tmp_path / "relocated"
    root.rename(relocated)
    assert campaign.fingerprint({**row, "waveform_dir": str(relocated)}, {}) != changed


def test_fresh_files_protect_sidecars_and_raw_roots(tmp_path):
    good = tmp_path / "results/LOCAL_P/catalogs/catalog.csv"
    assert fresh_file(good) == good
    good.parent.mkdir(parents=True)
    good.write_text("existing")
    with pytest.raises(FileExistsError):
        fresh_file(good)
    with pytest.raises(ValueError):
        local_p_output(Path.home() / "Desktop/Data_Bruit_FM/results/LOCAL_P/unsafe")
    with pytest.raises(ValueError):
        local_p_output(tmp_path / "results/LOCAL_P/../../wrong")


def test_old_protocol_cannot_export_as_current_candidate():
    assert (
        export.result_status({"outcome": "conditional_gain_candidate"}, {})
        == "legacy_protocol_unvalidated"
    )
    assert (
        export.result_status(
            {
                "outcome": "conditional_gain_candidate",
                "protocol_version": lc.PROTOCOL_VERSION,
            },
            {},
        )
        == "conditional_gain_candidate"
    )


def test_parent_failure_does_not_stop_next_pair(tmp_path, monkeypatch):
    out = tmp_path / "results/LOCAL_P/campaign"
    rows = [
        {
            "network": "YV",
            "station": s,
            "pressure_seed_id": f"YV.{s}..BDH",
            "vertical_seed_id": f"YV.{s}..BHZ",
            "plan_status": "eligible",
        }
        for s in ["A", "B"]
    ]
    monkeypatch.setattr(campaign, "plan_network", lambda *a: rows)
    calls = []

    def run(row, *a):
        calls.append(row["station"])
        if row["station"] == "A":
            raise OSError("parent failure")
        return (
            {"execution_status": "completed", "outcome": "unresolved"},
            None,
            None,
            "ran",
        )

    monkeypatch.setattr(campaign, "run_pair", run)
    monkeypatch.setattr(
        "sys.argv", ["campaign", "--campaign-dir", str(out), "--networks", "YV"]
    )
    campaign.main()
    report = json.loads((out / "campaign_summary.json").read_text())
    assert calls == ["A", "B"]
    assert report["execution_counts"] == {"failed": 1, "completed": 1}


def test_epoch_specific_results_override_pooled_station_numbers():
    s = {
        "pressure_epochs_in_coverage": [["a", "b"], ["c", "d"]],
        "pressure_epochs_accepted": [["a", "b"], ["c", "d"]],
        "bands": [{"gain_median": 99}],
        "epoch_results": [
            {
                "pressure_epoch": ["a", "b"],
                "outcome": "diagnostic_only",
                "bands": [{"gain_median": 1.2, "event_ids": ["e1", "e2"]}],
            }
        ],
    }
    scoped = export.epoch_summary(s, ["a", "b"])
    assert scoped["bands"][0]["gain_median"] == 1.2
    assert scoped["event_reason_counts"]["accepted"] == 2
    assert export.epoch_summary(s, ["c", "d"])["bands"] == []


def test_catalogue_sidecar_collision_refused_before_network(tmp_path, monkeypatch):
    from scripts import fetch_local_catalog as fetch

    out = tmp_path / "results/LOCAL_P/catalogs/cat.csv"
    out.parent.mkdir(parents=True)
    out.with_suffix(".json").write_text("preserve")
    monkeypatch.setattr(
        fetch, "urlopen", lambda *a, **kw: pytest.fail("network called")
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "fetch",
            "--inventory",
            "unused.xml",
            "--start",
            "2020-01-01",
            "--end",
            "2021-01-01",
            "--min-magnitude",
            "3",
            "--output",
            str(out),
        ],
    )
    with pytest.raises(FileExistsError):
        fetch.main()
    assert out.with_suffix(".json").read_text() == "preserve"
    assert not out.exists()


def test_summary_existing_output_refused_before_input_read(tmp_path, monkeypatch):
    from scripts import summarize_local_p_campaign as summarize

    out = tmp_path / "results/LOCAL_P/table.csv"
    out.parent.mkdir(parents=True)
    out.write_text("existing")
    monkeypatch.setattr(
        "sys.argv",
        [
            "summarize",
            "--campaign-dir",
            "missing",
            "--csv",
            str(out),
            "--markdown",
            str(out.with_suffix(".md")),
        ],
    )
    with pytest.raises(FileExistsError):
        summarize.main()
    assert out.read_text() == "existing"


def test_drift_fit_refuses_degenerate_inliers():
    from obspy import UTCDateTime

    from scripts.fit_local_p_clock_drift import fit

    t = UTCDateTime("2020-01-01")
    with pytest.raises(ValueError, match="distinct"):
        fit([t, t + 86400, t + 172800], [0.0, 5.0, 0.0], t, -1.0)
