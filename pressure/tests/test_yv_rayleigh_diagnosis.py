"""Saved-input-only YV Rayleigh diagnosis: response hypotheses and safe outputs."""

import json
from pathlib import Path

import numpy as np
import pytest
from obspy import UTCDateTime
from obspy.core.inventory import Channel, Inventory, Network, Site, Station
from obspy.core.inventory.response import Response

from compy_pressure.calibration_experiment import pole_factor
from scripts.diagnose_yv_rayleigh import (
    alternative_stages,
    coherent_event,
    conjugate_repair,
    event_category,
    hypothesis_terms,
    reflect_positive,
    run,
    snapshot_channel,
    stage_ratio,
    unpaired_roots,
    wrap,
)

P = -0.074016 + 0.07347j
Q = -502.65 + 596.9j
GOOD_Z = [([P, P.conjugate(), Q, Q.conjugate()], [0j, 0j])]
BAD_Z = [([P, P, Q, Q], [0j, 0j])]


def test_unpaired_roots_flags_same_sign_duplicates_only():
    assert unpaired_roots([P, P.conjugate(), -3.0, 0j]) == []
    assert unpaired_roots([P, P]) == [P, P]
    assert len(unpaired_roots([P, P, Q, Q.conjugate()])) == 2


def test_conjugate_repair_restores_pairs_and_refuses_ambiguity():
    repaired = conjugate_repair([P, P, Q, Q])
    assert unpaired_roots(repaired) == []
    assert sorted(repaired, key=lambda r: (r.real, r.imag)) == sorted(
        [P, P.conjugate(), Q, Q.conjugate()], key=lambda r: (r.real, r.imag)
    )
    assert conjugate_repair(GOOD_Z[0][0]) == GOOD_Z[0][0]
    with pytest.raises(ValueError):
        conjugate_repair([P, -1.0])


def test_positive_complex_pressure_pole_is_not_silently_made_real():
    with pytest.raises(ValueError, match="real pressure pole"):
        reflect_positive([0.1 + 0.2j])


def test_reflection_ratio_matches_completed_trial_factor():
    f = np.linspace(0.03, 0.07, 9)
    pressure = [([0.0943 + 0j], [0j])]
    reflected = alternative_stages(pressure, GOOD_Z)["H1"][0]
    assert np.allclose(stage_ratio(f, pressure, reflected), pole_factor(f, [0.0943]))


def synthetic_phase(f, pressure_true, vertical_true, pressure_meta, vertical_meta):
    """Phase saved by nominal processing when the physical P/A relation is real."""
    r_p = stage_ratio(f, pressure_meta, pressure_true)
    r_z = stage_ratio(f, vertical_meta, vertical_true)
    return wrap(-np.angle(r_p, deg=True) + np.angle(r_z, deg=True)), abs(r_z / r_p)


def test_vertical_repair_explains_phase_and_amplitude_of_defective_metadata():
    f = np.linspace(0.03, 0.07, 9)
    meta_p = [([0.1026 + 0j], [0j])]
    true_p = [([-0.1026 + 0j], [0j])]
    observed, amplitude = synthetic_phase(f, true_p, GOOD_Z, meta_p, BAD_Z)
    terms = hypothesis_terms(f, meta_p, BAD_Z)
    assert np.allclose(wrap(observed - terms["H2"][0]), 0, atol=1e-9)
    assert np.allclose(terms["H2"][1], amplitude)
    assert np.min(np.abs(wrap(observed - terms["H1"][0]))) > 90
    # Reversed polarity is exactly 180 degrees away from H2.
    assert np.allclose(abs(wrap(observed - terms["H3"][0])), 180)


def test_pole_divided_by_ten_is_a_fixed_alternative_not_a_fit():
    f = np.linspace(0.03, 0.07, 9)
    meta_p = [([0.747 + 0j], [0j])]
    true_p = [([-0.0747 + 0j], [0j])]
    observed, _ = synthetic_phase(f, true_p, GOOD_Z, meta_p, GOOD_Z)
    terms = hypothesis_terms(f, meta_p, GOOD_Z)
    assert np.allclose(wrap(observed - terms["H4"][0]), 0, atol=1e-9)
    assert np.all(np.abs(wrap(observed - terms["H1"][0])) > 40)
    # Hypotheses without a positive pole or a defect leave nominal unchanged.
    insu = hypothesis_terms(f, [([-0.0182 + 0j], [0j])], GOOD_Z)
    for name in ("H1", "H2", "H4"):
        assert np.allclose(insu[name][0], insu["H0"][0])


def test_coherent_event_rule_and_outcome_categories_follow_production_labels():
    assert coherent_event([0.96] * 9, 0.95)
    assert not coherent_event([0.99] * 5 + [0.5] * 4, 0.95)  # mean below
    assert event_category({"status": "excluded_cluster"}) == "selection_excluded_24h"
    assert event_category({"status": "rejected_qc"}) == "coherence_failed"
    assert event_category({"status": "polarity_review"}) == "phase_review"
    assert (
        event_category(
            {"status": "unresolved", "reason": "Grid minimum lies on a boundary"}
        )
        == "fit_bound_unresolved"
    )
    assert (
        event_category({"status": "unresolved", "reason": "Missing, gapped window"})
        == "input_refused"
    )


def write_snapshot(path, code, poles, zeros, units):
    start = UTCDateTime("2012-10-01")
    response = Response.from_paz(
        zeros=zeros, poles=poles, stage_gain=1.0, input_units=units, output_units="V"
    )
    channel = Channel(
        code=code,
        location_code="00",
        latitude=0,
        longitude=0,
        elevation=-4000,
        depth=0,
        sample_rate=50,
        start_date=start,
        end_date=UTCDateTime("2012-10-03"),
        response=response,
        dip=-90 if code == "BHZ" else 0,
        azimuth=0,
    )
    station = Station("RR99", 0, 0, -4000, channels=[channel], site=Site("RR99"))
    Inventory([Network("YV", stations=[station])]).write(str(path), "STATIONXML")


def synthetic_saved_products(tmp_path):
    f = np.linspace(0.03, 0.07, 9)
    meta_p = ([0.1026 + 0j], [0j])
    phase, amplitude = synthetic_phase(
        f, [([-0.1026 + 0j], [0j])], GOOD_Z, [meta_p], BAD_Z
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    snapshots = {}
    for code, (poles, zeros), units in (
        ("BDH", meta_p, "PA"),
        ("BHZ", BAD_Z[0], "M/S"),
    ):
        folder = tmp_path / f"YV.RR99.00.{code}__20121001T000000__x"
        folder.mkdir()
        write_snapshot(folder / "inventory.xml", code, poles, zeros, units)
        snapshots[str(folder / "inventory.xml")] = None
    from scripts.diagnose_yv_rayleigh import sha256

    snapshots = {k: sha256(k) for k in snapshots}
    events, trial_events = [], []
    for i in range(3):
        spectra = run_dir / f"e{i}_spectra.npz"
        np.savez(
            spectra,
            frequency_hz=f,
            phase_deg=phase,
            coherence_squared=np.full(9, 0.99),
            psd_ratio_pa=amplitude * 2.0,
            expected_pa=np.full(9, 2.0),
        )
        events.append(
            {
                "id": f"e{i}",
                "status": "polarity_review",
                "files": {"spectra": str(spectra)},
                "tilt": {"vertical_component": "YV.RR99.00.BHZ"},
            }
        )
        trial_events.append(
            {
                "id": f"e{i}",
                "reflected_phase_deg": list(
                    wrap(phase + np.angle(pole_factor(f, [0.1026]), deg=True))
                ),
            }
        )
    events.append({"id": "x", "status": "excluded_cluster"})
    epoch = {
        "Seed_id": "YV.RR99.00.BDH",
        "channel_start_time": "2012-10-01T00:00:00.000000Z",
        "channel_end_time": "2012-10-03T00:00:00.000000Z",
        "inventory_file": str(tmp_path / "source.xml"),
    }
    vertical_epoch = dict(epoch, Seed_id="YV.RR99.00.BHZ")
    comparison = run_dir / "comparison.json"
    comparison.write_text(
        json.dumps(
            {
                "station": "YV.RR99",
                "epoch_keys": [epoch, vertical_epoch],
                "reason": "Fewer than 3 QC-accepted independent events",
                "candidate_events": 0,
                "settings": {
                    "coherence_squared_threshold": 0.95,
                    "depth_m": 4000,
                    "phase_limit_deg": 20.0,
                    "band_hz": [0.03, 0.07],
                    "pressure_sign": -1,
                    "tilt_method": "none",
                },
                "response_snapshot_sha256": snapshots,
                "magnitude_search": [
                    {
                        "minimum_magnitude": 7.0,
                        "eligible_catalog_events": 4,
                        "attempted_events": 3,
                        "outcome": "eligible_events_exhausted",
                    }
                ],
                "events": events,
            }
        )
    )
    workflow = tmp_path / "rayleigh.json"
    workflow.write_text(
        json.dumps(
            {
                "execution_status": "complete",
                "options": {"band": [0.03, 0.07]},
                "entries": {
                    "k": {
                        "seed_id": "YV.RR99.00.BDH",
                        "status": "unresolved",
                        "epoch": epoch,
                        "update": {"rayleigh_result_file": str(comparison)},
                    }
                },
            }
        )
    )
    trial_root = tmp_path / "trial"
    (trial_root / "YV_RR99").mkdir(parents=True)
    (trial_root / "YV_RR99" / "trial.json").write_text(
        json.dumps(
            {
                "network": "YV",
                "station": "RR99",
                "baseline_result": str(comparison),
                "positive_real_pressure_poles_rad_s": [0.1026],
                "summaries": {"pole_absolute": {"candidate_events": 0}},
                "events": trial_events,
            }
        )
    )
    data = tmp_path / "data" / "RR99"
    for code in ("BDH", "BHZ"):
        (data / f"{code}.D").mkdir(parents=True)
        (data / f"{code}.D" / f"YV.RR99.00.{code}.D.2012.276").write_text("")
    inputs = [workflow, comparison, *(Path(p) for p in snapshots)]
    inputs.extend(run_dir.glob("*_spectra.npz"))
    (trial_root / "experiment.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "protected_inputs_unchanged": True,
                "source_sha256": {str(p.resolve()): sha256(p) for p in inputs},
            }
        )
    )
    return workflow, trial_root, tmp_path / "data"


def diagnosis_args(tmp_path):
    from argparse import Namespace

    workflow, trial_root, data = synthetic_saved_products(tmp_path)
    return Namespace(
        workflow=str(workflow),
        trial_root=str(trial_root),
        expect_trial_sha=None,
        microseism_run=None,
        data_root=str(data),
        output_dir=str(tmp_path / "out"),
    )


def test_snapshot_refuses_wrong_network(tmp_path):
    from obspy import read_inventory

    path = tmp_path / "wrong-network.xml"
    write_snapshot(path, "BDH", [-0.1], [0j], "PA")
    inventory = read_inventory(str(path))
    inventory[0].code = "XX"
    inventory.write(str(path), "STATIONXML")
    with pytest.raises(ValueError, match="exactly one"):
        snapshot_channel(path, "YV.RR99.00.BDH")


def test_snapshot_refuses_wrong_epoch(tmp_path):
    path = tmp_path / "snapshot.xml"
    write_snapshot(path, "BDH", [-0.1], [0j], "PA")
    epoch = {"channel_start_time": "2012-10-02", "channel_end_time": "2012-10-03"}
    with pytest.raises(ValueError, match="epoch mismatch"):
        snapshot_channel(path, "YV.RR99.00.BDH", epoch)


def test_run_refuses_changed_frozen_comparison(tmp_path):
    args = diagnosis_args(tmp_path)
    comparison = tmp_path / "run" / "comparison.json"
    with comparison.open("a") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="changed since"):
        run(args)
    assert not (tmp_path / "out").exists()


def test_unknown_event_outcome_is_not_mislabelled_as_missing_input():
    with pytest.raises(ValueError, match="Unknown saved event"):
        event_category({"status": "new_scientific_outcome"})


def test_invalid_spectral_grid_refused(tmp_path):
    from scripts.diagnose_yv_rayleigh import load_event

    path = tmp_path / "bad.npz"
    np.savez(
        path,
        frequency_hz=[0.03, 0.03, 0.05],
        phase_deg=[0, 0, 0],
        coherence_squared=[0.99] * 3,
        psd_ratio_pa=[1, 1, 1],
        expected_pa=[1, 1, 1],
    )
    with pytest.raises(ValueError, match="Invalid or mismatched"):
        load_event({"files": {"spectra": str(path)}})


def test_microseism_same_ids_but_wrong_epoch_refused(tmp_path, monkeypatch):
    import scripts.diagnose_yv_rayleigh as module

    args = diagnosis_args(tmp_path)
    micro = tmp_path / "micro"
    micro.mkdir()
    (micro / "summary.json").write_text("{}")
    args.microseism_run = str(micro)
    workflow = json.loads(Path(args.workflow).read_text())
    epoch = dict(workflow["entries"]["k"]["epoch"])
    wrong_epoch = dict(epoch, channel_start_time="2012-10-02T00:00:00.000000Z")
    paired = {
        "pressure_seed_id": "YV.RR99.00.BDH",
        "vertical_seed_id": "YV.RR99.00.BHZ",
        "pressure_epoch": wrong_epoch,
        "vertical_epoch": dict(epoch, Seed_id="YV.RR99.00.BHZ"),
    }
    monkeypatch.setattr(module, "microseism_bins", lambda *_: {"RR99": paired})
    with pytest.raises(ValueError, match="Microseism P/Z identity mismatch"):
        run(args)
    assert not Path(args.output_dir).exists()


def test_run_refuses_trial_for_different_baseline(tmp_path):
    args = diagnosis_args(tmp_path)
    path = next((tmp_path / "trial").glob("*/trial.json"))
    trial = json.loads(path.read_text())
    trial["baseline_result"] = str(tmp_path / "unrelated.json")
    path.write_text(json.dumps(trial))
    with pytest.raises(ValueError, match="baseline"):
        run(args)
    assert not (tmp_path / "out").exists()


def test_run_detects_input_change_before_writing(tmp_path, monkeypatch):
    import scripts.diagnose_yv_rayleigh as module

    args = diagnosis_args(tmp_path)
    original = module.station_record

    def changed_input(*positional, **keywords):
        record = original(*positional, **keywords)
        with Path(args.workflow).open("a") as stream:
            stream.write("\n")
        return record

    monkeypatch.setattr(module, "station_record", changed_input)
    with pytest.raises(ValueError, match="changed"):
        run(args)
    assert not (tmp_path / "out").exists()


def test_end_to_end_saved_products_only_and_refuses_existing_output(tmp_path):
    workflow, trial_root, data = synthetic_saved_products(tmp_path)

    class Args:
        pass

    args = Args()
    args.workflow, args.trial_root = str(workflow), str(trial_root)
    args.expect_trial_sha, args.microseism_run = None, None
    args.data_root, args.output_dir = str(data), str(tmp_path / "out")
    before = {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
    result = run(args)
    after = {p: p.read_bytes() for p in before}
    assert before == after  # inputs untouched
    record = result["stations"][0]
    assert record["event_outcomes"] == {"phase_review": 3, "selection_excluded_24h": 1}
    assert record["events_attempted"] == 3
    assert record["trial_reflection_reproduction_max_abs_deg"] < 1e-9
    assert len(record["response"]["vertical_unpaired_roots"]) == 4
    assert record["rayleigh"]["H2"]["phase_rms_deg"] < 1e-6
    assert record["rayleigh"]["H2"]["events_passing_phase_screen_only"] == 3
    assert record["rayleigh"]["H1"]["events_passing_phase_screen_only"] == 0
    assert abs(record["rayleigh"]["H2"]["median_log_ratio_slope"]) < 1e-9
    assert record["local_coverage"]["paired_epoch_day_files"] == 1
    assert record["local_coverage"]["epoch_days"] == 3
    assert result["gain_fitted"] is False and result["waveforms_read"] is False
    assert result["protected_inputs_unchanged"] is True
    assert (
        str((tmp_path / "run" / "e0_spectra.npz").resolve()) in result["input_sha256"]
    )
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == [
        "evidence.json",
        "stations.csv",
    ]
    with pytest.raises(SystemExit):
        run(args)
