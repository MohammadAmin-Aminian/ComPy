"""Saved-spectrum response pilot: trial replay, controls and refusals."""

import json
from types import SimpleNamespace

import numpy as np
import pytest
from obspy import UTCDateTime
from obspy.core.inventory import Channel, Inventory, Network, Site, Station
from obspy.core.inventory.response import Response

from compy_pressure.calibration_experiment import evaluate, pole_factor
from scripts.diagnose_yv_rayleigh import pz_value, sha256, stage_ratio
from scripts.run_response_pilot import control_checks, run

P = -0.074016 + 0.07347j
Q = -502.65 + 596.9j
GOOD_Z = [P, P.conjugate(), Q, Q.conjugate()]
BAD_Z = [P, P, Q, Q]
POLE = 0.1026
START = UTCDateTime("2012-10-20T18:27:45")
END = UTCDateTime("2013-06-29T20:55:33")
SETTINGS = {
    "pressure_sign": -1,
    "tilt_method": "none",
    "band_hz": [0.03, 0.07],
    "coherence_squared_threshold": 0.95,
    "phase_limit_deg": 20.0,
    "minimum_candidate_events": 3,
    "gain_grid_bounds": [0.1, 3.0],
}


def channel(code, poles, zeros, units, a0):
    response = Response.from_paz(
        zeros=zeros,
        poles=poles,
        stage_gain=1.0,
        input_units=units,
        output_units="V",
        normalization_factor=a0,
    )
    return Channel(
        code=code,
        location_code="00",
        latitude=0,
        longitude=0,
        elevation=-4000,
        depth=0,
        sample_rate=50,
        start_date=START,
        end_date=END,
        response=response,
    )


def a0(poles, zeros):
    return 1.0 / abs(pz_value(np.array([1.0]), poles, zeros)[0])


def write(path, channels):
    station = Station("RR99", 0, 0, -4000, channels=channels, site=Site("RR99"))
    Inventory([Network("YV", stations=[station])]).write(str(path), "STATIONXML")


def saved_products(tmp_path, g_true=0.8, events=3):
    """RR56-like saved run whose physical P/Z relation needs H2 to appear."""
    f = np.linspace(0.03, 0.07, 9)
    # Like RR56: pressure A0 from its stored pole, vertical A0 from paired roots.
    pressure = channel("BDH", [POLE + 0j], [0j], "PA", a0([POLE + 0j], [0j]))
    vertical = channel("BHZ", BAD_Z, [0j, 0j], "M/S", a0(GOOD_Z, [0j, 0j]))
    source = tmp_path / "YV.xml"
    write(source, [pressure, vertical])
    keys, snapshots = [], {}
    for item in (pressure, vertical):
        seed_id = f"YV.RR99.00.{item.code}"
        folder = tmp_path / f"{seed_id}__20121020T182745__abc"
        folder.mkdir()
        write(folder / "inventory.xml", [item])
        snapshots[str(folder / "inventory.xml")] = sha256(folder / "inventory.xml")
        keys.append(
            {
                "Seed_id": seed_id,
                "channel_start_time": str(START),
                "channel_end_time": str(END),
                "inventory_file": str(source),
            }
        )
    # Nominal processing of a true H2 world: P_saved = P_true / rp etc.
    rp = stage_ratio(f, [([POLE + 0j], [0j])], [([-POLE + 0j], [0j])])
    rz = stage_ratio(f, [(BAD_Z, [0j, 0j])], [(GOOD_Z, [0j, 0j])])
    expected = np.full(9, 2.0)
    comparison_events, trial_events = [], []
    for i in range(events):
        path = tmp_path / f"e{i}_spectra.npz"
        spectra = {
            "frequency_hz": f,
            "coherence_squared": np.full(9, 0.99),
            "cross_pa": (1 + 0.1 * i) / (rp * np.conj(rz)),
            "psd_ratio_pa": expected / g_true * abs(rz) / abs(rp),
            "expected_pa": expected,
            "segment_count": 9,
        }
        np.savez(path, **spectra)
        comparison_events.append(
            {
                "id": f"e{i}",
                "status": "polarity_review",
                "files": {"spectra": str(path)},
                "tilt": {"vertical_component": "YV.RR99.00.BHZ"},
            }
        )
        # The completed trial's own two variants, computed its way.
        nominal = evaluate(spectra, SETTINGS, False, np.ones(9), "absolute")
        pole = evaluate(spectra, SETTINGS, False, pole_factor(f, [POLE]), "absolute")
        trial_events.append(
            {
                "id": f"e{i}",
                "variants": {"nominal_absolute": nominal, "pole_absolute": pole},
            }
        )
    comparison_events.append({"id": "x", "status": "excluded_cluster"})
    trial_events.append({"id": "x", "status": "excluded_cluster"})
    comparison = tmp_path / "comparison.json"
    comparison.write_text(
        json.dumps(
            {
                "settings": SETTINGS,
                "epoch_keys": keys,
                "source_inventory_sha256": {str(source): sha256(source)},
                "response_snapshot_sha256": snapshots,
                "events": comparison_events,
            }
        )
    )
    workflow = tmp_path / "rayleigh.json"
    entry = {
        "seed_id": "YV.RR99.00.BDH",
        "update": {"rayleigh_result_file": str(comparison)},
    }
    workflow.write_text(json.dumps({"entries": {"k": entry}}))
    trial_root = tmp_path / "trial"
    (trial_root / "YV_RR99").mkdir(parents=True)
    empty = {"median_pressure_multiplier": None}
    (trial_root / "YV_RR99" / "trial.json").write_text(
        json.dumps(
            {
                "network": "YV",
                "station": "RR99",
                "baseline_result": str(comparison),
                "summaries": {"nominal_absolute": empty, "pole_absolute": empty},
                "events": trial_events,
            }
        )
    )
    (trial_root / "experiment.json").write_text("{}")
    return SimpleNamespace(
        workflow=str(workflow),
        workflow_sha256=sha256(workflow),
        trial_root=str(trial_root),
        trial_sha256=sha256(trial_root / "experiment.json"),
        out=str(tmp_path / "out"),
    )


PLAN = {"RR99": {"hypotheses": ["H2"], "role": "pilot", "label": "test"}}


def test_pilot_replays_trial_then_recovers_the_true_gain_under_h2(tmp_path):
    args = saved_products(tmp_path)
    report = run(args, PLAN)
    summaries = report["stations"]["RR99"]["summaries"]
    assert summaries["H0"]["status_counts"] == {"polarity_review": 3}
    assert summaries["H1"]["development_g"] is None
    assert summaries["H2"]["candidate_events"] == 3
    assert summaries["H2"]["development_g"] == pytest.approx(0.8)
    assert summaries["H2"]["phase_rms_deg"] == pytest.approx(0.0, abs=1e-9)
    effect = report["stations"]["RR99"]["if_a0_recomputed_from_hypothesis_roots"]
    # The recorded vertical A0 belongs to the paired roots, not the stored ones.
    assert effect["H0"]["vertical_a0_ratio"] < 0.97
    assert effect["H2"]["g_factor"] == pytest.approx(1.0, abs=1e-12)
    assert effect["H1"]["pressure_a0_ratio"] == pytest.approx(1.0, abs=1e-12)
    assert (tmp_path / "out" / "pilot.json").is_file()
    assert (tmp_path / "out" / "stations.csv").read_text().count("\n") == 4
    with pytest.raises(FileExistsError):
        run(args, PLAN)


def test_pilot_stops_when_the_trial_outcome_does_not_replay(tmp_path):
    args = saved_products(tmp_path)
    trial_path = tmp_path / "trial" / "YV_RR99" / "trial.json"
    trial = json.loads(trial_path.read_text())
    trial["events"][0]["variants"]["pole_absolute"]["status"] = "candidate"
    trial_path.write_text(json.dumps(trial))
    with pytest.raises(ValueError, match="does not replay"):
        run(args, PLAN)
    assert not (tmp_path / "out").exists()


def test_pilot_refuses_changed_inputs(tmp_path):
    args = saved_products(tmp_path)
    snapshot = next(tmp_path.glob("YV.RR99.00.BHZ__*/inventory.xml"))
    snapshot.write_text(snapshot.read_text() + "\n")
    with pytest.raises(ValueError, match="recorded value"):
        run(args, PLAN)
    (tmp_path / "second").mkdir()
    args = saved_products(tmp_path / "second")
    args.trial_sha256 = "0" * 64
    with pytest.raises(ValueError, match="recorded value"):
        run(args, PLAN)


def test_identity_and_negative_control_checks():
    def report(role, hypotheses, g, h1_rms=None, rms=None):
        summaries = {"H1": {"development_g": None, "phase_rms_deg": h1_rms}}
        for name in hypotheses:
            summaries[name] = {"development_g": g, "phase_rms_deg": rms}
        return {"role": role, "hypotheses": hypotheses, "summaries": summaries}

    plan = {
        "RR52": {"expected_g": 0.6787278140712512, "exact": True},
        "RR16": {"expected_g": 1.0, "exact": False},
    }
    exact = control_checks(
        {"RR52": report("identity_control", ["H2"], 0.6787278140712512)}, plan
    )
    assert exact[0]["passed"]
    off = control_checks(
        {"RR52": report("identity_control", ["H2"], 0.6787278140712513)}, plan
    )
    assert not off[0]["passed"]
    close = control_checks(
        {"RR16": report("identity_control", ["H2"], 1.0 + 1e-12)}, plan
    )
    assert close[0]["passed"]
    worse = control_checks(
        {"RR10": report("negative_control", ["H4"], None, 10.0, 40.0)}, plan
    )
    better = control_checks(
        {"RR10": report("negative_control", ["H4"], None, 10.0, 5.0)}, plan
    )
    assert worse[0]["passed"] and not better[0]["passed"]
