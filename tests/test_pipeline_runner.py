"""Failure propagation and station scoping across actual child processes."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

PIPELINE = Path(__file__).resolve().parents[1] / "pipeline"


@pytest.fixture
def runner(monkeypatch, tmp_path):
    monkeypatch.setenv("COMPY_RUN_DIR", str(tmp_path / "run"))
    monkeypatch.delenv("COMPY_PIPELINE_CONFIG", raising=False)
    monkeypatch.syspath_prepend(str(PIPELINE))
    monkeypatch.delitem(sys.modules, "config", raising=False)
    spec = importlib.util.spec_from_file_location(
        "run_pipeline", PIPELINE / "run_pipeline.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    scripts = tmp_path / "scripts"
    scripts.mkdir()
    monkeypatch.setattr(module, "HERE", scripts)
    yield module
    sys.modules.pop("config", None)


def test_failed_stage_stops_the_chain_and_keeps_its_log(runner):
    (runner.HERE / "b_clean.py").write_text(
        "print('bad waveform'); raise SystemExit(7)"
    )
    (runner.HERE / "d_calibrate.py").write_text("raise AssertionError('must not run')")
    assert runner.main(["OBS1", "--start-at", "B", "--stop-after", "D"]) == 1
    manifest = next((runner.config.ROOT / "logs").glob("pipeline_*.json"))
    record = json.loads(manifest.read_text())
    assert record["status"] == "failed"
    assert len(record["steps"]) == 1
    assert record["steps"][0]["returncode"] == 7
    assert "bad waveform" in Path(record["steps"][0]["log"]).read_text()


def test_all_children_receive_the_selected_stations_and_absolute_config(
    runner, monkeypatch
):
    override = runner.HERE / "custom.py"
    override.write_text("NETWORK = 'XX'\n")
    monkeypatch.setenv("COMPY_PIPELINE_CONFIG", str(override))
    (runner.HERE / "g_report.py").write_text(
        "import os\nfrom pathlib import Path\n"
        "assert os.environ['COMPY_STATIONS'] == 'OBS1,OBS2'\n"
        "assert Path(os.environ['COMPY_PIPELINE_CONFIG']).is_absolute()\n"
        "print('selected scope received')\n"
    )
    assert runner.main(["OBS1", "OBS2", "--start-at", "G"]) == 0
    manifest = next((runner.config.ROOT / "logs").glob("pipeline_*.json"))
    record = json.loads(manifest.read_text())
    assert record["status"] == "complete"
    assert record["stations"] == ["OBS1", "OBS2"]
    assert record["steps"][0]["returncode"] == 0


def test_compliance_ignores_a_gain_from_a_failed_calibration(runner):
    spec = importlib.util.spec_from_file_location(
        "c_compliance", PIPELINE / "c_compliance.py"
    )
    stage = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stage)
    folder = stage.WORK / "OBS1" / "calibration"
    folder.mkdir(parents=True)
    (folder / "compy_gain.json").write_text(
        json.dumps({"gain": 0.7, "error": "bad spectra"})
    )
    assert stage.gains_for("OBS1") == {"nominal": 1.0}


def test_compliance_window_levels_match_the_band_spectra(runner):
    import types

    import numpy as np
    from obspy import Stream, Trace, UTCDateTime

    from compy.compliance import spectra

    spec = importlib.util.spec_from_file_location(
        "c_compliance", PIPELINE / "c_compliance.py"
    )
    stage = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stage)
    assert stage.SAVE_WINDOW_LEVELS is True
    rng = np.random.default_rng(3)
    windows = []
    for i in range(3):
        t0 = UTCDateTime(2020, 1, 1) + i * stage.WINDOW_S
        p = rng.normal(0, 1.0, 2625)
        z = 1e-11 * p + rng.normal(0, 1e-12, 2625)
        windows.append(Stream([
            Trace(z, header={"channel": stage.Z, "sampling_rate": 1.25, "starttime": t0}),
            Trace(p, header={"channel": stage.P, "sampling_rate": 1.25, "starttime": t0}),
        ]))  # fmt: skip
    result = types.SimpleNamespace(
        window_starts=[windows[0][0].stats.starttime, windows[2][0].stats.starttime],
        band=(0.007, 0.02),
    )
    levels = stage.window_levels(windows, result, gain=0.8)
    assert levels.shape == (2, 4)
    f, zz, pp, zp, _ = spectra(windows[2], stage.NSEG)
    sel = (f >= 0.007) & (f <= 0.02)
    expected = np.log10(np.median(np.abs(0.8 * zp[sel]) / (0.64 * pp[sel])))
    assert levels[1, 0] == float(windows[2][0].stats.starttime)
    assert levels[1, 1] == pytest.approx(expected)
    assert levels[1, 2] == pytest.approx(10 * np.log10(np.mean(0.64 * pp[sel])))


def test_calibration_sign_from_the_compliance_phase(runner):
    import numpy as np

    pytest.importorskip("compy_pressure.rayleigh_network")
    spec = importlib.util.spec_from_file_location(
        "d_calibrate", PIPELINE / "d_calibrate.py"
    )
    stage = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stage)
    assert stage.pressure_sign("OBS1") == (None, "compy_pressure default")
    assert stage.pressure_sign("OBS1", 1) == (1, "config")
    assert stage.pressure_sign("OBS1", "phase")[0] is None  # no stage C yet
    folder = stage.RESULTS / "OBS1"
    folder.mkdir(parents=True)
    f = np.linspace(0.005, 0.03, 20)
    np.savez(folder / "OBS1_compliance_nominal.npz", freqs=f, band=np.array([0.007, 0.02]),
             phase_deg=np.full(20, 176.0))  # fmt: skip
    sign, source = stage.pressure_sign("OBS1", "phase")
    assert sign == 1 and "OBS1_compliance_nominal.npz" in source
