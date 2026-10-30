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
