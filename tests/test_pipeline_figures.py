"""Pipeline figure stage: results without a Rayleigh-wave calibration."""

import importlib
import sys
from pathlib import Path

PIPELINE = Path(__file__).resolve().parents[1] / "pipeline"


def _f_figures(tmp_path, monkeypatch):
    monkeypatch.setenv("COMPY_RUN_DIR", str(tmp_path))
    monkeypatch.delenv("COMPY_PIPELINE_CONFIG", raising=False)
    monkeypatch.syspath_prepend(str(PIPELINE))
    for name in ("config", "e_invert", "f_figures"):
        sys.modules.pop(name, None)
    return importlib.import_module("f_figures")


def test_result_label_follows_stage_e_order(tmp_path, monkeypatch):
    f = _f_figures(tmp_path, monkeypatch)
    both = {"results": {"compy": {"n_windows": 600}, "nominal": {"n_windows": 600}}}
    assert f.result_label(both) == "compy"
    # Calibration failed (no M >= 7 event): stage C falls back to the nominal gain.
    uncalibrated = {"results": {"nominal": {"n_windows": 430}}}
    assert f.result_label(uncalibrated) == "nominal"
    assert f.result_label({"results": {"compy": {"attempts": {}}}}) is None
    assert f.result_label(None) is None
