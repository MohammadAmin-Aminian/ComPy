"""Pipeline: one gain-selection rule for stages E, F and G; figures without calibration."""

import importlib
import sys
from pathlib import Path

PIPELINE = Path(__file__).resolve().parents[1] / "pipeline"


def _pipeline(tmp_path, monkeypatch, *names):
    monkeypatch.setenv("COMPY_RUN_DIR", str(tmp_path))
    monkeypatch.delenv("COMPY_PIPELINE_CONFIG", raising=False)
    monkeypatch.syspath_prepend(str(PIPELINE))
    for name in ("config", "e_invert", "f_figures", "g_report"):
        sys.modules.pop(name, None)
    return [importlib.import_module(n) for n in names]


def test_result_label_order(tmp_path, monkeypatch):
    (e_invert,) = _pipeline(tmp_path, monkeypatch, "e_invert")
    both = {"results": {"compy": {"n_windows": 600}, "nominal": {"n_windows": 600}}}
    assert e_invert.result_label(both) == "compy"
    # Calibration failed (no M >= 7 event): stage C falls back to the nominal gain.
    assert (
        e_invert.result_label({"results": {"nominal": {"n_windows": 430}}}) == "nominal"
    )
    assert e_invert.result_label({"results": {"compy": {"attempts": {}}}}) is None
    assert e_invert.result_label(None) is None


def test_stages_share_one_rule(tmp_path, monkeypatch):
    e_invert, f_figures, g_report = _pipeline(
        tmp_path, monkeypatch, "e_invert", "f_figures", "g_report"
    )
    assert f_figures.result_label is e_invert.result_label
    assert g_report.result_label is e_invert.result_label


def test_reference_figures_skip_missing_stage_output(tmp_path, monkeypatch):
    (f_figures,) = _pipeline(tmp_path, monkeypatch, "f_figures")
    # A station that stage B rejected has no clean segments and no tilt table.
    assert f_figures.fig03_psd("XX01") is None
    assert f_figures.fig04_tilt("XX01") is None


def test_compliance_module_imports_no_plotting_or_fdsn():
    import subprocess

    code = (
        "import sys, compy.compliance; "
        "bad = [m for m in ('matplotlib.pyplot', 'obspy.clients.fdsn') if m in sys.modules]; "
        "assert not bad, bad"
    )
    root = Path(__file__).resolve().parents[1]
    subprocess.run([sys.executable, "-c", code], check=True, cwd=root)
