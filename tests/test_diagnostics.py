"""Behavioral checks for independent-chain diagnostics and saved run integrity."""

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from compy_diagnostics import save_run, summarize_chains


def test_diagnostics_detect_location_scale_and_serial_correlation():
    iid = np.random.default_rng(42).normal(size=(4, 2000, 1))
    baseline = summarize_chains(iid, ["x"])["parameters"]["x"]
    assert baseline["rhat"] < 1.01
    assert baseline["ess_bulk"] > 4000
    shifted = iid.copy()
    shifted[0] += 3
    assert summarize_chains(shifted, ["x"])["parameters"]["x"]["rhat"] > 1.1
    scaled = iid.copy()
    scaled[0] *= 5
    assert summarize_chains(scaled, ["x"])["parameters"]["x"]["rhat"] > 1.1
    correlated = iid.copy()
    for draw in range(1, 2000):
        correlated[:, draw] += 0.95 * correlated[:, draw - 1]
    assert summarize_chains(correlated, ["x"])["parameters"]["x"]["ess_bulk"] < 1000


def test_burnin_exclusion_and_stuck_chain_warning():
    draws = np.random.default_rng(2).normal(size=(4, 100, 1))
    draws[:, :20] += 1000
    report = summarize_chains(draws, ["x"], burnin=20)
    assert report["parameters"] == summarize_chains(draws[:, 20:], ["x"])["parameters"]
    draws[0, 20:] = 0
    report = summarize_chains(draws, ["x"], burnin=20)
    assert report["parameters"]["x"]["rhat"] is None
    assert any("constant" in message for message in report["warnings"])
    json.dumps(report, allow_nan=False)
    identical = summarize_chains(np.repeat(draws[1:2], 4, axis=0), ["x"], burnin=20)
    assert any("Identical" in message for message in identical["warnings"])


@pytest.mark.parametrize(
    "draws,names,burnin",
    [
        (np.zeros((1, 10, 1)), ["x"], 0),
        (np.zeros((4, 10, 1)), ["x"], 7),
        (np.zeros((4, 10, 1)), ["x"], True),
        (np.full((4, 10, 1), np.nan), ["x"], 0),
        (np.zeros((4, 10, 2)), ["x", "x"], 0),
    ],
)
def test_invalid_diagnostics_rejected(draws, names, burnin):
    with pytest.raises(ValueError):
        summarize_chains(draws, names, burnin=burnin)


def test_run_roundtrip_provenance_checksums_and_overwrite_protection(tmp_path):
    output = tmp_path / "run"
    data = np.arange(12.0).reshape(3, 4)
    save_run(output, {"models": data}, settings={"seed": 3, "units": "m/s"})
    with np.load(output / "arrays.npz", allow_pickle=False) as saved:
        np.testing.assert_array_equal(saved["models"], data)
    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["settings"]["seed"] == 3
    assert metadata["arrays"]["models"]["shape"] == [3, 4]
    assert metadata["dependencies"]["numpy"] == np.__version__
    assert len(metadata["source_sha256"]["compy_diagnostics.py"]) == 64
    manifest = json.loads((output / "manifest.json").read_text())
    for name, expected in manifest.items():
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == expected
    with pytest.raises(FileExistsError):
        save_run(output, {"models": data + 1}, settings={})
    assert json.loads((output / "manifest.json").read_text()) == manifest


def test_invalid_metadata_and_failed_write_leave_no_partial_run(tmp_path, monkeypatch):
    output = tmp_path / "run"
    with pytest.raises(ValueError):
        save_run(output, {"x": [1.0]}, settings={"bad": float("nan")})
    assert not output.exists()
    with pytest.raises(ValueError):
        save_run(output, {"x": np.array([object()])}, settings={})
    assert not output.exists()

    def fail(*args, **kwargs):
        raise OSError("simulated disk failure")

    monkeypatch.setattr(np, "savez_compressed", fail)
    with pytest.raises(OSError):
        save_run(output, {"x": [1.0]}, settings={})
    assert not output.exists()
    assert list(tmp_path.iterdir()) == []


def test_offline_example_is_reproducible_and_preserves_forward_consistency(tmp_path):
    example = Path(__file__).resolve().parents[1] / "_Example/reproducible_inversion.py"
    spec = importlib.util.spec_from_file_location("reproducible_example", example)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    first, second = tmp_path / "first", tmp_path / "second"
    module.run(first, iterations=30, burnin=10, chains=2)
    module.run(second, iterations=30, burnin=10, chains=2)
    import inv_compy as inv

    with (
        np.load(first / "arrays.npz", allow_pickle=False) as a,
        np.load(second / "arrays.npz", allow_pickle=False) as b,
    ):
        for name in a.files:
            np.testing.assert_array_equal(a[name], b[name])
        np.testing.assert_allclose(
            a["predictions"][0, -1],
            inv.calc_norm_compliance(4000, a["frequency"], a["models"][0, :, :, -1]),
        )
    assert (first / "diagnostics.png").read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
