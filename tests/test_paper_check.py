"""Formulas checked against the papers in November 2026 (docs/paper_check_2026-11)."""

import importlib
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from compy.compliance import gravitational_attraction
from compy.inversion import auto, dram
from compy.numerics import GRAVITATIONAL_CONSTANT, GRAVITY, wavenumber
from compy.response_check import predicted_residual
from compy_pressure.rayleigh_calibration import water_ratio


@pytest.fixture
def stage_c(tmp_path, monkeypatch):
    monkeypatch.setenv("COMPY_RUN_DIR", str(tmp_path))
    monkeypatch.delenv("COMPY_PIPELINE_CONFIG", raising=False)
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "pipeline"))
    for name in ("config", "c_compliance"):
        sys.modules.pop(name, None)
    yield importlib.import_module("c_compliance")
    for name in ("config", "c_compliance"):
        sys.modules.pop(name, None)


def test_ruan_sediment_profile_follows_its_equation():
    # Ruan et al. (2014), eq. 5 and Table 2: 0.1 km/s at the seafloor.
    assert dram.ruan_vs(0.0) == pytest.approx(100.0)
    h = 0.3
    expected = 1000 * (0.02 * h**2 + 1.27 * h + 0.048) / (h + 0.48)
    assert dram.ruan_vs(300.0) == pytest.approx(expected)
    z = np.array([20.0, 50.0, 100.0, 300.0, 500.0])
    ratio = dram.ruan_vs(z) / dram.hamilton_vs(z)
    assert np.all(np.diff(dram.ruan_vs(z)) > 0)
    # inside the hard bounds of the sediment Vs scale (0.4-2.5) and within 1.5 sigma (0.3)
    assert np.all(np.abs(np.log(ratio)) < 0.45)


def test_seafloor_model_accepts_the_ruan_profile():
    reference = np.array([[2000.0, 2700.0, 5500.0, 3000.0], [5000.0, 2950.0, 6800.0, 3700.0],
                          [0.0, 3300.0, 8000.0, 4500.0]])  # fmt: skip
    model = auto.SeafloorModel(reference, sediment_vs=dram.ruan_vs)
    layers = model.layers(model.reference_theta(sediment=200.0))
    sediment = layers[: model.sediment_sublayers]
    assert sediment[:, 0].sum() == pytest.approx(200.0)
    assert sediment[0, 3] < dram.hamilton_vs(
        1.0
    )  # slower than Hamilton at the seafloor
    assert np.all(layers[:, 2] >= np.sqrt(2) * layers[:, 3])


def test_wave_attraction_is_crawford_1998_with_the_exponent_k_h():
    # a_w = 2 pi G rho exp(-kH) h_w (eq. 6) with h_w = p cosh(kH) / (rho g) (eq. 3).
    depth, f, rho = 4700.0, np.array([0.004, 0.008, 0.014]), 1025.0
    k = wavenumber(2 * np.pi * f, depth)
    ratio, _, h_w = gravitational_attraction(np.ones(3), depth, f, pw=rho)
    np.testing.assert_allclose(h_w, np.cosh(k * depth) / (rho * GRAVITY), rtol=1e-12)
    np.testing.assert_allclose(
        ratio,
        2 * np.pi * GRAVITATIONAL_CONSTANT * rho * np.exp(-k * depth) * h_w,
        rtol=1e-12,
    )


def test_rayleigh_pressure_reduces_to_rho_h_and_departs_from_it_in_deep_water():
    # An et al. (2017), eq. 6: p = rho h a where the water column is short against the wave.
    rho, c0 = 1028.0, 1500.0
    shallow = water_ratio([0.02], [4000.0], 350.0, rho, c0)[0]
    assert shallow == pytest.approx(rho * 350.0, rel=1e-3)
    # Deng et al. (2022), eq. 5 (Zha & Webb 2016, A5): rho h tan(x)/x with compressible water.
    depth, f, c = 4700.0, 0.06, 4000.0
    x = 2 * np.pi * f * depth * np.sqrt(c0**-2 - c**-2)
    deep = water_ratio([f], [c], depth, rho, c0)[0]
    assert deep == pytest.approx(rho * depth * np.tan(x) / x, rel=1e-12)
    assert (
        deep / (rho * depth) > 1.5
    )  # rho h alone is not usable at this depth and frequency


def test_stage_c_response_check_recovers_a_lower_leak_corner(stage_c):
    # Metadata: one pole at 2 mHz. Truth: 1 mHz (a longer leak time constant at depth).
    stages = [(np.array([0j]), np.array([-2 * np.pi * 0.002 + 0j]), True)]
    f = np.arange(0.005, 0.0151, 0.00122)
    phase = predicted_residual(stages, f, 0.5)
    result = SimpleNamespace(freqs=f, band=(f[0], f[-1]), phase_deg=phase)
    check = stage_c.response_check(result, stages)
    assert check["k"] == pytest.approx(0.5, rel=0.03)
    assert check["fitted_corner_hz"] == pytest.approx(0.001, rel=0.03)
    assert check["metadata_corner_hz"] == pytest.approx(0.002)
    assert check["polarity_deg"] == 0.0 and check["applied"] is False
    # deconvolving too high a corner inflates low-frequency pressure: compliance is low,
    # more at the low end of the band
    assert check["bias_percent_band_low"] < check["bias_percent_band_high"] < 0
    assert stage_c.response_check(result, []) is None


def test_stage_e_leak_correction_needs_a_usable_check(tmp_path, monkeypatch):
    monkeypatch.setenv("COMPY_RUN_DIR", str(tmp_path))
    monkeypatch.delenv("COMPY_PIPELINE_CONFIG", raising=False)
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "pipeline"))
    for name in ("config", "e_invert"):
        sys.modules.pop(name, None)
    stage_e = importlib.import_module("e_invert")
    f, eta = np.array([0.006, 0.01]), np.array([2e-11, 1.5e-11])
    try:
        assert stage_e.leak_correction("STA", None, f, eta) is None
        assert (
            stage_e.leak_correction("STA", {"k": 0.02, "at_grid_edge": True}, f, eta)
            is None
        )
    finally:
        for name in ("config", "e_invert"):
            sys.modules.pop(name, None)
