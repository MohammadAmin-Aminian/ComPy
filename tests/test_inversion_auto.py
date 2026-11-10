"""Tests of the ComPy 3 automatic inversion (compy.inversion.auto, .errors).

Each test states the finding of docs/inversion_overhaul_2026-11 that it protects.
"""

import numpy as np
import pytest

from compy.inversion import auto, dram
from compy.inversion.errors import (
    compliance_from_sums,
    default_covariance,
    log_frequency_ramp,
    stack_covariance,
)
from compy.numerics import wavenumber

REFERENCE = np.array([
    [1170.0, 2400.0, 4700.0, 2540.0],
    [1170.0, 2760.0, 6300.0, 3590.0],
    [4560.0, 2960.0, 6900.0, 3940.0],
    [0.0, 3190.0, 7730.0, 4310.0],
])  # fmt: skip
DEPTH = 3000.0
FREQS = np.linspace(0.007, 0.017, 9)


def synthetic_windows(rng, n=600, slow=0.0, days=120.0):
    """Windows of P and Z spectra with a known transfer function and vertical noise."""
    f = FREQS
    transfer = 8e-8 * (
        1 + 20 * f
    )  # |Z/P| in m/Pa: realistic size, so it outweighs the attraction
    starts = np.sort(rng.uniform(0, days * 86400, n))
    level = 1 + slow * np.sin(2 * np.pi * starts / (60 * 86400))
    pp = rng.gamma(4, 0.25, (n, len(f))) * 1e2
    noise = rng.normal(0, 0.15, (n, len(f))) + 1j * rng.normal(0, 0.15, (n, len(f)))
    zp = pp * transfer * level[:, None] * (1 + noise)
    return f, pp, zp, starts, transfer


# -- errors ---------------------------------------------------------------------
def test_stack_covariance_matches_scatter_of_independent_stacks():
    """The window bootstrap reproduces the scatter of independent realizations."""
    rng = np.random.default_rng(0)
    f, pp, zp, starts, _ = synthetic_windows(rng)
    stat = stack_covariance(f, pp, zp, starts, DEPTH, seed=1)
    etas = []
    for i in range(200):
        f, p2, z2, _, _ = synthetic_windows(np.random.default_rng(100 + i))
        etas.append(compliance_from_sums(f, DEPTH, z2.sum(0), p2.sum(0)))
    truth = np.std(np.array(etas) / np.mean(etas, 0), axis=0)
    assert np.all(np.linalg.eigvalsh(stat["cov_rel"]) > 0)
    assert np.allclose(stat["sd_rel"], truth, rtol=0.35)
    assert stat["slow_level"] < 0.5 * np.median(
        stat["sd_rel"]
    )  # stationary: no slow term


def test_stack_covariance_finds_a_slow_level_change():
    """A level that drifts over weeks appears as the slow level term, not as noise."""
    rng = np.random.default_rng(3)
    f, pp, zp, starts, _ = synthetic_windows(rng, n=3000, slow=0.03, days=360.0)
    stat = stack_covariance(f, pp, zp, starts, DEPTH, seed=1)
    assert stat["block_days"] == 30.0 and stat["n_blocks"] >= 11
    assert 0.003 < stat["slow_level"] < 0.02
    assert stat["slow_level"] > 2 * np.median(stat["sd_rel"])


def test_stack_covariance_rejects_bad_input():
    f, pp, zp, starts, _ = synthetic_windows(np.random.default_rng(0), n=30)
    with pytest.raises(ValueError):
        stack_covariance(f, pp[:10], zp[:10], starts[:10], DEPTH)
    with pytest.raises(ValueError):
        stack_covariance(f, pp, zp[:, :-1], starts, DEPTH)


def test_default_covariance_has_the_measured_structure():
    eta = np.linspace(2e-11, 4e-11, 9)
    sigma = 0.002 * eta
    c = default_covariance(FREQS, eta, sigma)
    sd = np.sqrt(np.diag(c))
    corr = c / np.outer(sd, sd)
    assert np.all(np.linalg.eigvalsh(c) > 0)
    assert np.all(sd > 1.36 * 0.002) and np.all(sd < 0.006)
    assert corr[0, 1] > corr[0, 2] > 0  # neighbouring bins correlated
    assert log_frequency_ramp(FREQS)[[0, -1]] == pytest.approx([-0.5, 0.5])


def test_tilt_and_smooth_terms_of_the_covariance():
    """A tilt of the curve costs little chi-square, a zigzag of the same size much."""
    eta = dram.forward(DEPTH, FREQS, REFERENCE)
    white = (0.002 * eta) ** 2 * np.eye(len(eta))
    data = dram.ComplianceData.build(FREQS, eta, DEPTH, random_cov=white, gain_rel=0, floor_rel=0,
                                     tilt_rel=0.04, smooth_rel=0.005)  # fmt: skip
    assert set(data.components) == {"random", "tilt", "smooth"}
    ramp = log_frequency_ramp(FREQS)
    tilt = eta * (1 + 0.04 * ramp)
    zigzag = eta * (1 + 0.02 * (-1) ** np.arange(len(eta)))
    assert data.chi2(tilt) < 2.0
    assert data.chi2(zigzag) > 100.0
    rel = np.sqrt(np.diag(data.components["tilt"])) / eta
    assert rel[[0, -1]] == pytest.approx([0.02, 0.02])


# -- rock relations ---------------------------------------------------------------
def test_oceanic_vp_is_continuous_and_follows_the_mafic_line():
    vs = np.linspace(300, 4800, 2000)
    vp = auto.oceanic_vp(vs)
    assert np.all(np.diff(vp) > 0)
    assert np.max(np.abs(np.diff(vp))) < 15.0  # no step between the pieces
    assert auto.oceanic_vp(2880.0) == pytest.approx(5250.0)
    assert auto.oceanic_vp(3920.0) == pytest.approx(7250.0)
    ratio = auto.oceanic_vp(np.array([2880.0, 3400.0, 3920.0, 4500.0])) / np.array(
        [2880.0, 3400.0, 3920.0, 4500.0]
    )
    assert np.all((ratio > 1.82) & (ratio < 1.86))
    assert auto.oceanic_vp(800.0) / 800.0 > 2.5  # slow rock has a high Vp/Vs


def test_density_relations():
    assert auto.carlson_raskin_density(5000.0) == pytest.approx(2742.0)
    assert auto.carlson_raskin_density(6900.0) == pytest.approx(
        1000 * (3.81 - 5.99 / 6.9)
    )
    assert 1600 < auto.nafe_drake_density(1510.0) < 1700  # soft sediment, not 2000
    vp = np.linspace(1500, 8200, 2000)
    rho = auto.crust_density(vp)
    step = np.diff(rho)
    # continuous at the switch; the published relation itself steps down 0.7 % at 6.65 km/s
    assert np.max(step) < 3.0 and np.sum(step < 0) == 1 and np.min(step) > -25.0
    assert 6640 < vp[np.argmin(step)] < 6660
    assert rho[0] > 1600 and rho[-1] < 3100


def test_half_space_compliance_is_the_plane_strain_modulus():
    """Quasi-static compliance of a half-space is 1/M, M = mu/(1 - nu) (Crawford 2004, eq. 5)."""
    row = np.array([[0.0, 2900.0, 6800.0, 3800.0]])
    eta = dram.forward(4000.0, np.array([0.004, 0.008]), row)
    assert eta * auto.plane_strain_modulus(row)[0] == pytest.approx(1.0, abs=0.004)


def test_thin_soft_layer_acts_through_its_p_modulus():
    """For kh << 1 sediment adds k h/(rho Vp^2) to compliance, whatever its Vs."""
    rock = np.array([[400, 2500, 4600, 2400], [1500, 2750, 6000, 3300],
                     [4500, 2950, 6900, 3850], [0, 3300, 8000, 4500]], float)  # fmt: skip
    f = np.array([0.004, 0.006])
    k = wavenumber(2 * np.pi * f, DEPTH)
    h, rho, vp = 300.0, 1800.0, 1650.0
    added = {}
    for vs in (300.0, 600.0):
        model = np.vstack([[h, rho, vp, vs], rock])
        added[vs] = dram.forward(DEPTH, f, model) - dram.forward(DEPTH, f, rock)
    uniaxial = k * h / (rho * vp**2)
    assert added[300.0] / uniaxial == pytest.approx(1.0, abs=0.05)
    assert added[300.0] / added[600.0] == pytest.approx(
        1.0, abs=0.05
    )  # Vs doubled: 3 %
    assert np.all(added[300.0] > 0.5 * dram.forward(DEPTH, f, rock))  # and it is large


# -- model ------------------------------------------------------------------------
def test_seafloor_model_layout_and_moho_cut():
    model = auto.SeafloorModel(REFERENCE)
    assert model.names[:3] == [
        "sediment_thickness",
        "sediment_vs_scale",
        "sediment_modulus",
    ]
    assert model.names[-4:] == ["moho_depth", "vp_upper", "vp_lower", "density"]
    assert model.fixed.sum() == pytest.approx(model.moho_bounds[1])
    assert model.reference_vs[0] == pytest.approx(2540.0)
    assert model.reference_vs[-1] == pytest.approx(
        3940.0
    )  # crust continued below the Moho
    theta = model.reference_theta(sediment=0.0)
    theta[model._index["moho_depth"]] = 2000.0
    layers = model.layers(theta)
    assert layers[:-1, 0].sum() == pytest.approx(2000.0)
    assert layers[-1, 3] == theta[model._index["vs_mantle"]]
    theta[model._index["sediment_thickness"]] = 250.0
    layers = model.layers(theta)
    assert layers[:-1, 0].sum() == pytest.approx(2250.0)
    assert len(layers) == model.sediment_sublayers + 5 + 1
    assert model.vs_at(theta, [3000.0])[0] == theta[model._index["vs_mantle"]]
    assert model.valid(theta)
    theta[model._index["vs_crust_5"]] = (
        4400.0  # last layer above the Moho, faster than the mantle
    )
    assert not model.valid(theta)
    theta[model._index["vs_crust_5"]] = 3700.0
    theta[model._index["vs_crust_7"]] = 4700.0  # below the Moho: has no effect
    assert model.valid(theta)


def test_compiled_layers_equal_the_reference():
    """The compiled layer builder and the NumPy one give the same model."""
    rng = np.random.default_rng(4)
    for kwargs in ({}, dict(sediment_max=0.0), dict(upper_crust=300.0)):
        model = auto.SeafloorModel(REFERENCE, **kwargs)
        for _ in range(200):
            theta = model.prior_draw(rng)
            if "sediment_thickness" in model._index and rng.random() < 0.2:
                theta[0] = rng.uniform(0, 1.0)  # no sediment rows below 1 m
            a, b = model.layers(theta), model.layers_reference(theta)
            assert a.shape == b.shape
            assert np.allclose(a, b, rtol=1e-12, atol=1e-9)


def test_prior_draws_follow_the_stated_prior():
    model = auto.SeafloorModel(REFERENCE)
    rng = np.random.default_rng(0)
    draws = np.array([model.prior_draw(rng) for _ in range(3000)])
    assert all(model.valid(t) for t in draws[:200])
    ix = model._index
    i = ix["vs_crust_1"]
    dev = np.log(draws[:, i : i + model.n_crust] / model.reference_vs)
    # shallow layers are far from the bounds and from the Moho rule
    assert np.std(dev[:, 0]) == pytest.approx(model.vs_sigma, rel=0.1)
    assert np.corrcoef(dev[:, 0], dev[:, 1])[0, 1] == pytest.approx(
        model.vs_correlation, abs=0.06
    )
    assert np.std(draws[:, ix["vp_upper"]]) == pytest.approx(model.vp_sigma, rel=0.1)
    assert np.std(draws[:, ix["density"]]) == pytest.approx(
        model.density_sigma, rel=0.1
    )
    lo, hi = model.bounds[:, 0], model.bounds[:, 1]
    assert np.all(draws >= lo) and np.all(draws <= hi)


def test_log_prior_is_the_density_of_prior_draws():
    """Metropolis sampling of log_prior alone reproduces the direct prior draws."""
    model = auto.SeafloorModel(REFERENCE, sediment_max=0.0)
    eta = dram.forward(DEPTH, FREQS, REFERENCE)
    data = dram.ComplianceData.build(
        FREQS, eta, DEPTH, random=0.01 * eta, gain_rel=0, floor_rel=0
    )
    rng = np.random.default_rng(1)
    direct = np.array([model.prior_draw(rng) for _ in range(4000)])
    starts = [np.append(t, 1.0) for t in direct[:4]]
    post = dram.sample(data, model, prior_only=True, starts=starts, n_chains=4, n_iter=40000,
                       n_burn=10000, thin=10, processes=1, gain_bounds=(0.8, 1.2), gain_sigma=0.05,
                       seed=2)  # fmt: skip
    flat = post.theta.reshape(-1, post.theta.shape[-1])
    for name in ("vs_crust_1", "vs_crust_3", "moho_depth", "vp_lower", "density"):
        j = model._index[name]
        a, b = (
            np.quantile(flat[:, j], [0.16, 0.5, 0.84]),
            np.quantile(direct[:, j], [0.16, 0.5, 0.84]),
        )
        scale = b[2] - b[0]
        assert np.all(np.abs(a - b) < 0.15 * scale), name
    assert np.std(flat[:, -1]) == pytest.approx(0.05, rel=0.12)  # Gaussian gain prior
    assert np.mean(flat[:, -1]) == pytest.approx(1.0, abs=0.01)


def test_reference_from_crust1_leaves_out_sediment():
    column = dict(layers=[("upper_sediments", 200.0, 1800.0, 400.0, 1900.0),
                          ("upper_crust", 1170.0, 4700.0, 2540.0, 2400.0),
                          ("lower_crust", 4560.0, 6900.0, 3940.0, 2960.0)],
                  mantle=(7730.0, 4310.0, 3190.0))  # fmt: skip
    ref = auto.reference_from_crust1(column)
    assert ref.shape == (3, 4) and ref[0, 3] == 2540.0 and ref[-1, 3] == 4310.0


# -- data -------------------------------------------------------------------------
def test_compliance_data_caps_the_band_and_omits_the_gain():
    f = np.linspace(0.005, 0.03, 21)
    eta = dram.forward(DEPTH, f, REFERENCE)
    data = auto.compliance_data(f, eta, DEPTH, sigma=0.002 * eta)
    fc = np.sqrt(9.79329 / (2 * np.pi * DEPTH))
    assert data.freqs.max() <= 0.8 * fc + 1e-9 < f.max()
    assert "gain" not in data.components and "floor" not in data.components
    assert set(data.components) == {"random", "tilt", "smooth"}
    with pytest.raises(ValueError):
        auto.compliance_data(f, eta, DEPTH)
    with pytest.raises(ValueError):
        auto.compliance_data(f, eta, DEPTH, sigma=0.002 * eta, band=(0.028, 0.03))


# -- inversion --------------------------------------------------------------------
@pytest.fixture(scope="module")
def synthetic_result():
    model = auto.SeafloorModel(REFERENCE)
    rng = np.random.default_rng(5)
    theta = model.reference_theta(sediment=120.0)
    i = model._index["vs_crust_1"]
    theta[i : i + 4] *= [0.8, 0.8, 0.85, 0.9]  # a slower upper crust than the reference
    eta = dram.forward(DEPTH, FREQS, model.layers(theta))
    cov_rel = default_covariance(FREQS, eta, 0.0015 * eta)
    d = eta + np.linalg.cholesky(cov_rel * np.outer(eta, eta)) @ rng.standard_normal(
        len(eta)
    )
    data = auto.compliance_data(FREQS, d, DEPTH, cov_rel=cov_rel, tilt_rel=0.02)
    kwargs = dict(gain_sigma=0.05, processes=1, n_chains=4, n_iter=12000, n_burn=5000, thin=5,
                  max_rounds=1, n_prior=800, seed=3)  # fmt: skip
    return model, theta, data, auto.invert(data, model, **kwargs), kwargs


def test_invert_recovers_a_synthetic_model(synthetic_result):
    model, theta, data, res, _ = synthetic_result
    s = res.summary()
    assert s["max_rhat"] < 1.1
    assert s["chi2"]["median"] < s["quality"]["chi2_limit"]
    prof = res.profile(np.array([200.0, 500.0, 1000.0, 2000.0]))
    truth = model.vs_at(theta, prof["z"])
    inside = (truth >= prof["vs"][0]) & (truth <= prof["vs"][-1])  # 5-95 %
    assert inside.sum() >= 3
    sed = s["parameters"]["sediment_thickness"]
    assert sed["p05"] <= 120.0 <= sed["p95"]
    assert sed["p95"] - sed["p05"] < 0.5 * 2000.0  # the data constrain the sediment


def test_profile_marks_resolved_depths(synthetic_result):
    model, _, _, res, _ = synthetic_result
    prof = res.profile()
    assert prof["resolved"].any() and not prof["resolved"].all()
    deep = prof["z"] > prof["reach"]
    assert not prof["resolved"][deep].any()
    res = res.parameter_resolution()
    best = min(res[f"vs_crust_{k}"] for k in range(1, model.n_crust + 1))
    assert (
        best < 0.6 < res[f"vs_crust_{model.n_crust}"]
    )  # some layer resolved, the deepest not
    assert res["sediment_thickness"] < 0.3 and res["density"] > 0.8
    assert prof["modulus"].shape == prof["vs"].shape


def test_invert_is_reproducible_and_saves(synthetic_result, tmp_path):
    model, _, data, res, kwargs = synthetic_result
    again = auto.invert(data, model, **kwargs)
    assert np.array_equal(again.posterior.theta, res.posterior.theta)
    res.save(tmp_path / "sta_auto")
    saved = np.load(tmp_path / "sta_auto.npz")
    assert (
        "prior_theta" in saved.files
        and saved["theta"].shape == res.posterior.theta.shape
    )
    table = np.genfromtxt(tmp_path / "sta_auto_profile.csv", delimiter=",", names=True)
    assert "vs_p50" in table.dtype.names and "resolved" in table.dtype.names
    import json

    record = json.loads((tmp_path / "sta_auto.json").read_text())
    assert record["quality"]["label"] in ("good", "usable", "poor")
    assert record["model"]["model"] == "SeafloorModel"
