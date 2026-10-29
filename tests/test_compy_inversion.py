"""Tests of the Bayesian inversion workflow in compy_inversion."""

import numpy as np
import pytest

import compy_inversion as ci
import inv_compy as inv


def random_model(rng):
    n = rng.integers(2, 12)
    vs = rng.uniform(80, 4600, n)
    vp = vs * rng.uniform(1.6, 8, n)
    rho = rng.uniform(1400, 3300, n)
    thick = rng.uniform(5, 3000, n)
    thick[-1] = 0
    vs[-1] = rng.uniform(1500, 4600)
    vp[-1] = 1.8 * vs[-1]
    return np.column_stack([thick, rho, vp, vs])


def test_forward_matches_legacy_propagator():
    rng = np.random.default_rng(1)
    f = np.geomspace(0.002, 0.04, 12)
    compared = 0
    for _ in range(120):
        model, depth = random_model(rng), rng.uniform(800, 5500)
        try:
            legacy = inv.calc_norm_compliance(depth, f, model)
        except ValueError:
            with pytest.raises(ValueError):
                ci.forward(depth, f, model)
            continue
        # Agreement is limited by cancellation inside the propagator when the
        # infragravity phase speed is far below the seismic velocities.
        np.testing.assert_allclose(ci.forward(depth, f, model), legacy, rtol=2e-6)
        compared += 1
    assert compared > 100


def test_forward_half_space_static_limit():
    rho, vp, vs = 2800.0, 6000.0, 3400.0
    mu, lam = rho * vs**2, rho * vp**2 - 2 * rho * vs**2
    static = (lam + 2 * mu) / (2 * mu * (lam + mu))
    model = np.array([[500.0, rho, vp, vs], [0.0, rho, vp, vs]])
    eta = ci.forward(4000.0, np.array([0.005, 0.01, 0.02]), model)
    # The dynamic correction is of order (c / Vs)**2 < 1e-3 here.
    np.testing.assert_allclose(eta, static, rtol=3e-3)
    assert np.all(np.diff(np.abs(eta / static - 1)) < 0)  # shrinks as c falls


def test_forward_rejects_half_space_slower_than_forcing():
    model = np.array([[100.0, 2000.0, 1600.0, 300.0], [0.0, 2000.0, 1500.0, 50.0]])
    with pytest.raises(ValueError):
        ci.forward(4000.0, np.array([0.005]), model)


def test_attraction_term_matches_gravity_correction():
    from compy_processing import gravity_corrected_compliance

    f, depth = np.array([0.006, 0.01, 0.015]), 4000.0
    zero = gravity_corrected_compliance(f, depth, np.zeros(3), np.ones(3))
    np.testing.assert_allclose(ci.attraction_compliance(f, depth), zero, rtol=1e-12)


def test_gain_vector_is_derivative_of_corrected_compliance():
    from compy_processing import gravity_corrected_compliance

    f, depth = np.array([0.006, 0.01, 0.015]), 4000.0
    ratio = np.array([3e-9, 1.5e-9, 1e-9])
    eta = gravity_corrected_compliance(f, depth, ratio, np.ones(3))
    shifted = gravity_corrected_compliance(f, depth, ratio / 1.001, np.ones(3))
    np.testing.assert_allclose((eta - shifted) / 1e-3,
                               eta - ci.attraction_compliance(f, depth), rtol=2e-3)  # fmt: skip


def test_covariance_terms_and_chi2():
    f = np.array([0.008, 0.01, 0.012, 0.014])
    eta = np.array([2e-11, 2.3e-11, 2.6e-11, 3e-11])
    rng = np.random.default_rng(0)
    subsets = eta * (1 + 0.02 * rng.standard_normal((12, 4)))
    data = ci.ComplianceData.build(f, eta, 3000.0, random=0.002 * eta, subsets=subsets,
                                   gain_rel=0.05, floor_rel=0.01)  # fmt: skip
    assert set(data.components) == {"random", "subsets", "gain", "floor"}
    np.testing.assert_allclose(sum(data.components.values()), data.covariance)
    g = eta - ci.attraction_compliance(f, 3000.0)
    np.testing.assert_allclose(data.components["gain"], 0.05**2 * np.outer(g, g))
    # A pure gain error costs almost nothing; the same size independent error does.
    assert data.chi2(eta - 0.05 * g) < 1.5
    r = np.array([1, -1, 1, -1]) * 0.05 * eta
    assert data.chi2(eta + r) > 20
    resid = rng.standard_normal(4) * 1e-12
    expected = resid @ np.linalg.solve(data.covariance, resid)
    np.testing.assert_allclose(data.chi2(eta - resid), expected, rtol=1e-10)


def test_from_npz_restricts_to_band(tmp_path):
    f = np.arange(1, 11) * 1e-3
    path = tmp_path / "c.npz"
    np.savez(path, freqs=f, compliance=np.linspace(1, 2, 10) * 1e-11,
             uncertainty=np.full(10, 1e-13), band=np.array([0.003, 0.007]), depth=3000.0)  # fmt: skip
    data = ci.ComplianceData.from_npz(path)
    np.testing.assert_allclose(data.freqs, [0.003, 0.004, 0.005, 0.006, 0.007])
    assert data.depth == 3000.0


def test_crustal_model_geometry_and_elastic_relations():
    model = ci.CrustalModel(moho_depth=6900.0, sediment_max=1000.0)
    np.testing.assert_allclose(model.fixed, [100, 200, 400, 800, 1600])
    assert model.names[0] == "sediment_thickness" and model.names[-1] == "vs_mantle"
    theta = np.array([250.0, 1800, 2200, 2900, 3500, 3800, 3950, 4300])
    layers = model.layers(theta)
    assert layers[:-1, 0].sum() == pytest.approx(6900.0)
    assert layers[:8, 0].sum() == pytest.approx(250.0)  # sediment sublayers
    assert np.all(layers[:, 2] > layers[:, 3])
    np.testing.assert_allclose(layers[8:, 2], 1.85 * layers[8:, 3])
    np.testing.assert_allclose(layers[8:, 1], 1850 + 0.165 * layers[8:, 2])
    np.testing.assert_allclose(model.vs_at(theta, [0, 260, 349, 351, 7000]),
                               [layers[0, 3], 1800, 1800, 2200, 4300])  # fmt: skip
    assert model.roughness(theta) == pytest.approx(np.sum(np.diff(theta[1:7], 2) ** 2))
    no_sed = model.layers(np.r_[0.0, theta[1:]])
    assert no_sed[0, 3] == 1800 and no_sed[:-1, 0].sum() == pytest.approx(6900.0)


def test_crustal_model_optional_parameters():
    model = ci.CrustalModel(moho_depth=7000.0, sediment_scale=(0.5, 2.0),
                            moho_shift=(-1000.0, 1000.0), vp_vs=None)  # fmt: skip
    assert model.names[:2] == ["sediment_thickness", "sediment_vs_scale"]
    assert model.names[-1] == "moho_shift"
    theta = np.array([300.0, 2.0] + [2500.0] * model.n_crust + [4400.0, -800.0])
    layers = model.layers(theta)
    assert layers[:-1, 0].sum() == pytest.approx(6200.0)
    np.testing.assert_allclose(layers[:8, 3], 2 * ci.hamilton_vs(np.linspace(18.75, 281.25, 8)))
    assert np.all(layers[:, 2] > 1.4 * layers[:, 3])
    with pytest.raises(ValueError):
        ci.CrustalModel(moho_depth=1200.0, sediment_max=1000.0)


def test_sediment_profiles_are_continuous():
    for z in (36.0, 120.0):
        assert ci.hamilton_vs(z - 1e-9) == pytest.approx(ci.hamilton_vs(z), abs=1.1)
    assert ci.hamilton_vp(1500.0) == pytest.approx(2667.0, abs=1.0)


def test_crust1_column(tmp_path):
    lines = {key: ["0 0 0 0 0 0 0 0 0"] * 64800 for key in ("bnds", "vp", "vs", "rho")}
    row = 112 * 360 + 233  # cell centred on -22.5, 53.5
    lines["bnds"][row] = "0.00 -4.42 -4.42 -4.68 -4.68 -4.68 -5.34 -6.80 -11.31"
    lines["vp"][row] = "1.50 3.81 1.75 0.00 0.00 5.00 6.50 7.10 8.12"
    lines["vs"][row] = "0.00 1.94 0.34 0.00 0.00 2.70 3.70 4.05 4.51"
    lines["rho"][row] = "1.02 0.92 1.82 0.00 0.00 2.55 2.85 3.05 3.34"
    for key, rows in lines.items():
        (tmp_path / f"crust1.{key}").write_text("\n".join(rows) + "\n")
    col = ci.crust1_column(-22.7152, 53.1595, tmp_path)
    assert col["cell"] == (-22.5, 53.5)
    assert col["water_depth"] == pytest.approx(4420.0)
    assert col["sediment_thickness"] == pytest.approx(260.0)
    assert col["moho_depth"] == pytest.approx(6890.0)
    model = ci.crust1_model(col)
    np.testing.assert_allclose(model[:, 3], [340, 2700, 3700, 4050, 4510])
    np.testing.assert_allclose(model[:-1, 0], [260, 660, 1460, 4510])


class _Gaussian:
    """Correlated 2-D Gaussian restricted to a box, as a dram_chain target."""

    def __init__(self, bounds):
        self.bounds = bounds
        self.precision = np.linalg.inv(np.array([[1.0, 0.8], [0.8, 1.0]]))

    def __call__(self, x):
        if np.any(x < self.bounds[:, 0]) or np.any(x > self.bounds[:, 1]):
            return None
        return -0.5 * x @ self.precision @ x, 0.0, x


@pytest.mark.parametrize("delayed_rejection", [True, False])
def test_dram_samples_truncated_gaussian(delayed_rejection):
    bounds = np.array([[-1.0, 3.0], [-3.0, 3.0]])
    target = _Gaussian(bounds)
    out = ci.dram_chain(target, np.array([2.5, -2.5]), bounds, n_iter=120000, n_burn=10000,
                        seed=3, delayed_rejection=delayed_rejection)  # fmt: skip
    g = np.stack(np.meshgrid(np.linspace(-1, 3, 401), np.linspace(-3, 3, 601), indexing="ij"), -1)
    w = np.exp(-0.5 * np.einsum("...i,ij,...j", g, target.precision, g))
    w /= w.sum()
    mean = np.einsum("ij,ijk->k", w, g)
    var = np.einsum("ij,ijk->k", w, (g - mean) ** 2)
    np.testing.assert_allclose(out["theta"].mean(0), mean, atol=0.04)
    np.testing.assert_allclose(out["theta"].var(0), var, rtol=0.06)
    assert 0.1 < out["acceptance"] < 0.5


def test_hierarchical_smoothing_weight_marginal():
    # One parameter in [-1, 1] with roughness theta**2 (k = 1) and alpha on
    # [0.5, 50]: the joint alpha**(k/2 - 1) exp(-alpha theta**2 / 2) has the
    # theta marginal integral alpha**(-1/2) exp(-alpha theta**2 / 2) dalpha.
    from scipy.integrate import quad

    def target(x):
        return (0.0, float(x[0] ** 2), x) if abs(x[0]) <= 1 else None

    out = ci.dram_chain(target, np.array([0.9]), np.array([[-1.0, 1.0]]), n_iter=150000,
                        n_burn=10000, seed=5, n_rough=1, alpha_bounds=(0.5, 50.0))  # fmt: skip
    t = np.linspace(-1, 1, 801)
    dens = np.array([quad(lambda a: a**-0.5 * np.exp(-a * x * x / 2), 0.5, 50)[0] for x in t])
    dens /= np.trapezoid(dens, t)
    expected = np.trapezoid(np.abs(t) * dens, t)
    assert np.mean(np.abs(out["theta"][:, 0])) == pytest.approx(expected, rel=0.04)
    assert np.all((out["alpha"] >= 0.5) & (out["alpha"] <= 50))


def test_truncated_gamma_respects_bounds_and_shape():
    rng = np.random.default_rng(0)
    x = np.array([ci.truncated_gamma(rng, 2.0, 1.0, 0.5, 3.0) for _ in range(20000)])
    assert x.min() >= 0.5 and x.max() <= 3.0
    from scipy.stats import gamma

    d = gamma(2.0)
    m = (d.expect(lambda v: v, lb=0.5, ub=3.0)) / (d.cdf(3.0) - d.cdf(0.5))
    assert x.mean() == pytest.approx(m, rel=0.02)
    far = [ci.truncated_gamma(rng, 2.0, 1e6, 1.0, 2.0) for _ in range(10)]
    assert all(1.0 <= v <= 2.0 for v in far)
    flat = np.array([ci.truncated_gamma(rng, 1.0, 0.0, 1.0, 3.0) for _ in range(20000)])
    assert flat.mean() == pytest.approx(2.0, rel=0.02)


def test_synthetic_recovery_and_prior_only():
    depth = 2918.0
    f = np.arange(6, 15) * 1.220703125e-3
    model = ci.CrustalModel(moho_depth=6900.0, sediment_max=1000.0)
    truth = np.array([150.0, 1800, 2200, 2900, 3500, 3800, 3950, 4300])
    clean = ci.forward(depth, f, model.layers(truth))
    data = ci.ComplianceData.build(f, clean, depth, random=0.003 * clean, gain_rel=0.05,
                                   floor_rel=0.01, station="synthetic")  # fmt: skip
    post = ci.sample(data, model, n_chains=4, n_iter=30000, n_burn=10000, thin=5, seed=2,
                     processes=1)  # fmt: skip
    s = post.summary()
    assert s["max_rhat"] < 1.1
    sed = s["parameters"]["sediment_thickness"]
    assert sed["p05"] < 150 < sed["p95"] and sed["p95"] - sed["p05"] < 200
    for name, value in zip(model.names, truth):
        assert s["parameters"][name]["p05"] < value < s["parameters"][name]["p95"]
    assert s["chi2"]["median"] < 3 * len(f)
    prior = ci.sample(data, model, n_chains=2, n_iter=6000, n_burn=2000, thin=5, seed=3,
                      prior_only=True, processes=1)  # fmt: skip
    ratio = ci.resolution_ratio(post, prior, [75.0, 1000.0])
    assert ratio[0] < 0.6  # the sediment column is resolved
    q = post.vs_quantiles([100.0, 5000.0])
    assert q.shape == (5, 2) and np.all(np.diff(q, axis=0) >= 0)


def test_posterior_save_roundtrip(tmp_path):
    f = np.arange(6, 12) * 1.220703125e-3
    model = ci.CrustalModel(moho_depth=6000.0, sediment_max=0.0)
    truth = np.array([2000.0, 2500, 3000, 3500, 3800, 3900, 4300])[: len(model.names)]
    data = ci.ComplianceData.build(f, ci.forward(3500.0, f, model.layers(truth)), 3500.0,
                                   floor_rel=0.02)  # fmt: skip
    post = ci.sample(data, model, n_chains=2, n_iter=3000, n_burn=1000, thin=10, processes=1)
    post.save(tmp_path / "run")
    stored = np.load(tmp_path / "run.npz")
    assert stored["theta"].shape == post.theta.shape
    import json

    summary = json.loads((tmp_path / "run.json").read_text())
    assert summary["model"]["names"] == model.names
