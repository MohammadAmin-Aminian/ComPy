import numpy as np
import pytest
import inv_compy as inv


@pytest.fixture
def problem():
    model = np.array(
        [
            [100.0, 2200.0, 3000.0, 1500.0],
            [1000.0, 2800.0, 6000.0, 3500.0],
            [0.0, 3300.0, 8000.0, 4500.0],
        ]
    )
    f = np.array([0.005, 0.007, 0.01, 0.015])
    return model, f, inv.calc_norm_compliance(4000, f, model)


def test_initial_state_and_every_saved_prediction_match_model(problem):
    model, f, data = problem
    chain, profiles, prior, misfits, predictions, likelihood, rate = (
        inv.invert_compliance_beta(
            data,
            f,
            4000,
            starting_model=model,
            s=np.full(4, 1e-12),
            iteration=25,
            sigma_v=50,
            sigma_h=50,
            alpha=0,
            seed=4,
        )
    )
    np.testing.assert_equal(chain[:, :, 0], model)
    np.testing.assert_allclose(predictions[0], data)
    assert misfits[0, 0] == 0
    assert likelihood[0, 0] == 1
    for i in range(25):
        np.testing.assert_allclose(
            predictions[i], inv.calc_norm_compliance(4000, f, chain[:, :, i])
        )
        np.testing.assert_equal(chain[-1, :, i], model[-1])
        assert np.isclose(chain[:-1, 0, i].sum(), model[:-1, 0].sum())
        assert np.all(chain[:-1, 0, i] > 0)
    assert profiles.shape == (25, 1100, 1)
    assert prior.shape == (1, 1100, 1)
    assert 0 <= rate <= 1


def test_reproducible_sampler_does_not_reset_global_rng(problem):
    model, f, data = problem
    kwargs = dict(starting_model=model, s=1e-12, iteration=6, return_profiles=False)
    np.random.seed(8)
    expected = np.random.random(2)
    np.random.seed(8)
    first = inv.invert_compliance(data, f, 4000, **kwargs)
    np.testing.assert_equal(np.random.random(2), expected)
    second = inv.invert_compliance(data, f, 4000, **kwargs)
    np.testing.assert_equal(first[0], second[0])
    assert first[1] is None


def test_required_uncertainties_and_invalid_frequency(problem):
    model, f, data = problem
    with pytest.raises(ValueError, match="s is required"):
        inv.invert_compliance(data, f, 4000, starting_model=model, iteration=2)
    with pytest.raises(ValueError, match="frequencies"):
        inv.calc_norm_compliance(4000, [0], model)
    with pytest.raises(ValueError):
        inv.Model_V2(1, n_layer=4)


def test_forward_layer_split_invariance(problem):
    model, f, _ = problem
    split = np.vstack([model[0], model[0], model[1:]])
    split[:2, 0] = model[0, 0] / 2
    np.testing.assert_allclose(
        inv.calc_norm_compliance(4000, f, model),
        inv.calc_norm_compliance(4000, f, split),
        rtol=1e-8,
    )


def test_homogeneous_halfspace_static_elastic_limit():
    # Analytic normalized compliance magnitude = (1-nu)/mu.
    rho, vs, nu = 3000.0, 3000.0, 0.25
    vp = vs * np.sqrt((1 - nu) / (0.5 - nu))
    model = np.array([[100.0, rho, vp, vs], [0.0, rho, vp, vs]])
    result = inv.calc_norm_compliance(0.1, np.array([1e-4]), model)
    np.testing.assert_allclose(result, (1 - nu) / (rho * vs**2), rtol=2e-7)


def test_scalar_and_equivalent_per_layer_steps_give_identical_chains(problem):
    model, f, data = problem
    kwargs = dict(starting_model=model, s=0.03 * data, iteration=200, alpha=0, seed=1)
    scalar = inv.invert_compliance_beta(data, f, 4000, sigma_v=40, sigma_h=9, **kwargs)
    layered = inv.invert_compliance_beta(
        data, f, 4000, sigma_v=[40, 40], sigma_h=[9], **kwargs
    )
    np.testing.assert_array_equal(scalar[0], layered[0])
    with pytest.raises(ValueError, match="sigma_v must be a scalar or have 2"):
        inv.invert_compliance_beta(data, f, 4000, sigma_v=[40, 40, 40], **kwargs)
    with pytest.raises(ValueError, match="sigma_h must be finite and positive"):
        inv.invert_compliance_beta(data, f, 4000, sigma_h=[-1.0], **kwargs)


def test_adaptation_tunes_steps_then_freezes(problem):
    model, f, data = problem
    kwargs = dict(starting_model=model, s=0.03 * data, alpha=0, seed=2, sigma_v=5)
    short, long = {}, {}
    inv.invert_compliance_beta(
        data, f, 4000, iteration=600, adapt_until=400, diagnostics=short, **kwargs
    )
    inv.invert_compliance_beta(
        data, f, 4000, iteration=900, adapt_until=400, diagnostics=long, **kwargs
    )
    # Same seed, same adaptation window: frozen scales do not depend on run length.
    np.testing.assert_array_equal(short["sigma_v"], long["sigma_v"])
    np.testing.assert_array_equal(short["sigma_h"], long["sigma_h"])
    # A 5 m/s step is far too small for this weakly constrained layer.
    assert short["sigma_v"][0] > 5
    assert np.all((short["acceptance_vs"] >= 0) & (short["acceptance_vs"] <= 1))


def test_initial_state_starts_chain_but_keeps_prior_from_starting_model(problem):
    model, f, data = problem
    state = model.copy()
    state[0, 3], state[0, 0], state[1, 0] = 900.0, 150.0, 950.0
    chain, _, prior, *_ = inv.invert_compliance_beta(
        data, f, 4000, starting_model=model, s=0.03 * data, iteration=5,
        alpha=0, seed=0, initial_state=state,
    )
    np.testing.assert_array_equal(chain[:, :, 0], state)
    assert prior[0, 0, 0] == model[0, 3]
    outside = state.copy()
    outside[0, 3] = 1900.0  # above the beta prior bound 1.2 x 1500, below Vp
    with pytest.raises(ValueError, match="initial_state"):
        inv.invert_compliance_beta(
            data, f, 4000, starting_model=model, s=0.03 * data, iteration=5,
            initial_state=outside,
        )
    thicker = state.copy()
    thicker[1, 0] += 10  # changes the conserved total finite thickness
    with pytest.raises(ValueError, match="initial_state"):
        inv.invert_compliance_beta(
            data, f, 4000, starting_model=model, s=0.03 * data, iteration=5,
            initial_state=thicker,
        )


def _reference_compliance(H, freq, model, g=9.79329):
    """Independent P-SV solver (Aki & Richards eq. 7.28, matrix exponentials)."""
    from scipy.linalg import expm, qr
    from scipy.optimize import brentq

    def system(k, w, rho, vp, vs):
        mu = rho * vs * vs
        lam = rho * vp * vp - 2 * mu
        l2m = lam + 2 * mu
        return np.array([
            [0.0, k, 1 / mu, 0.0],
            [-k * lam / l2m, 0.0, 0.0, 1 / l2m],
            [k * k * 4 * mu * (lam + mu) / l2m - w * w * rho, 0.0, 0.0, k * lam / l2m],
            [0.0, -w * w * rho, -k, 0.0],
        ])

    out = []
    for w in 2 * np.pi * np.asarray(freq, float):
        k = brentq(lambda x: g * x * np.tanh(x * H) - w * w, 1e-12, 10 * w * w / g + 10, xtol=1e-16)
        val, vec = np.linalg.eig(system(k, w, *model[-1, 1:]))
        R = np.real(vec[:, np.argsort(np.real(val))[:2]])  # decaying with depth
        for row in model[-2::-1]:
            R, _ = qr(expm(-system(k, w, *row[1:]) * row[0]) @ R, mode="economic")
        r = R @ np.array([R[2, 1], -R[2, 0]])  # zero shear traction at seafloor
        out.append(k * abs(r[1] / r[3]))
    return np.array(out)


@pytest.mark.parametrize("H", [1000.0, 4200.0])
def test_forward_model_matches_independent_propagator(H):
    crust_lvz = np.array([
        [300, 1900, 2000, 800], [900, 2500, 4200, 2100], [1300, 2800, 6200, 3500],
        [1500, 2800, 6000, 2700], [2500, 3000, 7000, 3900], [0, 3300, 8000, 4400],
    ], float)
    soft = np.array([[200, 1600, 1600, 150], [800, 2300, 3500, 1500],
                     [5000, 2900, 6500, 3700], [0, 3300, 8000, 4500]], float)
    f = np.geomspace(0.002, 1.5 * np.sqrt(9.81 / (2 * np.pi * H)), 9)
    for model in (crust_lvz, soft):
        np.testing.assert_allclose(
            inv.calc_norm_compliance(H, f, model), _reference_compliance(H, f, model), rtol=2e-5
        )


def test_resampled_model_has_identical_compliance(problem):
    model, f, data = problem
    fine, mid = inv.resample_model(model, 75, max_depth=3000)
    assert np.isclose(fine[:-1, 0].sum(), 3000) and len(mid) == len(fine) - 1
    np.testing.assert_allclose(inv.calc_norm_compliance(4000, f, fine), data, rtol=1e-8)


def test_sensitivity_kernels_integrate_to_uniform_perturbation(problem):
    model, f, data = problem
    z, K = inv.compliance_kernels(4000, f, model, dz=100, max_depth=12000)
    assert K.shape == (len(f), len(z))
    # A uniformly stiffer crust lowers compliance. Individual slices can have
    # small positive kernels: with Vp fixed, raising Vs also lowers lambda
    # (confirmed with the independent propagator above).
    assert np.all(K.sum(axis=1) < 0)
    fine, _ = inv.resample_model(model, 100, 12000)
    scaled = fine.copy()
    scaled[:-1, 3] *= 1.01
    direct = np.log(inv.calc_norm_compliance(4000, f, scaled) / data) / 0.01
    np.testing.assert_allclose((K * 0.1).sum(axis=1), direct, rtol=0.03)
    # Longer water waves sense deeper: below 6 km the lowest frequency has the
    # largest sensitivity and the highest frequency the smallest.
    deep = np.abs(K[:, z > 6000]).sum(axis=1)
    assert deep.argmax() == 0 and deep.argmin() == len(f) - 1


def test_occam_reaches_target_smoothly_from_different_starts():
    truth = np.array([[400, 2000, 2500, 900], [1500, 2600, 5000, 2600],
                      [3000, 2900, 6800, 3700], [0, 3300, 8000, 4400]], float)
    f = np.geomspace(0.004, 0.018, 10)
    clean = inv.calc_norm_compliance(4200, f, truth)
    s = 0.02 * clean
    data = clean + np.random.default_rng(5).normal(0, s)
    thick = 200 * 1.3 ** np.arange(9)
    mid = np.cumsum(thick) - thick / 2
    idx = np.searchsorted(np.cumsum(truth[:-1, 0]), mid, side="right")
    grid = np.vstack([np.column_stack([thick, truth[idx, 1], truth[idx, 2],
                                       np.full(9, 2000.0)]), truth[-1]])
    results = []
    for start in (1500.0, 3000.0):
        grid[:-1, 3] = np.minimum(start, 0.9 * grid[:-1, 2])
        res = inv.invert_occam(data, f, 4200, grid, s)
        assert res["chi2"] <= len(f) + 1e-9
        assert all(np.diff([h["roughness"] for h in res["history"]
                            if h["chi2"] <= len(f)]) <= 0)
        results.append(res["model"][:-1, 3])
    np.testing.assert_allclose(results[0], results[1], rtol=0.08)


def test_explicit_bounds_and_fixed_geometry(problem):
    model, f, data = problem
    kwargs = dict(starting_model=model, s=0.03 * data, iteration=400, alpha=0,
                  seed=3, return_profiles=False, sigma_v=300, adapt_until=200)
    chain = inv.invert_compliance_beta(
        data, f, 4000, vs_bounds=(1000, [2900, 4000]), vary_thickness=False, **kwargs
    )[0]
    assert np.all(chain[:, 0, :] == model[:, 0, None])  # geometry fixed
    assert chain[0, 3].min() >= 1000 and chain[0, 3].max() <= 2900
    assert chain[1, 3].max() <= 4000 and np.all(chain[-1, 3] == model[-1, 3])
    # The default beta prior would cap layer 2 at 1.1 x 3500; explicit bounds do not.
    assert chain[1, 3].max() > 3500 * 1.1 or chain[1, 3].min() < 3500 * 0.9
    for bad in [(3000, 1000), (0, 2000), ([1000], 3000), 5]:
        with pytest.raises(ValueError, match="vs_bounds"):
            inv.invert_compliance_beta(data, f, 4000, vs_bounds=bad, **kwargs)
