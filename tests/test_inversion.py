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
