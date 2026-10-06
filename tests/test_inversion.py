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
