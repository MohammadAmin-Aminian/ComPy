"""Forward model, published parameterizations, minimum-structure inversion and the Moho rule."""

import numpy as np
import pytest

from compy.inversion import dram, impedance, literature_mcmc, minimum_structure, occam
from compy.inversion.propagator import normalized_compliance

DEPTH = 4000.0
FREQS = np.array([0.005, 0.0073, 0.01, 0.015, 0.02])
MODELS = {
    "sediment over crust": [
        [300.0, 1900.0, 1700.0, 400.0],
        [1000.0, 2600.0, 5000.0, 2700.0],
        [4000.0, 2900.0, 6800.0, 3800.0],
        [0.0, 3300.0, 8000.0, 4500.0],
    ],
    "layer slower than the water wave": [
        [30.0, 1500.0, 1500.0, 60.0],
        [200.0, 2000.0, 1800.0, 300.0],
        [2000.0, 2800.0, 6000.0, 3400.0],
        [0.0, 3300.0, 8000.0, 4500.0],
    ],
    "low-velocity zone": [
        [500.0, 2500.0, 4500.0, 2500.0],
        [800.0, 2300.0, 3500.0, 900.0],
        [3000.0, 2900.0, 6800.0, 3800.0],
        [0.0, 3300.0, 8000.0, 4500.0],
    ],
}


def zha_webb():
    return literature_mcmc.ZhaWebb2016Model(
        7000.0, [0.0, 1200.0, 2400.0], [4700.0, 6300.0, 6900.0], 8000.0
    )


def test_wavenumber_solves_the_dispersion_relation():
    omega = 2 * np.pi * FREQS
    k = impedance.gravity_wave_wavenumber(omega, DEPTH, 9.8)
    np.testing.assert_allclose(9.8 * k * np.tanh(k * DEPTH), omega**2, rtol=1e-13)


@pytest.mark.parametrize(
    "vp, vs, rho",
    [(6000.0, 3400.0, 2800.0), (3000.0, 1200.0, 2200.0), (1600.0, 400.0, 1900.0)],
)
def test_homogeneous_half_space_matches_the_exact_formula(vp, vs, rho):
    forward = impedance.ImpedanceForward(DEPTH, FREQS)
    got = forward(np.array([0.0]), np.array([rho]), np.array([vp]), np.array([vs]))
    exact = impedance.exact_halfspace_compliance(forward.k, forward.omega, rho, vp, vs)
    np.testing.assert_allclose(got, exact, rtol=1e-10)


@pytest.mark.parametrize("name", sorted(MODELS))
def test_layered_models_match_the_matrix_exponential(name):
    model = np.array(MODELS[name])
    forward = impedance.ImpedanceForward(DEPTH, FREQS)
    got = impedance.compliance_by_impedance(DEPTH, FREQS, model)
    reference = [
        impedance.matrix_exponential_compliance(k, w, *model.T)
        for k, w in zip(forward.k, forward.omega)
    ]
    np.testing.assert_allclose(got, reference, rtol=1e-7)


@pytest.mark.parametrize("name", sorted(MODELS))
def test_both_forward_models_agree(name):
    model = np.array(MODELS[name])
    np.testing.assert_allclose(
        impedance.compliance_by_impedance(DEPTH, FREQS, model),
        normalized_compliance(DEPTH, FREQS, model),
        rtol=2e-5,
    )


def test_forward_rejects_invalid_models():
    with pytest.raises(ValueError):
        impedance.compliance_by_impedance(
            DEPTH, FREQS, np.array([[100.0, 2000.0, 3000.0]])
        )
    with pytest.raises(ValueError):
        impedance.compliance_by_impedance(
            DEPTH,
            FREQS,
            np.array([[100.0, -1.0, 3000.0, 1000.0], [0.0, 3000.0, 8000.0, 4500.0]]),
        )


def test_zha_webb_layers_follow_the_published_setup():
    model = zha_webb()
    assert model.fixed.tolist() == [100.0, 200.0, 400.0, 800.0, 1600.0]
    theta = np.r_[120.0, 1600.0, 2000.0, 2500.0, 3100.0, 3700.0, 4100.0, 4500.0]
    thickness, rho, vp, vs = model.layers(theta)
    assert thickness[:3].tolist() == [
        50.0,
        50.0,
        20.0,
    ]  # sublayers, the last one partial
    assert thickness[:-1].sum() == pytest.approx(7000.0)
    assert np.all(rho[:3] == 2000.0)
    np.testing.assert_allclose(vp[3:] / vs[3:], 1.85)
    # Density is fixed by the starting Vp, not by the sampled Vs.
    assert rho[3] == pytest.approx(1850 + 0.165 * 4700.0)
    assert model.layers(np.r_[120.0, 900.0, theta[2:]])[1][3] == rho[3]


def test_zero_sediment_has_no_sediment_rows():
    model = zha_webb()
    theta = np.r_[0.0, 1600.0, 2000.0, 2500.0, 3100.0, 3700.0, 4100.0, 4500.0]
    assert len(model.layers(theta)[0]) == model.n_crust + 1


def test_zheng_layers_follow_the_published_setup():
    model = literature_mcmc.Zheng2024Model(7000.0)
    theta = np.r_[250.0, 300.0, 900.0, 1500.0, 2500.0, 3500.0, 4200.0, 4600.0]
    assert model.is_valid(theta)
    thickness, rho, vp, vs = model.layers(theta)
    assert thickness[:-1].sum() == pytest.approx(7000.0)
    assert (
        np.all(rho[:8] == 2000.0) and np.all(rho[8:-1] == 2800.0) and rho[-1] == 3000.0
    )
    assert np.all(
        np.diff(vs) >= 0
    )  # monotonic, with a jump at the sediment base and the Moho
    assert vs[-2] == 4200.0 and vs[-1] == 4600.0
    assert not model.is_valid(
        np.array([250.0, 900.0, 300.0, *theta[3:]])
    )  # sediment Vs decreasing
    assert not model.is_valid(
        np.array([*theta[:3], 500.0, *theta[4:]])
    )  # negative jump at the sediment base


@pytest.mark.parametrize(
    "model, theta",
    [
        (
            zha_webb(),
            np.r_[200.0, 1600.0, 2000.0, 2500.0, 3100.0, 3700.0, 4300.0, 4200.0],
        ),
        (
            literature_mcmc.Zheng2024Model(7000.0),
            np.r_[250.0, 300.0, 900.0, 1500.0, 2500.0, 3500.0, 4300.0, 4200.0],
        ),
    ],
)
def test_published_models_reject_a_velocity_drop_at_the_moho(model, theta):
    assert not model.is_valid(theta)
    fixed = theta.copy()
    fixed[-1] = 4400.0
    assert model.is_valid(fixed)


def test_misfit_filters_follow_the_papers():
    chi2 = np.array([10.0, 15.0, 19.9, 20.1, 40.0])
    np.testing.assert_array_equal(
        zha_webb().retained(chi2, 10), [True, True, True, False, False]
    )
    # Minimum rms 1: below 2. Minimum rms 0.1: below 0.6.
    zheng = literature_mcmc.Zheng2024Model(7000.0)
    np.testing.assert_array_equal(
        zheng.retained(np.array([10.0, 39.0, 41.0]), 10), [True, True, False]
    )
    np.testing.assert_array_equal(
        zheng.retained(np.array([0.1, 3.5, 3.7]), 10), [True, True, False]
    )


def test_monthly_error_needs_three_curves():
    with pytest.raises(ValueError):
        literature_mcmc.monthly_standard_deviation([[1.0, 2.0], [1.1, 2.1]])
    np.testing.assert_allclose(
        literature_mcmc.monthly_standard_deviation([[1.0], [2.0], [3.0]]), [1.0]
    )


@pytest.mark.parametrize(
    "make", [zha_webb, lambda: literature_mcmc.Zheng2024Model(7000.0)]
)
def test_sampler_is_repeatable_and_keeps_the_moho_rule(make):
    model = make()
    freqs = FREQS[1:4]
    theta = model.initial(np.random.default_rng(0))
    while not model.is_valid(theta):
        theta = model.initial(np.random.default_rng(1))
    eta = impedance.ImpedanceForward(DEPTH, freqs)(*model.layers(theta))
    kwargs = dict(n_chains=2, n_iter=1200, n_burn=600, thin=10, seed=3, processes=1)
    a = literature_mcmc.sample_literature_model(
        model, DEPTH, freqs, eta, 0.05 * eta, **kwargs
    )
    b = literature_mcmc.sample_literature_model(
        model, DEPTH, freqs, eta, 0.05 * eta, **kwargs
    )
    np.testing.assert_array_equal(a.theta, b.theta)
    assert a.theta.shape == (2, 60, len(model.names))
    assert np.all(a.theta[:, :, -1] > a.theta[:, :, model.index_last_crust])
    assert a.summary()["moho_jump"]["minimum"] > 0


def test_sampler_checks_its_input():
    with pytest.raises(ValueError):
        literature_mcmc.sample_literature_model(
            zha_webb(), DEPTH, FREQS, FREQS[:2], FREQS[:2]
        )
    with pytest.raises(ValueError):
        literature_mcmc.sample_literature_model(
            zha_webb(), DEPTH, FREQS, np.ones(5), np.ones(5), n_iter=10, n_burn=10
        )


def test_thin_layers_follow_crawford_2004():
    reference = np.array(
        [
            [1200.0, 2400.0, 4700.0, 2500.0],
            [5800.0, 2960.0, 6900.0, 3900.0],
            [0.0, 3200.0, 8000.0, 4400.0],
        ]
    )
    model, edges = minimum_structure.crawford_thin_layers(reference)
    assert model[0, 0] == 50.0 and model[1, 0] == pytest.approx(55.0)
    assert edges[-1] == 7000.0 and model[:-1, 0].sum() == pytest.approx(7000.0)
    np.testing.assert_allclose(model[:-1, 3], model[:-1, 2] / 1.87)
    assert model[0, 2] == 4700.0 and model[-2, 2] == 6900.0
    np.testing.assert_array_equal(model[-1], reference[-1])


def thin_problem():
    truth = np.array(
        [
            [400.0, 2000.0, 2500.0, 900.0],
            [1500.0, 2600.0, 5000.0, 2600.0],
            [3000.0, 2900.0, 6800.0, 3700.0],
            [0.0, 3300.0, 8000.0, 4400.0],
        ]
    )
    freqs = np.geomspace(0.004, 0.018, 10)
    clean = impedance.compliance_by_impedance(4200.0, freqs, truth)
    sigma = 0.02 * clean
    data = clean + np.random.default_rng(5).normal(0, sigma)
    grid, _ = minimum_structure.crawford_thin_layers(truth)
    return data, freqs, grid, sigma


def test_minimum_structure_reaches_the_target_from_two_starts():
    data, freqs, grid, sigma = thin_problem()
    a = minimum_structure.invert_minimum_structure(data, freqs, 4200.0, grid, sigma)
    b = minimum_structure.invert_minimum_structure(
        data, freqs, 4200.0, grid, sigma, start=np.full(len(grid) - 1, 2000.0)
    )
    assert a.reached_target and b.reached_target
    np.testing.assert_allclose(a.model[:-1, 3], b.model[:-1, 3], rtol=0.03)
    fits = [h["roughness"] for h in a.history if h["chi2"] <= len(freqs)]
    assert all(np.diff(fits) <= 0)


def test_minimum_structure_keeps_the_moho_rule():
    data, freqs, grid, sigma = thin_problem()
    grid = grid.copy()
    grid[-1, 3] = 3000.0  # a half-space slow enough to bind
    held = minimum_structure.invert_minimum_structure(data, freqs, 4200.0, grid, sigma)
    assert held.model[:-1, 3].max() < held.model[-1, 3]
    free = minimum_structure.invert_minimum_structure(
        data, freqs, 4200.0, grid, sigma, moho_jump=False
    )
    assert free.model[-2, 3] > free.model[-1, 3]


def test_invert_occam_keeps_the_moho_rule_only_on_request():
    data, freqs, grid, sigma = thin_problem()
    grid = grid.copy()
    grid[-1, 3] = 3000.0
    held = occam.invert_occam(data, freqs, 4200.0, grid, sigma, moho_jump=True)
    assert held["model"][:-1, 3].max() <= grid[-1, 3]
    free = occam.invert_occam(data, freqs, 4200.0, grid, sigma)
    assert free["model"][-2, 3] > grid[-1, 3]


def test_crustal_model_applies_the_moho_rule_only_when_asked():
    freqs = FREQS[1:4]
    theta_ok = None
    for flag in (False, True):
        model = dram.CrustalModel(moho_depth=7000.0, require_moho_jump=flag)
        theta = np.mean(model.bounds, axis=1)
        theta[model._index["vs_mantle"]] = 4000.0
        theta[model._index["vs_crust_1"] + model.n_crust - 1] = 4200.0
        eta = dram.forward(DEPTH, freqs, model.layers(theta))
        data = dram.ComplianceData.build(
            freqs, eta, DEPTH, random=0.05 * eta, gain_rel=0.0
        )
        state = dram._Target(data, model, False)(theta)
        assert (state is None) == flag
        theta_ok = theta
    theta_ok[model._index["vs_mantle"]] = 4400.0
    assert dram._Target(data, model, False)(theta_ok) is not None


def test_crust1_start_obeys_the_moho_rule():
    column = dict(
        layers=[
            ("upper_crust", 1000.0, 5000.0, 2700.0, 2500.0),
            ("lower_crust", 6000.0, 7100.0, 4400.0, 3000.0),
        ],
        mantle=(7800.0, 4100.0, 3200.0),
        moho_depth=7000.0,
        sediment_thickness=0.0,
    )
    model = dram.CrustalModel(
        moho_depth=7000.0, require_moho_jump=True, vp_vs_range=(1.7, 2.4)
    )
    theta = dram.crust1_start(model, column)
    last = theta[model._index["vs_crust_1"] + model.n_crust - 1]
    assert theta[model._index["vs_mantle"]] > last
    assert model.in_bounds(theta)
