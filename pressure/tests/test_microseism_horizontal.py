import numpy as np
import pytest

from compy_pressure.microseism_horizontal import conditional_products, corrected_fourier
from compy_pressure.paired_spectra import spectral_terms


def test_fourier_normalization_reproduces_existing_pair_and_upward_sign():
    rng = np.random.default_rng(12)
    samples = rng.normal(size=(4, 7200))
    rate = 2
    f = np.fft.rfftfreq(7200, 1 / rate)
    mask = (f >= 0.02) & (f <= 0.1)
    full = np.array([np.full_like(f, gain, dtype=complex) for gain in (4, 3, 5, 8)])
    full *= np.exp(1j * f)
    x = corrected_fourier(samples, full[:, mask], rate, mask, upward_sign=-1)
    pp, zz, zp, _ = spectral_terms(
        samples[0], samples[3], rate, full[0], full[3], sign=-1
    )
    np.testing.assert_allclose(abs(x[0]) ** 2, pp[mask], rtol=1e-12)
    np.testing.assert_allclose(abs(x[3]) ** 2, zz[mask], rtol=1e-12)
    np.testing.assert_allclose(x[3] * x[0].conj(), zp[mask], rtol=1e-12)


def random_complex(rng, shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def test_horizontal_independence_preserves_pz_relationship():
    rng = np.random.default_rng(15)
    x = random_complex(rng, (5000, 4, 5))
    x[:, 3] = 2 * x[:, 0] + 0.1 * x[:, 3]
    result = conditional_products(x)
    assert result["valid_conditional"].all()
    assert result["z_power_predictable_from_h"].max() < 0.003
    assert result["p_power_predictable_from_h"].max() < 0.003
    assert result["partial_coherence_squared_pz_given_h"].min() > 0.99
    np.testing.assert_allclose(
        abs(result["residual_z_p_transfer_ratio"]), 1, atol=0.003
    )


def test_projecting_real_seismic_horizontal_motion_is_not_tilt_estimation():
    rng = np.random.default_rng(16)
    x = random_complex(rng, (5000, 4, 5))
    # A common seismic mode, with independent sensor noise, and NO tilt term.
    signal = random_complex(rng, (5000, 5))
    x[:, 0] = signal + 0.05 * x[:, 0]
    x[:, 1] = signal + 0.05 * x[:, 1]
    x[:, 3] = 2 * signal + 0.1 * x[:, 3]
    result = conditional_products(x)
    assert result["z_power_predictable_from_h"].min() > 0.99
    assert abs(result["residual_z_p_transfer_ratio"]).max() < 0.01
    # Removing the H projection would destroy >99% of the desired transfer,
    # even though none of it was caused by tilt.


def test_schur_complement_matches_direct_covariance_and_is_unit_invariant():
    rng = np.random.default_rng(17)
    x = random_complex(rng, (100, 4, 3))
    x[:, 3] += 0.5 * x[:, 0] + 0.4 * x[:, 1]
    result = conditional_products(x)
    s = result["cross_matrix"][0]
    remaining = [0, 3]
    expected = (
        s[np.ix_(remaining, remaining)]
        - s[np.ix_(remaining, [1, 2])]
        @ np.linalg.inv(s[1:3, 1:3])
        @ s[np.ix_([1, 2], remaining)]
    )
    coh = abs(expected[1, 0]) ** 2 / (expected[0, 0].real * expected[1, 1].real)
    assert result["partial_coherence_squared_pz_given_h"][0] == pytest.approx(coh)
    x *= np.array([1e6, 1e-9, 1e9, 1e-6])[None, :, None]
    rescaled = conditional_products(x)
    for key in (
        "partial_coherence_squared_pz_given_h",
        "z_power_predictable_from_h",
        "residual_z_p_transfer_ratio",
    ):
        np.testing.assert_allclose(result[key], rescaled[key], rtol=1e-12, atol=1e-14)


def test_collinear_horizontals_are_masked_not_regularized_into_a_result():
    rng = np.random.default_rng(18)
    x = random_complex(rng, (50, 4, 3))
    x[:, 2] = x[:, 1]
    result = conditional_products(x)
    assert not result["valid_conditional"].any()
    assert np.isnan(result["partial_coherence_squared_pz_given_h"]).all()
    assert np.isnan(result["residual_z_p_transfer_ratio"]).all()


def test_too_few_averages_and_bad_responses_are_refused():
    with pytest.raises(ValueError, match=">=6 windows"):
        conditional_products(np.ones((2, 4, 20)))
    samples = np.ones((4, 20))
    f = np.fft.rfftfreq(20, 0.5)
    mask = f > 0
    with pytest.raises(ValueError, match="nonzero responses"):
        corrected_fourier(samples, np.zeros((4, mask.sum())), 2, mask)


def test_finite_sample_horizontal_predictability_is_not_tilt_power():
    rng = np.random.default_rng(19)
    null6 = conditional_products(random_complex(rng, (6, 4, 6000)))
    null24 = conditional_products(random_complex(rng, (24, 4, 6000)))
    # Independent random channels still have apparent in-sample predictability.
    assert 0.30 < np.mean(null6["z_power_predictable_from_h"]) < 0.36
    assert 0.07 < np.mean(null24["z_power_predictable_from_h"]) < 0.10


@pytest.mark.parametrize(
    "hours,starts", [(5, ["2013-06-10"]), (25, ["2013-06-10"]), (6, ["2013-06-10"] * 3)]
)
def test_pilot_size_is_bounded_before_any_file_read(hours, starts):
    from scripts.pilot_microseism_horizontal import run

    with pytest.raises(ValueError, match="6..24 hourly"):
        run(None, None, None, "XX.SYN", 0, starts, hours, ["H1", "H2"], None)


def test_pilot_window_csv_preserves_sample_grid_refusal_reason(tmp_path):
    import csv

    from compy_pressure.godin_pilot import write_csv
    from scripts.pilot_microseism_horizontal import window_rows

    rows = window_rows(
        [
            {"start": "00", "status": "retained", "actual_start": "00.001"},
            {"start": "23", "status": "sample_grid_refused", "reason": "offset"},
        ]
    )
    path = tmp_path / "windows.csv"
    write_csv(path, rows)
    with path.open() as stream:
        saved = list(csv.DictReader(stream))
    assert saved[1]["reason"] == "offset"
    assert saved[0]["reason"] == ""
