"""P/H1/H2/Z conditional spectra: diagnostic projections, NOT tilt correction.

Horizontal seismic motion may be a physical part of the same Rayleigh wave.
Projecting it out can remove the desired P/Z signal even with zero tilt.
"""

import numpy as np
from scipy.signal.windows import hann


def corrected_fourier(samples, responses, rate, frequency_mask, upward_sign=1):
    """One-sided density-normalized coefficients, order P,H1,H2,Z.

    All seismic responses must be counts/displacement. Pressure stays nominal;
    only Z receives the recorded upward-orientation convention.
    """
    samples = np.asarray(samples, dtype=float)
    responses = np.asarray(responses, dtype=complex)
    if samples.ndim != 2 or samples.shape[0] != 4 or samples.shape[1] < 4:
        raise ValueError("Expected four synchronized sample arrays")
    if not np.isfinite(samples).all() or not np.isfinite(rate) or rate <= 0:
        raise ValueError("Samples and sample rate must be finite")
    if upward_sign not in (-1, 1):
        raise ValueError("Upward orientation sign must be +/-1")
    taper = hann(samples.shape[1], sym=False)
    transformed = np.fft.rfft((samples - samples.mean(axis=1, keepdims=True)) * taper)
    f = np.fft.rfftfreq(samples.shape[1], 1 / rate)
    selected = np.asarray(frequency_mask, dtype=bool)
    if selected.shape != f.shape or not selected.any() or np.any(f[selected] <= 0):
        raise ValueError("Select positive bins on the full FFT grid")
    transformed = transformed[:, selected]
    if responses.shape != transformed.shape or np.any(
        ~np.isfinite(responses) | (abs(responses) <= 0)
    ):
        raise ValueError("Four finite nonzero responses must match selected bins")
    scale = np.full(f.shape, 2 / (rate * np.sum(taper**2)))
    scale[0] /= 2
    if samples.shape[1] % 2 == 0:
        scale[-1] /= 2
    coefficients = transformed / responses * np.sqrt(scale[selected])
    coefficients[3] *= upward_sign
    return coefficients


def conditional_products(window_coefficients, maximum_horizontal_condition=1e6):
    """Schur complement of the horizontal block after averaging cross spectra.

    Returns the P/Z coherence, H-predictable fractions of P and Z power, and
    conditional coherence. The complex transfer ratio is (Z minus its linear
    H projection)/P divided by Z/P. It is NOT a pressure-gain correction.
    Ill-conditioned horizontal blocks or vanishing residual powers are masked.
    """
    x = np.asarray(window_coefficients, dtype=complex)
    if x.ndim != 3 or x.shape[0] < 6 or x.shape[1] != 4 or not np.isfinite(x).all():
        raise ValueError("Need >=6 windows of finite P,H1,H2,Z coefficients")
    if (
        not np.isfinite(maximum_horizontal_condition)
        or maximum_horizontal_condition <= 1
    ):
        raise ValueError("Horizontal conditioning limit must exceed one")
    matrix = np.einsum("wif,wjf->fij", x, x.conj()) / len(x)
    powers = np.real(np.diagonal(matrix, axis1=1, axis2=2))
    if np.any(powers <= 0):
        raise ValueError("All four spectral powers must be positive")
    normal = matrix / np.sqrt(powers[:, :, None] * powers[:, None, :])
    hh = normal[:, 1:3, 1:3]
    condition = np.linalg.cond(hh)
    valid = np.isfinite(condition) & (condition < maximum_horizontal_condition)
    # pinv is used only after conditioning screening; invalid bins remain NaN.
    inverse = np.linalg.pinv(hh)
    z_h, h_p = normal[:, 3, 1:3], normal[:, 1:3, 0]
    projection_zp = np.einsum("fi,fij,fj->f", z_h, inverse, h_p)
    z_fraction = np.real(np.einsum("fi,fij,fj->f", z_h, inverse, z_h.conj()))
    p_h = normal[:, 0, 1:3]
    p_fraction = np.real(np.einsum("fi,fij,fj->f", p_h, inverse, h_p))
    rp, rz = 1 - p_fraction, 1 - z_fraction
    residual = normal[:, 3, 0] - projection_zp
    valid &= (rp > 1e-10) & (rz > 1e-10) & (abs(normal[:, 3, 0]) > 1e-12)
    partial = np.full(len(matrix), np.nan)
    ratio = np.full(len(matrix), np.nan + 1j * np.nan)
    partial[valid] = abs(residual[valid]) ** 2 / (rp[valid] * rz[valid])
    ratio[valid] = residual[valid] / normal[valid, 3, 0]
    return {
        "cross_matrix": matrix,
        "coherence_squared_pz": abs(normal[:, 3, 0]) ** 2,
        "horizontal_condition": condition,
        "valid_conditional": valid,
        "p_power_predictable_from_h": np.where(
            valid, np.clip(p_fraction, 0, 1), np.nan
        ),
        "z_power_predictable_from_h": np.where(
            valid, np.clip(z_fraction, 0, 1), np.nan
        ),
        "partial_coherence_squared_pz_given_h": np.clip(partial, 0, 1),
        "residual_z_p_transfer_ratio": ratio,
    }
