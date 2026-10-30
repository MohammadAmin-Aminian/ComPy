"""Independent count-domain Welch audit; never writes scientific products."""

import numpy as np
from scipy.signal import csd


def count_domain_products(arrays, rate, response, weights, band, pressure_sign=-1):
    """Form count cross-spectra first, then correct/project the spectral matrix.

    ``response(f)`` returns full complex counts/Pa and counts/(m/s) responses
    in pressure/H1/H2/Z order. ``weights`` include the upward-Z orientation.
    This uses SciPy CSD, independently of the production segment FFT functions.
    The response/model metadata are still shared, so this is not field certification.
    """
    x = np.asarray(arrays, dtype=float)
    b = np.asarray(weights, dtype=float)
    n = round(200 * rate)
    if (
        x.ndim != 2
        or x.shape[0] != 4
        or n < 4
        or x.shape[1] < 2 * n
        or not np.isfinite(x).all()
        or b.shape != (3,)
        or not np.isfinite(b).all()
        or not np.isclose(np.linalg.norm(b), 1)
        or pressure_sign not in (-1, 1)
        or not 0 < band[0] < band[1] < rate / 2
    ):
        raise ValueError("Invalid synchronized counts, projection or band")
    options = dict(
        fs=rate,
        window="hann",
        nperseg=n,
        noverlap=0,
        detrend="constant",
        scaling="density",
        return_onesided=True,
    )
    frequency, _ = csd(x[0], x[0], **options)
    mask = (frequency >= band[0]) & (frequency <= band[1])
    f = frequency[mask]
    r = np.asarray(response(f))
    if r.shape != (4, len(f)) or not np.all(np.isfinite(r) & (abs(r) > 0)):
        raise ValueError("Invalid response matrix")
    matrix = np.empty((4, 4, len(f)), complex)
    for i in range(4):
        for j in range(i + 1):
            # scipy.csd(u,v)=conj(U)*V: reverse indices for Xi*conj(Xj).
            _, cross = csd(x[j], x[i], **options)
            matrix[i, j] = cross[mask] / (r[i] * r[j].conj())
            matrix[j, i] = matrix[i, j].conj()
    pp = matrix[0, 0].real
    velocity = np.einsum("i,ijf,j->f", b, matrix[1:, 1:], b).real
    omega = 2 * np.pi * f
    aa = omega**2 * velocity
    pa = pressure_sign * (-1j * omega) * np.einsum("if,i->f", matrix[0, 1:], b)
    if not np.all((pp > 0) & (aa > 0)):
        raise ValueError("Nonpositive projected spectra")
    return {
        "frequency_hz": f,
        "pressure_psd": pp,
        "acceleration_psd": aa,
        "cross_pa": pa,
        "psd_ratio_pa": np.sqrt(pp / aa),
        "h1_pa": pa / aa,
        "coherence_squared": np.clip(abs(pa) ** 2 / (pp * aa), 0, 1),
        "phase_deg": np.angle(pa, deg=True),
    }
