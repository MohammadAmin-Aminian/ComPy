"""Time-resolved spectra, coherence and tilt direction from plain arrays.

These functions take NumPy arrays sampled on a common grid (no ObsPy), so
they run on cleaned MiniSEED read with any reader, on synthetic data and in
tests. They feed the diagnostic figures of ``compy.plotting.diagnostics``:

* :func:`window_spectra` cuts records into windows and returns, per window,
  Welch auto- and cross-spectra of every channel pair (the same estimator,
  segment length, overlap and Hann window as ``compy.compliance.spectra``);
* :func:`coherence` and :func:`admittance` derive amplitude coherence and
  |S_xy|/S_xx from them;
* :func:`tilt_direction` finds, per window, the horizontal azimuth whose
  rotated component is most coherent with the vertical (the ATaCR/OBStools
  approach of Bell et al. 2015), with that coherence and admittance;
* :func:`band_average` reduces a spectrogram to one value per window.

Cross-spectra follow SciPy: ``S_xy = E[conj(X) Y]``, so S_xy/S_xx is the
transfer function from X to Y. Coherence is amplitude coherence
|S_xy|/sqrt(S_xx S_yy), as everywhere in ComPy.
"""

from itertools import combinations

import numpy as np
from scipy import signal

from compy.numerics import integer, positive


def _windows(x, width, step):
    view = np.lib.stride_tricks.sliding_window_view(x, width)[::step]
    return view


def window_spectra(channels, fs, window_s, nseg, *, step_s=None, t0=0.0):
    """Welch spectra of every channel and channel pair in each window.

    Args:
        channels (dict): ``{name: 1-D array}`` on a common sample grid, e.g.
            ``{"Z": z, "P": p, "1": h1, "2": h2}``. Non-finite samples make
            the windows that contain them NaN.
        fs (float): Sampling rate [Hz].
        window_s (float): Window length [s].
        nseg (int): Welch segment length [samples], 50 % overlap, Hann.
        step_s (float or None): Window step [s]; defaults to ``window_s``.
        t0 (float): Time of the first sample [s], e.g. a POSIX time.

    Returns:
        dict: ``times`` (window centres [s]), ``f`` [Hz], and spectra keyed
        ``"ZZ"``, ``"PP"``, ``"PZ"`` (= E[conj(P) Z]) ... with shape
        (n_windows, n_freqs). Pairs follow the order of ``channels``.
    """
    fs = positive(fs, "sampling rate")
    nseg = integer(nseg, "nseg", minimum=2)
    names = list(channels)
    data = {k: np.asarray(v, float) for k, v in channels.items()}
    n = min(len(v) for v in data.values())
    width = int(round(positive(window_s, "window_s") * fs))
    step = int(round((window_s if step_s is None else positive(step_s, "step_s")) * fs))
    if width < nseg:
        raise ValueError("window shorter than one Welch segment")
    if n < width:
        raise ValueError("records shorter than one window")
    kw = dict(fs=fs, nperseg=nseg, noverlap=nseg // 2, window="hann",
              detrend="constant", axis=-1)  # fmt: skip
    cut = {k: _windows(v[:n], width, step) for k, v in data.items()}
    n_win = len(next(iter(cut.values())))
    out = {"times": t0 + (np.arange(n_win) * step + width / 2) / fs}
    for k in names:
        f, out[k + k] = signal.welch(cut[k], **kw)
    for a, b in combinations(names, 2):
        _, out[a + b] = signal.csd(cut[a], cut[b], **kw)
    out["f"] = f
    bad = np.zeros(n_win, bool)
    for v in cut.values():
        bad |= ~np.all(np.isfinite(v), axis=1)
    for key, value in out.items():
        if key not in ("times", "f") and bad.any():
            value[bad] = np.nan
    return out


def _cross(spec, a, b):
    if a + b in spec:
        return spec[a + b]
    return np.conj(spec[b + a])


def coherence(spec, a, b):
    """Amplitude coherence |S_ab|/sqrt(S_aa S_bb) per window and frequency."""
    sab, saa, sbb = _cross(spec, a, b), spec[a + a], spec[b + b]
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.clip(np.abs(sab) / np.sqrt(saa * sbb), 0, 1)


def admittance(spec, a, b):
    """|S_ab|/S_aa: the gain of the transfer function from ``a`` to ``b``."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.abs(_cross(spec, a, b)) / spec[a + a]


def phase(spec, a, b):
    """Phase of the transfer function from ``a`` to ``b`` [deg]."""
    return np.degrees(np.angle(_cross(spec, a, b)))


def band_average(f, values, band, reducer=np.nanmedian):
    """Reduce ``values`` (windows x freqs) over ``band`` = (f_lo, f_hi) [Hz]."""
    f = np.asarray(f)
    sel = (f >= band[0]) & (f <= band[1])
    if not sel.any():
        raise ValueError("band contains no frequency bins")
    return reducer(np.asarray(values)[..., sel], axis=-1)


def tilt_direction(spec, band, *, vertical="Z", h1="1", h2="2", n_azimuths=360):
    """Azimuth of the horizontal that best predicts the vertical, per window.

    The horizontal rotated by azimuth phi (clockwise from H1) is
    ``H = cos(phi) H1 + sin(phi) H2``. Rotation keeps the horizontal power,
    so the band-mean cross-spectral power |S_HZ(phi)|² peaks at the
    direction whose motion leaks into Z: the tilt azimuth (the principal
    direction of the vector (S_1Z, S_2Z)). Coherence with Z is reported at
    that azimuth, as is the admittance |S_HZ|/S_HH, which for acceleration
    channels is the sine of the tilt angle. The 180° ambiguity is resolved
    by the sign of Re(S_HZ) (positive: Z in phase with H). ATaCR/OBStools
    pick the maximum coherence instead (Bell et al. 2015); with quiet
    horizontals coherence is flat in azimuth, which this estimator avoids.

    Args:
        spec (dict): Output of :func:`window_spectra` with Z, 1 and 2.
        band (tuple): (f_lo, f_hi) [Hz], e.g. (0.005, 0.035) as in ATaCR.
        n_azimuths (int): Grid size over [0, 180).

    Returns:
        dict: ``azimuth`` [deg, 0-360], ``coherence``, ``admittance``, and
        ``coherence_grid`` (windows x azimuths), per window.
    """
    f = spec["f"]
    sel = (f >= band[0]) & (f <= band[1])
    if not sel.any():
        raise ValueError("band contains no frequency bins")
    phi = np.radians(np.arange(n_azimuths) * 180.0 / n_azimuths)
    c, s = np.cos(phi)[None, :, None], np.sin(phi)[None, :, None]
    s11 = spec[h1 + h1][:, None, sel]
    s22 = spec[h2 + h2][:, None, sel]
    s12 = _cross(spec, h1, h2)[:, None, sel]
    s1z = _cross(spec, h1, vertical)[:, None, sel]
    s2z = _cross(spec, h2, vertical)[:, None, sel]
    szz = spec[vertical + vertical][:, None, sel]
    shh = c * c * s11 + s * s * s22 + 2 * c * s * s12.real
    shz = c * s1z + s * s2z
    with np.errstate(invalid="ignore", divide="ignore"):
        power = np.nanmean(np.abs(shz) ** 2, axis=-1)
        coh = np.nanmean(np.abs(shz) / np.sqrt(shh * szz), axis=-1)
        adm = np.nanmean(np.abs(shz) / shh, axis=-1)
        sign = np.sign(np.nanmean(shz.real, axis=-1))
    best = np.argmax(np.where(np.isfinite(power), power, -1), axis=1)
    rows = np.arange(len(best))
    azimuth = np.degrees(phi[best])
    azimuth = np.where(sign[rows, best] < 0, azimuth + 180.0, azimuth) % 360.0
    return dict(
        azimuth=azimuth,
        coherence=coh[rows, best],
        admittance=adm[rows, best],
        coherence_grid=coh,
    )
