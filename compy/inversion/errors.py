"""Errors of a stacked compliance estimate, measured from its own windows.

A compliance curve is a ratio of spectra stacked over thousands of windows.
Its Bendat & Piersol error is a fraction of a per cent, but that number is not
the whole error and it is not independent between frequencies:

* Hann windows with 50 % overlap correlate neighbouring bins (about 0.5).
* The sea state, the instrument and the window selection change over weeks,
  which moves the whole curve up or down and tilts it.

``stack_covariance`` measures both by resampling the windows (bootstrap): once
window by window, which gives the random error with its correlation between
bins, and once in blocks of days, which adds what varies more slowly. The slow
part is kept as two terms, a common level and a tilt in log frequency, because
a year of data has too few independent blocks to estimate more.

``default_covariance`` builds the same structure from a plain list of errors,
for compliance curves that come without their windows.

Measured on 172 stations of the 2026 archive (docs/inversion_overhaul_2026-11):
the bootstrap error is 1.36 times Bendat & Piersol (median), the correlation of
neighbouring bins 0.48 and of bins two apart 0.11, and the slow terms add a
level of 0.3 % and a tilt of 0.5 %.

All quantities are SI. Covariances returned here are relative (divided by the
outer product of compliance), so they can be scaled to any gain.
"""

from __future__ import annotations

import numpy as np

from compy.numerics import (
    FREE_AIR_GRADIENT,
    GRAVITATIONAL_CONSTANT,
    GRAVITY,
    positive,
    wavenumber,
)

#: Ratio of the bootstrap error to the Bendat & Piersol error (archive median).
BOOTSTRAP_OVER_BENDAT_PIERSOL = 1.36
#: Correlation of the stacked estimate between bins one and two apart
#: (Hann, 50 % overlap; archive medians).
BIN_CORRELATION = (0.48, 0.11)
#: Slow level and tilt (relative standard deviations) assumed when the windows
#: are not available (archive medians of the measured values).
DEFAULT_SLOW = dict(level=0.003, tilt=0.005)


def log_frequency_ramp(freqs):
    """Linear function of ln(frequency) from -0.5 at the lowest to +0.5 at the highest."""
    x = np.log(np.asarray(freqs, float))
    if len(x) < 2 or x[-1] == x[0]:
        return np.zeros_like(x)
    return (x - 0.5 * (x[0] + x[-1])) / (x[-1] - x[0])


def attraction_ratio(freqs, depth):
    """Wave attraction per unit seafloor pressure, a_w/p [m/s**2/Pa].

    a_w = 2 pi G rho_w h_w exp(-kH) with h_w = p cosh(kH)/(rho_w g)
    (Crawford et al. 1998, as quoted by Zha et al. 2014, eq. 2).
    """
    kh = wavenumber(2 * np.pi * np.asarray(freqs, float), depth) * depth
    return 2 * np.pi * (GRAVITATIONAL_CONSTANT / GRAVITY) * (1 + np.exp(-2 * kh)) / 2


def compliance_from_sums(freqs, depth, zp_sum, pp_sum):
    """Gravity-corrected normalized compliance from summed S_pz and S_pp.

    Equal to ``compy.compliance.gravity_corrected_compliance`` of the stack;
    ``zp_sum`` and ``pp_sum`` may carry leading axes (bootstrap replicates).
    """
    f = np.asarray(freqs, float)
    omega2 = (2 * np.pi * f) ** 2
    k = wavenumber(2 * np.pi * f, depth)
    ratio = np.abs(zp_sum) / pp_sum
    return k * (omega2 * ratio + attraction_ratio(f, depth)) / (omega2 + FREE_AIR_GRADIENT)


def _level_tilt(dev, ramp):
    """Least-squares level and tilt coefficients of each row of ``dev``."""
    basis = np.column_stack([np.ones_like(ramp), ramp])
    coef, *_ = np.linalg.lstsq(basis, dev.T, rcond=None)
    return coef  # (2, replicates)


def stack_covariance(freqs, pp, zp, starts, depth, *, gain=1.0, block_days=(30.0, 5.0, 1.0),
                     min_blocks=6, n_boot=400, seed=0):  # fmt: skip
    """Relative covariance of stacked compliance from its windows.

    Args:
        freqs (array): Frequencies [Hz] of the band, length F.
        pp (array): Pressure autospectra of the accepted windows, shape (W, F),
            before gain.
        zp (array): Complex cross-spectra S_pz of the same windows, (W, F).
        starts (array): Window start times [s], length W.
        depth (float): Water depth [m].
        gain (float): Pressure gain applied to the stack.
        block_days (tuple): Block lengths tried in turn for the slow terms;
            the first that gives at least ``min_blocks`` blocks is used.
        n_boot (int): Bootstrap replicates.
        seed (int): Random seed.

    Returns:
        dict: ``eta`` (stacked compliance), ``cov_rel`` (F x F relative
        covariance: window bootstrap plus slow level and tilt), ``random_rel``
        (window bootstrap alone), ``sd_rel`` (its standard deviations),
        ``slow_level`` and ``slow_tilt`` (relative standard deviations of the
        slow terms), ``block_days`` and ``n_blocks`` used (None and 0 when the
        record is too short), ``corr_adjacent``, and ``half_level``,
        ``half_tilt`` (log ratio of the first to the second half of the
        record: its mean and its change across the band).

    Raises:
        ValueError: With fewer than 20 windows or mismatched shapes.
    """
    f = np.asarray(freqs, float)
    pp = np.asarray(pp, float) * positive(gain, "gain") ** 2
    zp = np.asarray(zp, complex) * gain
    t = np.asarray(starts, float)
    if pp.shape != zp.shape or pp.shape != (len(t), len(f)):
        raise ValueError("pp and zp must have shape (windows, frequencies)")
    if len(t) < 20:
        raise ValueError("at least 20 windows are needed for a bootstrap covariance")
    depth = positive(depth, "water depth")
    rng = np.random.default_rng(seed)
    eta = compliance_from_sums(f, depth, zp.sum(0), pp.sum(0))
    ramp = log_frequency_ramp(f)

    def replicate(groups):
        ids, inverse = np.unique(groups, return_inverse=True)
        gz = np.zeros((len(ids), len(f)), complex)
        gp = np.zeros((len(ids), len(f)))
        np.add.at(gz, inverse, zp)
        np.add.at(gp, inverse, pp)
        w = rng.multinomial(len(ids), np.full(len(ids), 1 / len(ids)), size=n_boot)
        return compliance_from_sums(f, depth, w @ gz, w @ gp) / eta - 1, len(ids)

    dev, _ = replicate(np.arange(len(t)))
    random_rel = np.atleast_2d(np.cov(dev, rowvar=False))
    sd = np.sqrt(np.diag(random_rel))
    corr = random_rel / np.outer(sd, sd)
    coef_random = _level_tilt(dev, ramp)
    slow_level = slow_tilt = 0.0
    used, n_blocks = None, 0
    for days in block_days:
        groups = np.floor((t - t.min()) / (86400.0 * days)).astype(int)
        if len(np.unique(groups)) >= min_blocks:
            dev_b, n_blocks = replicate(groups)
            coef = _level_tilt(dev_b, ramp)
            slow_level = float(np.sqrt(max(np.var(coef[0]) - np.var(coef_random[0]), 0.0)))
            slow_tilt = float(np.sqrt(max(np.var(coef[1]) - np.var(coef_random[1]), 0.0)))
            used = float(days)
            break
    cov_rel = random_rel + slow_level**2 + slow_tilt**2 * np.outer(ramp, ramp)
    first = t < np.median(t)
    half = np.log(compliance_from_sums(f, depth, zp[first].sum(0), pp[first].sum(0))
                  / compliance_from_sums(f, depth, zp[~first].sum(0), pp[~first].sum(0)))  # fmt: skip
    half_tilt = float(_level_tilt(half[None, :], ramp)[1, 0]) if len(f) > 1 else 0.0
    return dict(
        eta=eta, cov_rel=cov_rel, random_rel=random_rel, sd_rel=sd,
        slow_level=slow_level, slow_tilt=slow_tilt, block_days=used, n_blocks=int(n_blocks),
        corr_adjacent=float(np.mean(np.diag(corr, 1))) if len(f) > 1 else np.nan,
        half_level=float(np.mean(half)), half_tilt=half_tilt,
    )  # fmt: skip


def default_covariance(freqs, eta, sigma, *, inflate=BOOTSTRAP_OVER_BENDAT_PIERSOL,
                       correlation=BIN_CORRELATION, slow_level=DEFAULT_SLOW["level"],
                       slow_tilt=DEFAULT_SLOW["tilt"]):  # fmt: skip
    """Relative covariance for a compliance curve that comes without its windows.

    The random part is ``inflate`` times the given errors, correlated between
    neighbouring bins as a Hann-window Welch stack is. This assumes uniformly
    spaced Welch bins; pass ``correlation=()`` for data that were averaged over
    frequency. Slow level and tilt terms are the archive medians.

    Args:
        freqs, eta, sigma (array): Frequencies [Hz], compliance and its
            one-sigma error [1/Pa].
        inflate (float): Factor on ``sigma``.
        correlation (tuple): Correlation at bin lags 1, 2, ...
        slow_level, slow_tilt (float): Relative standard deviations.

    Returns:
        ndarray: Relative covariance, shape (F, F).
    """
    f, eta, sigma = (np.asarray(v, float) for v in (freqs, eta, sigma))
    if not (f.shape == eta.shape == sigma.shape) or f.ndim != 1:
        raise ValueError("freqs, eta and sigma must be vectors of equal length")
    rel = inflate * sigma / eta
    n = len(f)
    corr = np.eye(n)
    for lag, c in enumerate(correlation, start=1):
        if lag < n:
            corr += c * (np.eye(n, k=lag) + np.eye(n, k=-lag))
    ramp = log_frequency_ramp(f)
    return np.outer(rel, rel) * corr + slow_level**2 + slow_tilt**2 * np.outer(ramp, ramp)


def usable_bins(sd_rel, factor=2.5):
    """Bins of a band whose scatter shows a healthy estimate; the rest are cut from the edges.

    Towards the infragravity cut-off the pressure signal no longer reaches the
    seafloor: the estimate turns noisy and, with noise now on the pressure
    channel, biased low. The bootstrap scatter shows where this starts well
    before the coherence gate does (at shelf station 7D.J09B it rises from
    0.2 % to 1-6 % above 35 mHz while the stacked coherence stays above 0.8).
    Bins are kept from the quietest bin outwards until the scatter exceeds
    ``factor`` times the median scatter of the quieter half of the band.

    Args:
        sd_rel (array): Relative standard deviations from ``stack_covariance``.
        factor (float): Limit in units of that median.

    Returns:
        ndarray of bool: True for the bins to keep (one contiguous run).
    """
    sd = np.asarray(sd_rel, float)
    limit = factor * np.median(np.sort(sd)[: max(len(sd) // 2, 1)])
    best = int(np.argmin(sd))
    lo = hi = best
    while lo > 0 and sd[lo - 1] <= limit:
        lo -= 1
    while hi < len(sd) - 1 and sd[hi + 1] <= limit:
        hi += 1
    keep = np.zeros(len(sd), bool)
    keep[lo : hi + 1] = True
    return keep
