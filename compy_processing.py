"""Welch spectra and window selection for response-corrected OBS records."""

from dataclasses import dataclass, field

import numpy as np
from scipy import signal

from compy_numerics import (
    FREE_AIR_GRADIENT,
    GRAVITY,
    integer,
    positive,
    welch_averages,
)
from compy_streams import cut_stream_with_overlap, validate_stream


def channel_pair(stream):
    """Select one vertical and one pressure channel, independently of order."""
    vertical, pressure = stream.select(component="Z"), stream.select(component="H")
    if len(vertical) != 1 or len(pressure) != 1:
        raise ValueError(
            "exactly one vertical (*Z) and pressure (*H, e.g. BDH) trace are required"
        )
    validate_stream(vertical + pressure)
    return vertical[0], pressure[0]


def spectra(stream, nseg=4096):
    """Consistent mean-Welch auto/cross spectra and amplitude coherence.

    Means are used for all spectra so |S_zp|^2/(S_zz*S_pp) follows the same
    estimator, without mixing median autospectra and mean cross-spectra.
    The cross-spectrum is S_pz = E[conj(P) Z], so S_pz / S_pp is Z/P.
    """
    z, p = channel_pair(stream)
    nseg = integer(nseg, "nseg", minimum=2)
    if len(z) != len(p) or z.stats.starttime != p.stats.starttime:
        raise ValueError(
            "vertical and pressure traces must have identical time coverage"
        )
    if len(z) < nseg:
        raise ValueError("processing window is shorter than nseg")
    kw = dict(
        fs=z.stats.sampling_rate,
        nperseg=nseg,
        noverlap=nseg // 2,
        window="hann",
        detrend="constant",
    )
    f, zz = signal.welch(z.data, **kw)
    _, pp = signal.welch(p.data, **kw)
    _, zp = signal.csd(p.data, z.data, **kw)
    coh = np.zeros_like(zz)
    np.divide(np.abs(zp), np.sqrt(zz * pp), out=coh, where=(zz * pp) > 0)
    return f, zz, pp, zp, np.clip(coh, 0, 1)


def gravity_corrected_compliance(f, depth, displacement_per_pressure, pressure_amp):
    """Normalized compliance from a measured |Z/P| ratio, with gravity terms.

    A seismometer measures a_d - a_w + a_e: the wave attraction a_w opposes
    the loading acceleration and the free-air term a_e follows it (Crawford et
    al. 1998; Zha 2014, eq. 2.2). The true displacement is therefore
    (omega² |Z/P| + a_w/p) / (omega² + FREE_AIR_GRADIENT) per pascal.

    Args:
        f (array): Frequencies [Hz].
        depth (float): Water depth [m].
        displacement_per_pressure (array): Coherent |Z/P| [m/Pa], already
            multiplied by coherence if it comes from an autospectral ratio.
        pressure_amp (array): Pressure amplitude used for a_w [Pa/sqrt(Hz)].

    Returns:
        array: Normalized compliance k Z/P [1/Pa].
    """
    from compy import gravitational_attraction
    from compy_numerics import wavenumber

    f = np.asarray(f, float)
    k = wavenumber(2 * np.pi * f, depth)
    gravity_ratio, _, _ = gravitational_attraction(pressure_amp, depth, f)
    omega2 = (2 * np.pi * f) ** 2
    return (
        k
        * (omega2 * displacement_per_pressure + gravity_ratio)
        / (omega2 + FREE_AIR_GRADIENT)
    )


@dataclass(frozen=True)
class SelectionCriteria:
    """Per-window quality gates evaluated in the compliance band.

    Coherence is amplitude coherence |S_zp|/sqrt(S_zz S_pp). PSD limits are
    band means in dB of calibrated pressure [Pa²/Hz] and vertical acceleration
    [(m/s²)²/Hz]; ``None`` disables a limit.

    Attributes:
        min_median_coherence (float): Lower bound on the band median coherence.
        pressure_psd_db (tuple or None): (low, high) pressure PSD limits.
        acceleration_psd_db (tuple or None): (low, high) acceleration limits.
        name (str): Label written to run records.
    """

    min_median_coherence: float = 0.8
    pressure_psd_db: tuple = (10.0, 50.0)
    acceleration_psd_db: tuple = (-170.0, -150.0)
    name: str = "gji2025"

    def check(self, f, zz, pp, coh, band):
        """Return the reason a window fails, or ``None`` if it passes."""
        if np.median(coh[band]) < self.min_median_coherence:
            return "coherence"
        if self.pressure_psd_db is not None:
            level = 10 * np.log10(np.mean(pp[band]))
            if not self.pressure_psd_db[0] <= level <= self.pressure_psd_db[1]:
                return "pressure_psd"
        if self.acceleration_psd_db is not None:
            level = 10 * np.log10(np.mean(zz[band] * (2 * np.pi * f[band]) ** 4))
            low, high = self.acceleration_psd_db
            if not low <= level <= high:
                return "acceleration_psd"
        return None


#: Aminian et al. (2025, GJI) section 4.2: median coherence > 0.8 (0.9 for
#: better stations), pressure PSD 10–50 dB, acceleration PSD -170 to -150 dB.
GJI2025 = SelectionCriteria()


@dataclass
class ComplianceResult:
    """Stacked compliance estimate and the windows that produced it.

    Attributes:
        freqs (ndarray): Frequencies [Hz] of the returned band.
        compliance (ndarray): Normalized compliance from stacked spectra [1/Pa].
        uncertainty (ndarray): Bendat & Piersol one-sigma error [1/Pa].
        coherence (ndarray): Amplitude coherence of the stacked spectra.
        phase_deg (ndarray): Phase of Z/P [degrees]; ±180 means downward
            displacement under positive pressure.
        window_mean (ndarray): Mean of per-window compliance curves [1/Pa].
        window_std (ndarray): Standard deviation across windows [1/Pa].
        n_windows (int): Accepted windows.
        n_averages (float): Effective independent spectral averages.
        window_starts (list): UTCDateTime start of each accepted window.
        rejected (dict): Rejected window counts by reason.
        band (tuple): (f_min, f_max) band used for window gates [Hz].
        gain (float): Pressure gain applied.
        depth (float): Water depth [m].
        criteria (SelectionCriteria): Gates used.
        all_freqs (ndarray): Full frequency axis [Hz].
        mean_coherence_all (ndarray): Mean per-window coherence, full axis.
        stack (dict): Stacked S_zz, S_pp, S_pz on the full axis (SI units).
    """

    freqs: np.ndarray
    compliance: np.ndarray
    uncertainty: np.ndarray
    coherence: np.ndarray
    phase_deg: np.ndarray
    window_mean: np.ndarray
    window_std: np.ndarray
    n_windows: int
    n_averages: float
    window_starts: list
    rejected: dict
    band: tuple
    gain: float
    depth: float
    criteria: SelectionCriteria
    all_freqs: np.ndarray = field(repr=False)
    mean_coherence_all: np.ndarray = field(repr=False)
    stack: dict = field(repr=False)


def stacked_compliance(
    windows,
    depth,
    gain,
    band,
    criteria=GJI2025,
    nseg=1024,
    *,
    output_band=None,
    gravity_correction=True,
):
    """Compliance from cross-spectra stacked over quality-selected windows.

    Each window must contain one *Z displacement trace [m] and one *H
    pressure trace [Pa, before gain]. Windows passing ``criteria`` in
    ``band`` are stacked: S_zz, S_pp and complex S_pz are averaged, then
    eta = k |S_pz| / S_pp, gravity corrected. The uncertainty follows Bendat
    & Piersol (1993) with n the effective number of independent Welch
    averages over all accepted windows (non-overlapping windows assumed).

    Args:
        windows (list of Stream): Equal-length, non-overlapping windows.
        depth (float): Water depth [m].
        gain (float): Factor multiplying nominal pressure.
        band (tuple): (f_min, f_max) [Hz] used for the window gates.
        criteria (SelectionCriteria): Window gates.
        nseg (int): Welch segment length in samples.
        output_band (tuple or None): Frequency range returned; defaults to band.
        gravity_correction (bool): Apply wave-attraction and free-air terms.

    Returns:
        ComplianceResult: Stacked estimate and per-window statistics.

    Raises:
        ValueError: If no window passes or the band holds no frequency bins.
    """
    from compy_numerics import compliance_uncertainty, wavenumber

    depth, gain = positive(depth, "water depth"), positive(gain, "gain factor")
    f_min, f_max = band
    if not 0 < f_min < f_max:
        raise ValueError("band must satisfy 0 < f_min < f_max")
    rejected = {"coherence": 0, "pressure_psd": 0, "acceleration_psd": 0}
    zz_sum = pp_sum = zp_sum = coh_sum = None
    curves, starts, n_samples = [], [], None
    for window in windows:
        f, zz, pp, zp, coh = spectra(window, nseg)
        sel = (f >= f_min) & (f <= f_max)
        if not sel.any():
            raise ValueError("band contains no spectral bins; increase nseg")
        pp, zp = pp * gain**2, zp * gain
        reason = criteria.check(f, zz, pp, coh, sel)
        if reason is not None:
            rejected[reason] += 1
            continue
        if n_samples is None:
            n_samples = len(window[0])
        elif len(window[0]) != n_samples:
            raise ValueError("all windows must have the same length")
        if zz_sum is None:
            zz_sum, pp_sum = np.zeros_like(zz), np.zeros_like(pp)
            zp_sum, coh_sum = np.zeros_like(zp), np.zeros_like(coh)
        zz_sum += zz
        pp_sum += pp
        zp_sum += zp
        coh_sum += coh
        ratio = np.zeros_like(f)
        np.divide(np.abs(zp), pp, out=ratio, where=pp > 0)
        curves.append(_compliance(f, depth, ratio, pp, gravity_correction))
        starts.append(window[0].stats.starttime)
    if not curves:
        raise ValueError(f"no windows passed compliance quality thresholds {rejected}")
    n = len(curves)
    szz, spp, spz = zz_sum / n, pp_sum / n, zp_sum / n
    coherence = np.zeros_like(szz)
    np.divide(np.abs(spz), np.sqrt(szz * spp), out=coherence, where=szz * spp > 0)
    coherence = np.clip(coherence, 0, 1)
    ratio = np.zeros_like(szz)
    np.divide(np.abs(spz), spp, out=ratio, where=spp > 0)
    eta = _compliance(f, depth, ratio, spp, gravity_correction)
    averages = n * welch_averages(n_samples, nseg)
    lo, hi = output_band if output_band is not None else band
    keep = (f >= lo) & (f <= hi)
    omega2 = (2 * np.pi * f[keep]) ** 2
    # Only the measured transfer function is random; gravity terms are scaled.
    measured = (
        wavenumber(2 * np.pi * f[keep], depth)
        * omega2
        * ratio[keep]
        / (omega2 + (FREE_AIR_GRADIENT if gravity_correction else 0))
    )
    curves = np.asarray(curves)
    return ComplianceResult(
        freqs=f[keep],
        compliance=eta[keep],
        uncertainty=compliance_uncertainty(measured, coherence[keep], averages),
        coherence=coherence[keep],
        phase_deg=np.degrees(np.angle(spz[keep])),
        window_mean=curves[:, keep].mean(axis=0),
        window_std=curves[:, keep].std(axis=0, ddof=1 if n > 1 else 0),
        n_windows=n,
        n_averages=averages,
        window_starts=starts,
        rejected=rejected,
        band=(f_min, f_max),
        gain=gain,
        depth=depth,
        criteria=criteria,
        all_freqs=f,
        mean_coherence_all=coh_sum / n,
        stack={"S_zz": szz, "S_pp": spp, "S_pz": spz},
    )


def _compliance(f, depth, ratio, pp, gravity_correction):
    if gravity_correction:
        return gravity_corrected_compliance(f, depth, ratio, np.sqrt(pp))
    from compy_numerics import wavenumber

    return wavenumber(2 * np.pi * f, depth) * ratio


def coherent_band(f, coherence, threshold, search):
    """Widest contiguous band inside ``search`` with coherence >= threshold.

    Implements the band definition of Aminian et al. (2025, section 4.2).

    Args:
        f (array): Frequencies [Hz].
        coherence (array): Mean amplitude coherence on ``f``.
        threshold (float): Coherence threshold.
        search (tuple): (low, high) frequencies [Hz] to search within.

    Returns:
        tuple or None: (f_low, f_high) [Hz], or None if no bin qualifies.
    """
    f, coherence = np.asarray(f), np.asarray(coherence)
    inside = (f >= search[0]) & (f <= search[1])
    ok = inside & (coherence >= threshold)
    best, run = None, None
    for i, flag in enumerate(ok):
        if flag:
            run = (run[0], i) if run else (i, i)
            if best is None or run[1] - run[0] > best[1] - best[0]:
                best = run
        else:
            run = None
    if best is None:
        return None
    return float(f[best[0]]), float(f[best[1]])


def _covered_samples(windows, width):
    """Distinct samples covered by possibly overlapping, time-ordered windows."""
    starts = sorted(w[0].stats.starttime for w in windows)
    fs = windows[0][0].stats.sampling_rate
    total = 0
    for a, b in zip(starts, starts[1:]):
        total += min(width, int(round((b - a) * fs)))
    return total + width


def estimate(stream, depth, gain, hours, f_min, f_max, beta, nseg=4096):
    """Return selected compliance spectra and matching window metadata."""
    from compy_numerics import compliance_uncertainty

    depth, gain = positive(depth, "water depth"), positive(gain, "gain factor")
    hours = positive(hours, "time_window")
    if not np.isfinite(f_min) or not np.isfinite(f_max) or not 0 < f_min < f_max:
        raise ValueError("frequency band must satisfy 0 < f_min < f_max")
    # Historical step lengths are retained: beta 5 minutes, main 1 minute.
    step = min(hours * 3600, 300 if beta else 60)
    windows = cut_stream_with_overlap(stream, hours * 3600, hours * 3600 - step)
    if not windows:
        raise ValueError("no complete processing windows")
    selected, coherences, accepted, pressure_psds, vertical_psds = [], [], [], [], []
    for window in windows:
        f, zz, pp, _, coh = spectra(window, nseg)
        band = (f >= f_min) & (f <= f_max)
        secondary = (f >= (0.02 if beta else 0.03)) & (f <= (0.05 if beta else 0.08))
        if not np.any(band) or not np.any(secondary):
            raise ValueError(
                "frequency bands contain no spectral bins; increase nseg or adjust bands"
            )
        calibrated_pp = pp * gain**2
        if beta:
            acc = zz[band] * (2 * np.pi * f[band]) ** 4
            valid = (
                np.median(coh[band]) > 0.9
                and np.mean(coh[band]) > 0.9
                and np.min(coh[band]) > 0.6
                and 100 < np.median(calibrated_pp[band]) < 1e5
                and np.median(calibrated_pp[secondary]) < 1e5
                and np.all(acc > 0)
                and -180 < np.mean(10 * np.log10(acc)) < -110
            )
        else:
            # Historical RR52-specific gates on uncalibrated pressure; see
            # SelectionCriteria/GJI2025 for the published, station-independent set.
            fc1 = np.sqrt(GRAVITY / (2 * np.pi * 0.5 * depth))
            fc2 = np.sqrt(GRAVITY / (2 * np.pi * 2 * depth))
            ig = (f >= fc2) & (f <= fc1)
            if not np.any(ig):
                raise ValueError("depth-dependent IG band contains no spectral bins")
            valid = (
                np.median(coh[band]) > 0.9
                and np.mean(coh[secondary]) > 0.8
                and 0.5 < np.mean(pp[secondary]) < 1
                and np.mean(zz[secondary]) > 1e-16
                and np.mean(zz[ig]) < 1e-11
            )
        if not valid:
            continue
        ratio = np.zeros_like(f)
        np.divide(
            np.sqrt(zz), np.sqrt(calibrated_pp), out=ratio, where=calibrated_pp > 0
        )
        # Coherence scales only the measured ratio: coh*sqrt(zz/pp) = |S_zp|/S_pp.
        # The deterministic gravity term must not be multiplied by coherence.
        compliance = gravity_corrected_compliance(
            f, depth, coh * ratio, np.sqrt(calibrated_pp)
        )
        selected.append(compliance)
        coherences.append(coh)
        accepted.append(window)
        pressure_psds.append(calibrated_pp)
        vertical_psds.append(zz)
    if not selected:
        raise ValueError("no windows passed compliance quality thresholds")
    keep = (f >= (0.001 if beta else 0.005)) & (f <= (0.1 if beta else 0.025))
    curves = np.asarray(selected)[:, keep]
    # Scatter is not the uncertainty on the mean and overlap is not independent.
    scatter = np.std(curves, axis=0) if beta else np.ptp(curves, axis=0)
    # Bendat & Piersol error of the averaged curve. Overlapping processing
    # windows share samples, so only the distinct samples they cover count.
    width = len(accepted[0][0])
    averages = welch_averages(_covered_samples(accepted, width), nseg)
    theory = compliance_uncertainty(
        curves.mean(axis=0), np.asarray(coherences)[:, keep].mean(axis=0), averages
    )
    return (
        curves,
        coherences,
        accepted,
        f[keep],
        f,
        scatter,
        theory,
        pressure_psds,
        vertical_psds,
    )
