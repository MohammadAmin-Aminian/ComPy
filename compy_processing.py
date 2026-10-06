"""Welch spectra and window selection for response-corrected OBS records."""

import numpy as np
from scipy import signal

from compy_numerics import integer, positive
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


def estimate(stream, depth, gain, hours, f_min, f_max, beta, nseg=4096):
    """Return selected compliance spectra and matching window metadata."""
    from compy import gravitational_attraction
    from compy_numerics import wavenumber, compliance_uncertainty

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
            fc1 = np.sqrt(9.8 / (2 * np.pi * 0.5 * depth))
            fc2 = np.sqrt(9.8 / (2 * np.pi * 2 * depth))
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
        k = wavenumber(2 * np.pi * f, depth)
        gravity_ratio, _, _ = gravitational_attraction(np.sqrt(calibrated_pp), depth, f)
        ratio = np.zeros_like(f)
        np.divide(
            np.sqrt(zz), np.sqrt(calibrated_pp), out=ratio, where=calibrated_pp > 0
        )
        omega2 = (2 * np.pi * f) ** 2
        # Convert wave-attraction acceleration/pressure to displacement/pressure
        # before adding it to the measured displacement ratio (SI units).
        compliance = k * coh * (omega2 * ratio + gravity_ratio) / (omega2 + 3.07e-6)
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
    averages = 1 + (len(accepted[0][0]) - nseg) // (nseg - nseg // 2)
    theory = compliance_uncertainty(
        curves[0], np.asarray(coherences)[0, keep], averages
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
