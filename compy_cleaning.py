"""Unattended OBS cleaning steps used before compliance estimation.

All functions take ObsPy traces on a common sample grid and return copies.
They follow Aminian et al. (2025, section 4.1) without interactive steps:

* :func:`periodic_template` / :func:`remove_periodic_transient`: stack-and-
  subtract removal of the hourly mass-position glitch (Deen et al. 2017).
* :func:`rotation_coefficients` / :func:`rotate_smoothly`: tilt rotation with
  hourly TiSKitPy angles, smoothed and interpolated per sample so no steps
  are introduced at window boundaries.
* :func:`sta_lta_spans`: local-event time spans from a recursive STA/LTA.
"""

from dataclasses import dataclass

import numpy as np
from obspy import UTCDateTime
from obspy.signal.rotate import rotate2zne
from obspy.signal.trigger import recursive_sta_lta, trigger_onset

from compy_numerics import positive


@dataclass
class PeriodicTemplate:
    """Stacked periodic transient for one trace and time block.

    Attributes:
        period (float): Transient period [s].
        reference (UTCDateTime): Time of template phase zero.
        template (ndarray): Mean transient sampled every ``delta`` seconds.
        delta (float): Template sample interval [s].
        strength (float): Fold-variance peak over its median across trial
            periods; values near 1 mean no periodic signal.
        n_cycles (int): Cycles stacked.
    """

    period: float
    reference: UTCDateTime
    template: np.ndarray
    delta: float
    strength: float
    n_cycles: int

    def evaluate(self, times):
        """Template value at absolute times (seconds since the epoch)."""
        phase = np.mod(np.asarray(times) - self.reference.timestamp, self.period)
        grid = np.arange(len(self.template)) * self.delta
        return np.interp(phase, grid, self.template, period=self.period)


def _fold_variance(x, t, period, n_bins):
    index = (np.mod(t, period) / period * n_bins).astype(int)
    count = np.bincount(index, minlength=n_bins)
    stack = np.bincount(index, weights=x, minlength=n_bins) / np.maximum(count, 1)
    return np.var(stack)


def periodic_template(
    trace,
    period_range=(3598.0, 3602.0),
    period_step=0.02,
    clip_mad=6.0,
    exclude=None,
):
    """Find the transient period by folding, then stack a robust template.

    Args:
        trace (Trace): Gap-free data, preferably high-passed (DC removed).
        period_range (tuple): Trial periods [s].
        period_step (float): Trial period spacing [s].
        clip_mad (float): Samples beyond this many MADs are clipped before
            folding so earthquakes do not dominate the stack.
        exclude (list of tuple): (start, end) UTCDateTime spans left out of
            the stack (earthquakes).

    Returns:
        PeriodicTemplate: Template at the best period.
    """
    fs = trace.stats.sampling_rate
    x = np.asarray(trace.data, float) - np.median(trace.data)
    t = trace.stats.starttime.timestamp + np.arange(len(x)) / fs
    keep = np.ones(len(x), bool)
    for start, end in exclude or []:
        keep &= ~((t >= start.timestamp) & (t <= end.timestamp))
    scale = 1.4826 * np.median(np.abs(x[keep])) or 1.0
    x = np.clip(x, -clip_mad * scale, clip_mad * scale)
    x, t = x[keep], t[keep]
    n_bins = int(round(np.mean(period_range) * fs))
    trial = np.arange(period_range[0], period_range[1] + period_step / 2, period_step)
    variance = np.array([_fold_variance(x, t - t[0], p, n_bins) for p in trial])
    best = trial[np.argmax(variance)]
    # Parabolic refinement around the peak.
    i = int(np.argmax(variance))
    if 0 < i < len(trial) - 1:
        a, b, c = variance[i - 1 : i + 2]
        denom = a - 2 * b + c
        if denom < 0:
            best += 0.5 * period_step * (a - c) / denom
    index = np.minimum((np.mod(t - t[0], best) / best * n_bins).astype(int), n_bins - 1)
    order = np.argsort(index, kind="stable")
    groups = np.split(x[order], np.cumsum(np.bincount(index, minlength=n_bins))[:-1])
    template = np.array([np.median(g) if len(g) else 0.0 for g in groups])
    template -= np.mean(template)
    return PeriodicTemplate(
        period=float(best),
        reference=UTCDateTime(t[0]),
        template=template,
        delta=best / n_bins,
        strength=float(variance.max() / np.median(variance)),
        n_cycles=int((t[-1] - t[0]) // best),
    )


def remove_periodic_transient(trace, template):
    """Subtract a periodic template from a trace (copy)."""
    out = trace.copy()
    fs = out.stats.sampling_rate
    times = out.stats.starttime.timestamp + np.arange(out.stats.npts) / fs
    out.data = np.asarray(out.data, float) - template.evaluate(times)
    return out


def rotation_coefficients(angle, azimuth):
    """Weights (cZ, c1, c2) of Z, 1 and 2 in TiSKitPy's rotated vertical.

    Reproduces ``SeisRotate.zrotate(angle, azimuth)`` for the vertical, which
    is linear in the three components.
    """
    dip_z = -90 + angle
    dip_n = angle * np.cos(np.radians(azimuth))
    dip_e = angle * np.sin(np.radians(azimuth))
    weights = []
    for unit in np.eye(3):
        z, _, _ = rotate2zne(
            np.array([unit[0]]), azimuth, dip_z,
            np.array([unit[1]]), 0, dip_n,
            np.array([unit[2]]), 90, dip_e,
        )  # fmt: skip
        weights.append(z[0])
    return tuple(weights)


def rotate_smoothly(stream, times, angles, azimuths, smooth_hours=24):
    """Rotate the vertical with smoothed, per-sample interpolated tilt.

    Hourly angle estimates are noisy (azimuth scatter of tens of degrees in
    Aminian et al. 2025, Fig. 4) and applying them as hourly steps puts
    offsets into displacement records. Here the tilt vector (dip_N, dip_E)
    is median filtered over ``smooth_hours`` and interpolated to each sample.

    Args:
        stream (Stream): Z, 1, 2 (and other) traces on a common grid.
        times (array): UTCDateTime centres of the hourly estimates.
        angles (array): Tilt angles [deg].
        azimuths (array): Tilt azimuths [deg].
        smooth_hours (float): Running-median length [h].

    Returns:
        Stream: Copy with the vertical replaced by the rotated vertical.
    """
    from scipy.ndimage import median_filter

    order = np.argsort([float(t) for t in times])
    centre = np.array([float(times[i]) for i in order])
    angle = np.asarray(angles, float)[order]
    azimuth = np.asarray(azimuths, float)[order]
    dip_n = angle * np.cos(np.radians(azimuth))
    dip_e = angle * np.sin(np.radians(azimuth))
    size = max(1, int(round(positive(smooth_hours, "smooth_hours"))))
    dip_n = median_filter(dip_n, size, mode="nearest")
    dip_e = median_filter(dip_e, size, mode="nearest")
    out = stream.copy()
    z = out.select(component="Z")[0]
    h1 = out.select(component="1")[0]
    h2 = out.select(component="2")[0]
    n = min(len(z), len(h1), len(h2))
    t = z.stats.starttime.timestamp + np.arange(n) / z.stats.sampling_rate
    tilt_n = np.interp(t, centre, dip_n)
    tilt_e = np.interp(t, centre, dip_e)
    tilt = np.hypot(tilt_n, tilt_e)
    az = np.degrees(np.arctan2(tilt_e, tilt_n))
    # Coefficients vary slowly; evaluate on a coarse grid and interpolate.
    knots = np.linspace(0, n - 1, max(2, n // 450))
    ki = knots.astype(int)
    coeff = np.array([rotation_coefficients(tilt[i], az[i]) for i in ki])
    c = [np.interp(np.arange(n), knots, coeff[:, j]) for j in range(3)]
    z.data = c[0] * z.data[:n] + c[1] * h1.data[:n] + c[2] * h2.data[:n]
    return out


def sta_lta_spans(
    trace, sta_s=60, lta_s=1800, on=3.5, off=1.5, pad_s=(300, 1800), band=(0.02, 0.2)
):
    """Time spans of local events detected by a recursive STA/LTA.

    Follows Withers et al. (1998) as used by Aminian et al. (2025): the
    trace is band-passed, triggered, and each trigger is padded to cover the
    P wave and later phases.

    Args:
        trace (Trace): Vertical channel.
        sta_s (float): Short-term window [s].
        lta_s (float): Long-term window [s].
        on (float): Trigger-on ratio.
        off (float): Trigger-off ratio.
        pad_s (tuple): Seconds added before and after each trigger.
        band (tuple): Band-pass corners [Hz] applied before triggering.

    Returns:
        list of tuple: (start, end) UTCDateTime spans.
    """
    tr = trace.copy()
    tr.detrend("linear")
    tr.filter("bandpass", freqmin=band[0], freqmax=band[1], zerophase=True)
    fs = tr.stats.sampling_rate
    cft = recursive_sta_lta(tr.data, int(sta_s * fs), int(lta_s * fs))
    cft[: int(lta_s * fs)] = 0  # LTA warm-up
    spans = []
    for a, b in trigger_onset(cft, on, off):
        start = tr.stats.starttime + a / fs - pad_s[0]
        end = tr.stats.starttime + b / fs + pad_s[1]
        if spans and start <= spans[-1][1]:
            spans[-1] = (spans[-1][0], max(end, spans[-1][1]))
        else:
            spans.append((start, end))
    return spans
