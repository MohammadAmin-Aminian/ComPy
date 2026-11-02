"""Tilt correction and coherent-noise removal for broadband OBS records.

``correct_tilt`` rotates the vertical onto the direction of least
horizontal noise window by window (tiskitpy ``CleanRotator``);
``remove_coherent_noise`` then subtracts what the horizontals predict on
the vertical (Crawford & Webb 2000). Aminian et al. (2025) use hourly
rotation followed by 12-24 h transfer functions.
"""

import warnings

import matplotlib.pyplot as plt
import numpy as np
import obspy
import tiskitpy as tiskit

from compy.numerics import positive
from compy.streams import split_stream


def correct_tilt(stream, time_window=1, *, plot=True, clean_horizontals=True):
    """Rotate and remove horizontal coherent noise in complete hour windows.

    Return (stream, azimuth_deg, angle_deg, variance_after_over_before).
    Processing failures raise with window context; input and seed metadata
    are preserved. A final incomplete window is omitted.

    With ``clean_horizontals=True`` (historical default) each window is also
    cleaned with a transfer function estimated from that window alone; a
    one-hour window gives only a few spectral averages, so coherence is biased
    high and part of the uncorrelated vertical signal is removed. Use
    ``clean_horizontals=False`` with :func:`remove_coherent_noise` to follow
    Aminian et al. (2025): hourly rotation, then 12–24 h transfer functions.
    """
    windows = split_stream(stream, positive(time_window, "time_window") * 3600)
    if not windows:
        raise ValueError("no complete rotation windows")
    azimuth, angles, variance_ratio, cleaned = [], [], [], []
    for i, window in enumerate(windows):
        try:
            before = np.var(window.select(component="Z")[0].data)
            if before <= 0:
                raise ValueError("vertical trace has zero variance")
            rotator = tiskit.CleanRotator(window, remove_eqs=False, verbose=False)
            rotated = rotator.apply(window, set_dtype=False)
            if clean_horizontals:
                cleaner = tiskit.DataCleaner(
                    rotated, ["*1", "*2"], n_to_reject=0, remove_eqs=False
                )
                rotated = cleaner.apply(rotated, set_dtype=False)
            after = np.var(rotated.select(component="Z")[0].data)
        except Exception as exc:
            raise RuntimeError(
                f"rotation failed in window {i} starting {window[0].stats.starttime}"
            ) from exc
        azimuth.append(rotator.azimuth % 360)
        angles.append(rotator.angle)
        variance_ratio.append(after / before)
        cleaned.extend(rotated)
    output = obspy.Stream(cleaned)
    output.merge(method=0)  # no interpolation of gaps or conflicting samples
    if plot:
        fig, axes = plt.subplots(2, 1, figsize=(12, 7), constrained_layout=True)
        for ax, values, label in zip(
            axes, (azimuth, angles), ("Azimuth [deg]", "Tilt [deg]")
        ):
            points = ax.scatter(
                np.arange(len(windows)) * time_window, values, c=variance_ratio
            )
            ax.set(xlabel="Hours since first complete window", ylabel=label)
            fig.colorbar(points, ax=ax, label="Variance after/before")
    return output, np.asarray(azimuth), np.asarray(angles), np.asarray(variance_ratio)


def remove_coherent_noise(
    stream,
    window_hours=24,
    *,
    noise_channels=("*1", "*2"),
    spectral_window_s=2048 / 1.25,
    avoid_spans=None,
    cover_all=False,
):
    """Remove horizontal-coherent noise from the vertical, window by window.

    Implements the transfer-function correction of Crawford & Webb (2000) as
    applied by Aminian et al. (2025, section 4.1.4): each 12–24 h window is
    cleaned with TiSKitPy's DataCleaner. The horizontals are removed one after
    the other, each from the vertical and from the remaining horizontals
    (conditioned spectra, Bendat & Piersol), so the result is the
    multiple-coherence correction and does not depend on their order; the
    paper's "dominant channel first" is therefore equivalent. Long windows give many spectral averages, avoiding
    the high-coherence bias of one-hour estimates. Pressure is never used as a
    noise channel, so the compliance signal itself is not removed.

    Args:
        stream (Stream): Rotated Z, 1, 2 and H traces on a common sample grid.
        window_hours (float): Cleaning window length [h].
        noise_channels (tuple): Horizontal channel patterns to remove.
        spectral_window_s (float): Welch segment length for the transfer
            functions [s].
        avoid_spans (TimeSpans or None): Earthquakes etc. excluded when
            estimating transfer functions (data inside are still cleaned).
        cover_all (bool): If True, use the smallest number of equal windows
            no longer than ``window_hours`` so the whole record is cleaned;
            otherwise a final incomplete window is omitted.

    Returns:
        Stream: Cleaned copy.
    """
    window_s = positive(window_hours, "window_hours") * 3600
    if cover_all:
        start = max(tr.stats.starttime for tr in stream)
        end = min(tr.stats.endtime for tr in stream)
        fs = stream[0].stats.sampling_rate
        n_samples = int(round((end - start) * fs)) + 1
        n_windows = max(1, int(np.ceil(n_samples / (window_s * fs))))
        window_s = (n_samples // n_windows) / fs
    windows = split_stream(stream, window_s)
    if not windows:
        raise ValueError("no complete cleaning windows")
    spans = []
    if avoid_spans is not None:
        spans = list(zip(avoid_spans.start_times, avoid_spans.end_times))
    cleaned = []
    for i, window in enumerate(windows):
        vertical = window.select(component="Z")
        if len(vertical) != 1:
            raise ValueError("exactly one vertical trace is required")
        noise = obspy.Stream(
            [tr for pattern in noise_channels for tr in window.select(channel=pattern)]
        )
        t0, t1 = window[0].stats.starttime, window[0].stats.endtime
        # Clip spans to this window: TiSKitPy rejects spans crossing its edges.
        inside = [(max(a, t0), min(b, t1)) for a, b in spans if a < t1 and b > t0]
        attempts = [{"avoid_spans": tiskit.TimeSpans(spans=inside)}] if inside else []
        attempts.append({})  # fall back to TiSKitPy's z-score window rejection
        z = None
        for extra in attempts:
            try:
                # DataCleaner corrects every channel it is given, so it only
                # sees Z and the noise channels; pressure passes through.
                cleaner = tiskit.DataCleaner(
                    vertical + noise,
                    list(noise_channels),
                    n_to_reject=0,
                    remove_eqs=False,
                    window_s=positive(spectral_window_s, "spectral_window_s"),
                    **extra,
                )
                out = cleaner.apply((vertical + noise).copy(), set_dtype=False)
                candidate = out.select(component="Z")[0]
                if not np.all(np.isfinite(candidate.data)):
                    raise ValueError("non-finite cleaned samples")
                z = candidate
                break
            except Exception as exc:
                error = exc
        if z is None:
            # Typically a window almost entirely inside event spans: keep the
            # vertical uncleaned so later quality gates decide on its use.
            warnings.warn(
                f"coherent-noise removal failed in window {i} starting {t0} "
                f"({error!r}); vertical left uncleaned"
            )
            z = vertical[0].copy()
        z.stats.location = vertical[0].stats.location
        z.stats.channel = vertical[0].stats.channel
        cleaned.append(z)
        cleaned.extend(tr for tr in window if tr.id != vertical[0].id)
    output = obspy.Stream(cleaned)
    output.merge(method=0)
    return output


def estimate_tilt_angles(stream, time_window=1):
    """Estimate angles and fractional variance reduction without applying cleaning."""
    windows = split_stream(stream, positive(time_window, "time_window") * 3600)
    if not windows:
        raise ValueError("no complete rotation windows")
    result = []
    for i, window in enumerate(windows):
        try:
            window.detrend("simple")
            rotator = tiskit.CleanRotator(window, remove_eqs=False)
        except Exception as exc:
            raise RuntimeError(f"angle estimation failed in window {i}") from exc
        result.append(
            (rotator.azimuth % 360, rotator.angle, rotator.variance_reduction)
        )
    return tuple(np.asarray(result).T)
