"""Teleseismic pressure-gain candidates; never modifies responses or masters."""

import argparse
import hashlib
import json
import re
from fractions import Fraction
from importlib.metadata import version
from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
from obspy import Stream, Trace, UTCDateTime, read, read_inventory
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel
from scipy.signal import resample_poly
from scipy.signal.windows import hann

from compy_pressure.archive import SampleGridError, discover_waveforms, read_window
from compy_pressure.common import complete_samples, write_csv
from compy_pressure.inventory import load_master_dataset
from compy_pressure.paired_spectra import (
    KEYS,
    checked_response,
    response_snapshot_path,
    select_pair,
    source_metadata,
)
from compy_pressure.ppsd import (
    atomic_json,
    is_pressure_unit,
    response_handling,
    select_epoch,
)
from compy_pressure.rayleigh_dispersion import DISPERSION_COLUMNS, load_dispersion

# Hash of the code as imported. Re-reading the file when each result is
# written would record an edit made while a long run was still executing
# the previously loaded code.
CODE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def real_array(value):
    if np.ma.isMaskedArray(value) and np.any(np.ma.getmaskarray(value)):
        raise ValueError("Masked values are not measurements")
    array = np.asarray(value)
    if array.dtype.kind not in "iuf" or not np.isfinite(array).all():
        raise ValueError("Expected finite real values")
    return array.astype(float)


def water_ratio(f, phase_velocity, depth, density, sound_speed):
    """Zha/Webb A5: P/upward-a, SI units; reject the first acoustic pole."""
    f, c = np.broadcast_arrays(real_array(f), real_array(phase_velocity))
    if any(not np.isfinite(v) or v <= 0 for v in (depth, density, sound_speed)):
        raise ValueError("Water parameters must be positive and finite")
    if (
        not np.isfinite(f).all()
        or not np.isfinite(c).all()
        or (f < 0).any()
        or (c <= 0).any()
    ):
        raise ValueError("Invalid frequency or phase velocity")
    q2 = sound_speed**-2 - c**-2
    x = 2 * np.pi * f * depth * np.sqrt(abs(q2))
    if np.any((q2 >= 0) & (x >= 0.9 * np.pi / 2)):
        raise ValueError("Band approaches first water-column resonance")
    result = np.ones_like(x)
    small = x < 1e-5
    result[small] = 1 + np.sign(q2[small]) * x[small] ** 2 / 3 + 2 * x[small] ** 4 / 15
    real = ~small & (q2 >= 0)
    imaginary = ~small & (q2 < 0)
    result[real] = np.tan(x[real]) / x[real]
    result[imaginary] = np.tanh(x[imaginary]) / x[imaginary]
    return density * depth * result


def fit_multiplier(observed, expected):
    """Fit g*observed=expected; g multiplies pressure, 1/g multiplies counts/Pa."""
    observed, expected = real_array(observed), real_array(expected)
    if observed.shape != expected.shape or observed.ndim != 1 or len(observed) < 3:
        raise ValueError("At least three corresponding frequency bins required")
    if (
        not (np.isfinite(observed).all() and np.isfinite(expected).all())
        or np.any(observed <= 0)
        or np.any(expected <= 0)
    ):
        raise ValueError("Fit ratios must be finite and positive")
    # Normalize first so SI scaling cannot overflow the objective.
    scale = np.max(expected)
    r, t = observed / scale, expected / scale
    g = float(np.dot(r, t) / np.dot(r, r))
    grid = np.linspace(0.1, 3.0, 2901)
    loss = np.mean((grid[:, None] * r - t) ** 2, axis=1)
    index = int(np.argmin(loss))
    if index in {0, len(grid) - 1}:
        raise ValueError(
            "Grid minimum lies on a boundary; expand justified search bounds"
        )
    if abs(g - grid[index]) > (grid[1] - grid[0]) / 2 + 1e-10:
        raise ValueError("Analytical and grid-search solutions disagree")
    return (
        {
            "pressure_multiplier": g,
            "sensitivity_counts_per_pa_multiplier": 1 / g,
            "grid_pressure_multiplier": float(grid[index]),
            "relative_rms": float(np.linalg.norm(g * r - t) / np.linalg.norm(t)),
        },
        grid,
        loss,
    )


def tilt_direction(covariance, max_angle_deg=10):
    """Quiet low-frequency minimum-variance axis; assumes orthogonal sensor axes."""
    cov = real_array(covariance)
    if cov.shape != (3, 3) or not np.isfinite(cov).all() or not np.allclose(cov, cov.T):
        raise ValueError("Invalid three-component covariance")
    _, vectors = np.linalg.eigh(cov)
    direction = vectors[:, 0]
    if direction[2] < 0:
        direction = -direction
    angle = float(np.degrees(np.arccos(np.clip(direction[2], -1, 1))))
    if angle > max_angle_deg:
        raise ValueError("Quiet-data optimum is not near vertical; tilt is unresolved")
    return direction, angle


def event_geometry(event, latitude, longitude, dispersion, duration=1800):
    """Use GROUP velocity for the packet, PHASE velocity only for the P/a model."""
    distance = gps2dist_azimuth(
        latitude, longitude, event["latitude"], event["longitude"]
    )[0]
    f = dispersion[:, 0]
    if f.min() > 0.04 or f.max() < 0.04:
        raise ValueError("Dispersion must bracket the 40 mHz packet")
    group = float(np.interp(0.04, f, dispersion[:, 2]))
    center = UTCDateTime(event["origin_utc"]) + distance / group
    return center - duration / 2, center + duration / 2, distance / 1000


def catalog_events(path):
    doc = json.loads(Path(path).read_text())
    if any(not re.fullmatch(r"[A-Za-z0-9_-]+", x["id"]) for x in doc["features"]):
        raise ValueError("Unsafe or malformed event identifier")
    return [
        {
            "id": x["id"],
            "origin_utc": str(UTCDateTime(x["properties"]["time"] / 1000)),
            "latitude": x["geometry"]["coordinates"][1],
            "longitude": x["geometry"]["coordinates"][0],
            "depth_km": x["geometry"]["coordinates"][2],
            "magnitude": x["properties"]["mag"],
            "magnitude_type": x["properties"]["magType"],
            "url": x["properties"]["url"],
        }
        for x in doc["features"]
    ]


def magnitude_schedule(
    events, summary, max_events, min_events, fallback, event_id, sampling="all"
):
    """Evaluate all eligible events, or explicitly select the old bounded pilot."""
    if sampling not in {"all", "bounded"}:
        raise ValueError("Unknown event sampling policy")
    bounded = sampling == "bounded"
    selected_origins = []
    levels = [(7.0, True), (6.5, False), (6.0, False)]
    if event_id:
        selected = [e for e in events if e["id"] == event_id]
        levels = [(selected[0]["magnitude"], False)] if selected else []
    elif not fallback:
        levels = levels[:1]
    for floor, strict in levels:
        if (
            bounded
            and sum(e["status"] == "candidate" for e in summary["events"]) >= min_events
        ):
            break
        eligible = [
            e
            for e in events
            if (e["magnitude"] > floor if strict else e["magnitude"] >= floor)
            and (not event_id or e["id"] == event_id)
        ]
        # An explicit event may be smaller than 7; no fallback/search is needed.
        if event_id:
            eligible = [e for e in events if e["id"] == event_id]
        record = {
            "minimum_magnitude": floor,
            "strictly_greater": strict,
            "eligible_catalog_events": len(eligible),
            "attempted_events": 0,
            "candidate_events_before": sum(
                e["status"] == "candidate" for e in summary["events"]
            ),
        }
        summary["magnitude_search"].append(record)
        print(f"Searching Mw {'>' if strict else '>='}{floor:g}", flush=True)
        seen = {e["id"] for e in summary["events"]}
        for event in sorted(eligible, key=lambda e: e["origin_utc"]):
            if event["id"] in seen:
                continue
            if (
                bounded
                and sum(e["status"] == "candidate" for e in summary["events"])
                >= min_events
            ):
                break
            if bounded and record["attempted_events"] >= max_events:
                break
            seen.add(event["id"])
            origin = UTCDateTime(event["origin_utc"])
            if any(abs(origin - previous) < 86400 for previous in selected_origins):
                summary["events"].append(
                    event
                    | {
                        "status": "excluded_cluster",
                        "reason": "Within 24 h of a previously attempted mainshock",
                        "selection_magnitude_floor": floor,
                    }
                )
                continue
            selected_origins.append(origin)
            record["attempted_events"] += 1
            yield event | {"selection_magnitude_floor": floor}
        record["candidate_events_after"] = sum(
            e["status"] == "candidate" for e in summary["events"]
        )
        record["outcome"] = (
            "target_reached"
            if bounded and record["candidate_events_after"] >= min_events
            else "attempt_limit_reached"
            if bounded and record["attempted_events"] >= max_events
            else "eligible_events_exhausted"
        )


def plot_paper_comparison(summary, dispersion, destination, show=False):
    """Figure 5 / thesis 3.13 layout, using only QC-accepted pilot events."""
    accepted = [e for e in summary["events"] if e["status"] == "candidate"]
    settings = summary["settings"]
    depth, density = settings["depth_m"], settings["density_kg_m3"]
    band = settings["band_hz"]
    curves = []
    for event in accepted:
        with np.load(event["files"]["spectra"]) as saved:
            curves.append(
                (
                    saved["paper_frequency_hz"],
                    saved["paper_coherence_squared"],
                    saved["paper_psd_ratio_pa"],
                )
            )
    fit = None
    fit_reason = None
    if curves:
        f = curves[0][0]
        for c in curves[1:]:
            f = np.intersect1d(f, c[0])
        indices = [np.searchsorted(c[0], f) for c in curves]
        median_coherence = np.median(
            [c[1][i] for c, i in zip(curves, indices, strict=True)], axis=0
        )
        median_ratio = np.median(
            [c[2][i] for c, i in zip(curves, indices, strict=True)], axis=0
        )
        fit_mask = (f >= band[0]) & (f <= band[1])
        fit_mask &= median_coherence >= settings["coherence_squared_threshold"]
        theory_fit = water_ratio(
            f[fit_mask],
            np.interp(f[fit_mask], dispersion[:, 0], dispersion[:, 1]),
            depth,
            density,
            settings["sound_speed_m_s"],
        )
        if fit_mask.sum() >= 3:
            try:
                fit, _, _ = fit_multiplier(median_ratio[fit_mask], theory_fit)
            except ValueError as exc:
                # This display diagnostic does not replace per-event station gains.
                fit_reason = str(exc)
    with plt.rc_context(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "font.size": 15,
            "axes.linewidth": 0.6,
            "grid.color": "#b0b0b0",
            "grid.linestyle": "--",
            "grid.linewidth": 0.5,
        }
    ):
        fig, axes = plt.subplots(2, 1, figsize=(10, 9))
        if curves:
            for frequency, coherence, ratio in curves:
                axes[0].semilogx(frequency, coherence, color="g", lw=0.6)
                axes[1].semilogx(
                    frequency, ratio / (density * depth), color="g", lw=0.6
                )
            axes[0].semilogx(
                f, median_coherence, color="blue", lw=2, label="Median of EQs"
            )
            axes[1].semilogx(
                f,
                median_ratio / (density * depth),
                color="blue",
                lw=2,
                label="Measured P/a ratio",
            )
        # Draw theory only where the supplied phase-velocity model is defined.
        theory_f = dispersion[:, 0]
        visible = (theory_f >= 0.005) & (theory_f <= 0.1)
        q2 = settings["sound_speed_m_s"] ** -2 - dispersion[:, 1] ** -2
        visible &= (q2 < 0) | (
            2 * np.pi * theory_f * depth * np.sqrt(abs(q2)) < 0.9 * np.pi / 2
        )
        theory = water_ratio(
            theory_f[visible],
            dispersion[visible, 1],
            depth,
            density,
            settings["sound_speed_m_s"],
        )
        axes[1].semilogx(
            theory_f[visible],
            theory / (density * depth),
            color="purple",
            lw=2,
            label="Theoretical P/a Ratio",
        )
        if fit:
            axes[1].semilogx(
                f,
                fit["pressure_multiplier"] * median_ratio / (density * depth),
                color="black",
                ls="--",
                lw=2,
                label="Corrected P/a ratio",
            )
        for i, ax in enumerate(axes):
            ax.set_xscale("log")
            ax.set_xlim(0.005, 0.1)
            ax.axvline(band[0], color="r", ls="--", lw=1.5, label="Frequency limits")
            ax.axvline(band[1], color="r", ls="--", lw=1.5)
            ax.grid(True)
            ax.minorticks_on()
            ax.text(
                0.01, 0.94, f"{'ab'[i]})", transform=ax.transAxes, fontweight="bold"
            )
        axes[0].set_ylim(-0.05, 1.05)
        axes[0].set_yticks(np.arange(0, 1.01, 0.2))
        axes[0].set_ylabel("Coherence")
        axes[0].set_title("High-Coherence Teleseismic Events")
        axes[0].legend(loc="lower right", fontsize=12)
        axes[1].set_ylim(-1, 10)
        axes[1].set_yticks(np.arange(0, 11, 2))
        axes[1].set_ylabel("Spectral ratio / ρh")
        axes[1].set_xlabel("Frequency [Hz]")
        station = summary["station"].split(".")[-1]
        gain_text = f"{fit['pressure_multiplier']:.2f}" if fit else "unresolved"
        axes[1].set_title(f"Spectral Ratio [Gain factor of {station} is {gain_text}]")
        axes[1].legend(loc="upper right", fontsize=12)
        fig.text(
            0.5,
            0.015,
            f"Provisional comparison • {len(accepted)} QC candidates • "
            "assumed pressure sign "
            f"{settings['pressure_sign']:+d} • "
            f"tilt method: {settings['tilt_method']} • "
            "no instrument correction applied",
            ha="center",
            fontsize=10,
        )
        fig.tight_layout(rect=(0, 0.035, 1, 1))
        fig.savefig(destination, dpi=180)
        if show:
            plt.show()
        plt.close(fig)
    summary["paper_comparison"] = {
        "fit_reason": fit_reason,
        "event_ids": [e["id"] for e in accepted],
        "coherence_kind": "magnitude_squared",
        "fit_kind": "gain of median event ratio; distinct from median event gains",
        "candidate_pressure_multiplier": fit["pressure_multiplier"] if fit else None,
        "frequency_axis_hz": [0.005, 0.1],
        "fit_band_hz": band,
        "theory_frequency_coverage_hz": [
            float(dispersion[0, 0]),
            float(dispersion[-1, 0]),
        ],
        "outside_fit_band": "display only; no gain calibration is claimed there",
        "reference_figures": "manuscript Figure 5; thesis Figure 3.13",
    }
    summary["files"]["paper_comparison"] = str(destination)


class SampleAlignmentError(ValueError):
    def __init__(self, rows, starts, rate):
        earliest = min(starts)
        self.timing = [
            {
                "seed_id": row.Seed_id,
                "actual_start_utc": str(start),
                "offset_from_earliest_seconds": float(start - earliest),
                "offset_samples": float(start - earliest) * rate,
            }
            for row, start in zip(rows, starts, strict=True)
        ]
        self.span_seconds = float(max(starts) - earliest)
        super().__init__(
            f"Unaligned sample grids: start spread {self.span_seconds:.9g} s"
        )


def native_resampling_samples(records, row, left, seconds, rate, tolerance_ppm):
    """Grid continuous finite native observations using their actual UTC times.

    Near-equal encoded rates require explicit ppm tolerance. Missing samples and
    conflicting overlaps are rejected before any interpolation or anti-aliasing.
    """
    paths = {
        r.path
        for r in records
        if r.seed_id == row.Seed_id and r.start < left + seconds and r.end >= left
    }
    traces = [
        t
        for p in sorted(paths)
        for t in read(
            str(p), format="MSEED", starttime=left, endtime=left + seconds + 2 / rate
        )
    ]
    times, values = [], []
    last_rate = rate
    for trace in sorted(traces, key=lambda t: t.stats.starttime):
        native_rate = float(trace.stats.sampling_rate)
        if abs(native_rate / rate - 1) * 1e6 > tolerance_ppm:
            raise ValueError("Native sample rate exceeds explicit ppm tolerance")
        data = real_array(trace.data)
        t = float(trace.stats.starttime - left) + np.arange(len(data)) / native_rate
        if not len(t):
            continue
        if times:
            last = times[-1][-1]
            if t[0] <= last + 1e-6:
                overlap = t <= last + 1e-6
                old_t, old_data = np.concatenate(times), np.concatenate(values)
                ix = np.searchsorted(old_t, t[overlap]).clip(0, len(old_t) - 1)
                ix = np.where(
                    (ix > 0)
                    & (abs(old_t[ix - 1] - t[overlap]) < abs(old_t[ix] - t[overlap])),
                    ix - 1,
                    ix,
                )
                if not (
                    np.all(abs(old_t[ix] - t[overlap]) <= 1e-6)
                    and np.array_equal(old_data[ix], data[overlap])
                ):
                    raise ValueError("Conflicting native overlap; resampling refused")
                t, data = t[~overlap], data[~overlap]
                if not len(t):
                    continue
            if t[0] - last > max(1 / native_rate, 1 / last_rate) + 1e-6:
                raise ValueError("Missing or gapped native window; resampling refused")
            if t[0] - last < min(1 / native_rate, 1 / last_rate) - 1e-6:
                raise ValueError("Fractional native-grid boundary; resampling refused")
        times.append(t)
        values.append(data)
        last_rate = native_rate
    if not times:
        raise ValueError("Missing native window")
    t, data = np.concatenate(times), np.concatenate(values)
    n = round(seconds * rate)
    grid = t[0] + np.arange(n) / rate
    if t[0] < -1 / rate - 1e-6 or t[0] > 1 / rate + 1e-6 or grid[-1] > t[-1]:
        raise ValueError("Missing native window coverage")
    if np.ptp(data) == 0:
        raise ValueError("Flat native window")
    return np.interp(grid, t, data), left + t[0]


def load_resampled(
    records, rows, channels, start, seconds, target_rate, native_rate_tolerance_ppm=0
):
    """Anti-alias working copies onto UTC grids; never fill gaps or shift clocks.

    Counts are resampled before frequency-domain response correction. Responses
    retain their native rates/epochs; no resampled inventory is fabricated.
    Padding protects the retained window from FIR/interpolation edge effects.
    """
    rates = [float(c.sample_rate) for c in channels]
    if not np.isfinite(target_rate) or not 0 < target_rate <= min(rates):
        raise ValueError("Event resample rate must be positive and not upsample")
    n = round(seconds * target_rate)
    if not np.isclose(n / target_rate, seconds, rtol=0, atol=1e-8):
        raise ValueError("Duration must contain an integer number of samples")
    padding = 32 / target_rate
    left, right = start - padding, start + seconds + padding
    arrays = []
    for row, channel, rate in zip(rows, channels, rates, strict=True):
        factor = Fraction(str(target_rate)) / Fraction(str(rate))
        if max(factor.numerator, factor.denominator) > 10000:
            raise ValueError("Unsupported exact resampling ratio")
        if (channel.start_date and left < channel.start_date) or (
            channel.end_date and right + 2 / rate > channel.end_date
        ):
            raise ValueError("Padded window crosses response epoch")
        native_n = round((seconds + 2 * padding) * rate)
        if not np.isclose(native_n / rate, seconds + 2 * padding, atol=1e-8):
            raise ValueError("Padded duration must contain native integer samples")
        array, actual = native_resampling_samples(
            records,
            row,
            left,
            seconds + 2 * padding,
            rate,
            native_rate_tolerance_ppm,
        )
        if array is None:
            raise ValueError("Missing, gapped, nonfinite or flat padded window")
        sampled = resample_poly(
            array - array.mean(),
            factor.numerator,
            factor.denominator,
            window=("kaiser", 12.0),
        )
        trace = Trace(sampled, header=dict(starttime=actual, sampling_rate=target_rate))
        # Native P/Z grids may start at distinct sub-sample UTC offsets. Evaluate
        # the continuous, filtered data on one grid; this is not clock correction.
        trace.interpolate(target_rate, method="lanczos", starttime=start, npts=n, a=12)
        arrays.append(real_array(trace.data))
    return np.array(arrays), start, target_rate


def load_aligned(
    records,
    rows,
    channels,
    start,
    seconds,
    event_resample_rate=None,
    native_rate_tolerance_ppm=0,
):
    if event_resample_rate is not None:
        return load_resampled(
            records,
            rows,
            channels,
            start,
            seconds,
            event_resample_rate,
            native_rate_tolerance_ppm,
        )
    rate = float(channels[0].sample_rate)
    if any(c.sample_rate != rate for c in channels):
        raise ValueError("Unequal sample rates; no implicit resampling")
    n = round(seconds * rate)
    if not np.isclose(n / rate, seconds, atol=1e-8, rtol=0):
        raise ValueError("Duration must contain an integer number of samples")
    arrays, starts = [], []
    for row, channel in zip(rows, channels, strict=True):
        if (channel.start_date and start < channel.start_date) or (
            channel.end_date and start + seconds > channel.end_date
        ):
            raise ValueError("Requested window crosses response epoch")
        rr = [r for r in records if r.seed_id == row.Seed_id]
        array, actual = complete_samples(
            read_window(rr, start, start + seconds + 2 / rate, strict_sample_grid=True),
            start,
            rate,
            n,
        )
        if array is None:
            raise ValueError("Missing, gapped, nonfinite or flat window")
        if channel.end_date and actual + seconds > channel.end_date:
            raise ValueError("Aligned window crosses response epoch")
        arrays.append(array)
        starts.append(actual)
    if max(starts) - min(starts) > 1e-6:
        raise SampleAlignmentError(rows, starts, rate)
    return np.array(arrays), starts[0], rate


def corrected_fft(arrays, rows, channels, rate, frequency_min, frequency_max):
    taper = hann(arrays.shape[1], sym=False)
    f = np.fft.rfftfreq(arrays.shape[1], 1 / rate)
    mask = (f >= frequency_min) & (f <= frequency_max)
    fft = np.fft.rfft((arrays - arrays.mean(axis=1, keepdims=True)) * taper, axis=1)[
        :, mask
    ]
    for i, (row, channel) in enumerate(zip(rows, channels, strict=True)):
        output = "DEF" if i == 0 else "VEL"
        response = checked_response(row, channel, f[mask], output)
        if not np.all(np.isfinite(response) & (abs(response) > 0)):
            raise ValueError("Invalid complex response in analysis band")
        fft[i] /= response
        if i == len(channels) - 1 and channel.dip > 0:
            fft[i] *= -1
    return f[mask], fft, 2 / (rate * np.sum(taper**2))


def quiet_axis(arrays, rows, channels, rate):
    """Experimental local axis; needs comparative field testing, not the default."""
    # Same physical velocity units on all three seismic components.
    _, coefficients, scale = corrected_fft(arrays, rows, channels, rate, 0.001, 0.01)
    v = coefficients[1:]
    covariance = np.real(v @ v.conj().T) * scale
    direction, angle = tilt_direction(covariance)
    reduction = float(direction @ covariance @ direction / covariance[2, 2])
    return direction, angle, reduction


def require_clean_rotator():
    try:
        from tiskitpy import CleanRotator
    except ImportError as exc:
        raise RuntimeError(
            "Install TiSKitPy with: python -m pip install -e '.[rayleigh]'"
        ) from exc
    return CleanRotator


def physical_quiet_stream(arrays, rows, channels, rate, start):
    """Response-correct temporary training traces once, then anti-alias to ≤1 Hz."""
    stream = Stream()
    target_rate = min(1.0, rate)
    factor = Fraction(str(target_rate / rate)).limit_denominator(10000)
    if not np.isclose(
        rate * factor.numerator / factor.denominator, target_rate, rtol=0, atol=1e-10
    ):
        raise ValueError("Quiet resampling ratio does not exactly represent the rate")
    for i in range(1, 4):
        row, channel = rows[i], channels[i]
        # The time-domain deconvolution uses the same snapshot already checked here.
        checked_response(row, channel, np.array([0.001, 0.005, 0.01]), "VEL")
        net, sta, loc, code = row.Seed_id.split(".")
        trace = Trace(
            arrays[i].copy(),
            header=dict(
                network=net,
                station=sta,
                location=loc,
                channel=code,
                starttime=start,
                sampling_rate=rate,
            ),
        )
        inventory = read_inventory(
            str(Path(row.ppsd_npz_file).parent / "inventory.xml")
        )
        trace.remove_response(
            inventory=inventory,
            output="VEL",
            water_level=None,
            pre_filt=(0.0003, 0.0005, 0.03, 0.05),
            taper=True,
            taper_fraction=0.05,
            zero_mean=True,
        )
        if i == 3 and channel.dip > 0:
            trace.data *= -1
        if target_rate < rate:
            trace.data = resample_poly(
                trace.data, factor.numerator, factor.denominator, window=("kaiser", 5.0)
            )
            trace.stats.sampling_rate = target_rate
        if not np.isfinite(trace.data).all():
            raise ValueError("Nonfinite physical quiet-training samples")
        stream += trace
    return stream


def clean_rotator_axis(stream):
    """Fit the real TiSKitPy class; recover its linear Z weights using apply()."""
    clean_rotator = require_clean_rotator()
    rotator = clean_rotator(
        stream,
        remove_eqs=False,
        save_eq_file=False,
        plot=False,
        verbose=False,
        filt_band=(0.001, 0.01),
        H_over_Z=1.0,
        N_E_comps=("1", "2"),
    )
    basis = stream.copy()
    for i, trace in enumerate(basis):
        trace.data = np.zeros(32, dtype=float)
        trace.data[i] = 1
    # tfs() is documented as unvalidated. Use the public apply() itself instead.
    rotated_basis = rotator.apply(basis, horiz_too=False, rot_limit=10, set_dtype=False)
    vertical = rotated_basis.select(component="Z")
    if len(vertical) != 1 or len(vertical[0]) != 32:
        raise ValueError("CleanRotator changed the synchronized basis shape")
    direction = np.asarray(vertical[0].data[:3], dtype=float)
    if not np.isfinite(direction).all() or not np.isclose(np.linalg.norm(direction), 1):
        raise ValueError("Invalid CleanRotator projection")
    return (
        direction,
        float(rotator.angle),
        1 - float(rotator.variance_reduction),
        {
            "method": "tiskitpy.CleanRotator",
            "package_version": version("tiskitpy"),
            "angle_deg": float(rotator.angle),
            "azimuth_sensor_basis_deg": float(rotator.azimuth),
            "variance_reduction": float(rotator.variance_reduction),
            "constructor_filter_band_hz": [0.001, 0.01],
            "optimizer_filter_band_hz": [0.001, 0.005],
            "earthquake_screen": "preceding interval screened by supplied catalogue",
            "implicit_catalog_download": False,
            "horizontal_rotation": False,
            "reference": "https://tiskitpy.readthedocs.io/latest/classes/clean_rotator.html",
        },
    )


def estimate_tilt(arrays, rows, channels, rate, start, method="tiskitpy"):
    if method == "local_eigen":
        direction, angle, ratio = quiet_axis(arrays, rows, channels, rate)
        return (
            direction,
            angle,
            ratio,
            {
                "method": "local_eigen",
                "validation_status": "needs_to_be_tested",
                "note": "Synthetic checks exist; comparative field testing is pending",
                "reference": "docs/tilt_methods.md#experimental-local-method",
            },
        )
    if method != "tiskitpy":
        raise ValueError("Unknown tilt method")
    stream = physical_quiet_stream(arrays, rows, channels, rate, start)
    direction, angle, ratio, metadata = clean_rotator_axis(stream)
    metadata.update(
        training_rate_hz=float(stream[0].stats.sampling_rate),
        response_output="VEL",
        response_prefilter_hz=[0.0003, 0.0005, 0.03, 0.05],
        anti_alias_resampling="scipy.signal.resample_poly; Kaiser beta=5",
    )
    return direction, angle, ratio, metadata


def event_spectra(
    arrays, rows, channels, rate, direction, segment_seconds, band, pressure_sign=1
):
    n = round(segment_seconds * rate)
    if not np.isclose(n / rate, segment_seconds, rtol=0, atol=1e-8):
        raise ValueError("Invalid segment duration")
    count = arrays.shape[1] // n
    if count < 8:
        raise ValueError("At least eight nonoverlapping segments required")
    pp, aa, ap, raw_aa, raw_ap = [], [], [], [], []
    for i in range(count):
        f, coefficients, scale = corrected_fft(
            arrays[:, i * n : (i + 1) * n], rows, channels, rate, *band
        )
        pressure = coefficients[0]
        velocity = direction @ coefficients[1:]
        acc = 2j * np.pi * f * velocity
        raw_acc = 2j * np.pi * f * coefficients[-1]
        pp.append(abs(pressure) ** 2 * scale)
        aa.append(abs(acc) ** 2 * scale)
        ap.append(pressure * pressure_sign * acc.conj() * scale)
        raw_aa.append(abs(raw_acc) ** 2 * scale)
        raw_ap.append(pressure * pressure_sign * raw_acc.conj() * scale)
    pp, aa, ap = np.mean(pp, axis=0), np.mean(aa, axis=0), np.mean(ap, axis=0)
    raw_aa, raw_ap = np.mean(raw_aa, axis=0), np.mean(raw_ap, axis=0)
    if not np.all(
        np.isfinite(pp) & np.isfinite(aa) & np.isfinite(ap) & (aa > 0) & (pp > 0)
    ):
        raise ValueError("Invalid event spectra")
    return {
        "frequency_hz": f,
        "pressure_psd": pp,
        "acceleration_psd": aa,
        "cross_pa": ap,
        "h1_pa": ap / aa,
        "psd_ratio_pa": np.sqrt(pp / aa),
        "coherence_squared": np.clip(abs(ap) ** 2 / (pp * aa), 0, 1),
        "phase_deg": np.angle(ap, deg=True),
        "nominal_phase_deg": np.angle(ap * pressure_sign, deg=True),
        "pressure_sign": pressure_sign,
        "raw_h1_pa": raw_ap / raw_aa,
        "raw_psd_ratio_pa": np.sqrt(pp / raw_aa),
        "raw_coherence_squared": np.clip(abs(raw_ap) ** 2 / (pp * raw_aa), 0, 1),
        "segment_count": count,
        "trailing_samples": arrays.shape[1] % n,
    }


def plot_arrival_review(
    records,
    rows,
    channels,
    origin,
    start,
    end,
    band,
    group_arrival_range,
    p_utc,
    s_utc,
    other_windows,
    pressure_sign,
    destination,
):
    """Wider, independently timed P/Z view; never interpolates or chooses a pick."""
    duration = np.ceil(float(end - origin + 900) / 60) * 60
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    for ax, index in zip(axes, (0, len(rows) - 1), strict=True):
        row, channel = rows[index], channels[index]
        rate = float(channel.sample_rate)
        if (channel.start_date and origin < channel.start_date) or (
            channel.end_date and origin + duration > channel.end_date
        ):
            raise ValueError("Overview crosses response epoch")
        rr = [r for r in records if r.seed_id == row.Seed_id]
        n = round(duration * rate)
        data, actual = complete_samples(
            read_window(rr, origin, origin + duration + 2 / rate), origin, rate, n
        )
        if data is None:
            plt.close(fig)
            raise ValueError("Overview has a gap, flat data or insufficient coverage")
        f = np.fft.rfftfreq(n, 1 / rate)
        mask = (f >= band[0]) & (f <= band[1])
        response = checked_response(
            row, channel, f[mask], "DEF" if index == 0 else "VEL"
        )
        if not np.all(np.isfinite(response) & (abs(response) > 0)):
            plt.close(fig)
            raise ValueError("Invalid overview response")
        fourier = np.zeros(len(f), complex)
        fourier[mask] = (
            np.fft.rfft((data - data.mean()) * hann(n, sym=False))[mask] / response
        )
        if index == 0:
            fourier *= pressure_sign
        else:
            fourier[mask] *= 2j * np.pi * f[mask] * (-1 if channel.dip > 0 else 1)
        wave = np.fft.irfft(fourier, n)
        step = max(1, round(rate))
        elapsed = float(actual - origin) + np.arange(0, n, step) / rate
        ax.plot(elapsed, wave[::step], lw=0.8)
        ax.axvspan(
            float(start - origin),
            float(end - origin),
            alpha=0.16,
            color="purple",
            label="Selected window",
        )
        ax.axvspan(
            *group_arrival_range,
            alpha=0.2,
            color="green",
            label="Model group-arrival range",
        )
        for stamp, label, color in [
            (p_utc, "Predicted P", "red"),
            (s_utc, "Predicted S", "orange"),
        ]:
            if stamp:
                ax.axvline(
                    float(UTCDateTime(stamp) - origin),
                    color=color,
                    ls="--",
                    label=label,
                )
        for i, (left, right) in enumerate(other_windows):
            if left < origin + duration and right > origin:
                ax.axvspan(
                    float(left - origin),
                    float(right - origin),
                    color="gray",
                    alpha=0.15,
                    label="Other predicted packet" if i == 0 else None,
                )
        ax.set_ylabel(
            "Pressure (Pa; sign assumed)" if index == 0 else "Upward a (m/s²)"
        )
        ax.grid(alpha=0.2)
    axes[0].legend(fontsize=8, loc="upper left")
    axes[1].set_xlabel("Seconds after event origin")
    fig.suptitle(
        "Arrival review — inspect packet and overlaps; no automatic pick\n"
        "Band-limited, Hann-tapered display; reference model is provisional"
    )
    fig.tight_layout()
    fig.savefig(destination, dpi=140)
    plt.close(fig)


def select_calibration_channels(
    master,
    network,
    station,
    pressure_epoch=None,
    tilt_method="tiskitpy",
    event_resample_rate=None,
):
    """Bind P/Z epochs, plus both horizontal epochs when rotation is requested."""
    if tilt_method not in {"none", "tiskitpy", "local_eigen"}:
        raise ValueError("Unknown tilt method")
    pressure = master[
        master.Seed_id.str.startswith(f"{network}.{station}.")
        & master.input_units.map(is_pressure_unit)
    ].drop_duplicates(KEYS)
    if pressure_epoch is not None:
        for key in KEYS:
            pressure = pressure[pressure[key].eq(str(pressure_epoch[key]))]
    if len(pressure) != 1:
        raise ValueError("Select one pressure epoch")
    prow = pressure.iloc[0]
    if event_resample_rate is not None and tilt_method != "none":
        raise ValueError("Event resampling currently requires tilt_method=none")
    pchannel, zrow, zchannel = select_pair(
        prow, master, allow_unequal_rates=event_resample_rate is not None
    )
    prefix = prow.Seed_id.rsplit(".", 1)[0]
    # Horizontals share the paired vertical's band/instrument codes (BHZ -> BH1/2).
    band = zrow.Seed_id.rsplit(".", 1)[1][:2]
    hrows, hchannels = [], []
    for suffix in () if tilt_method == "none" else (f".{band}1", f".{band}2"):
        matches = []
        for _, row in (
            master[master.Seed_id.eq(prefix + suffix)].drop_duplicates(KEYS).iterrows()
        ):
            _, channel = select_epoch(row)
            starts = [c.start_date for c in (pchannel, channel) if c.start_date]
            ends = [c.end_date for c in (pchannel, channel) if c.end_date]
            if not starts or not ends or max(starts) < min(ends):
                matches.append((row, channel))
        if len(matches) != 1:
            raise ValueError(
                f"Expected one overlapping {suffix} epoch, found {len(matches)}"
            )
        row, channel = matches[0]
        hrows.append(row)
        hchannels.append(channel)
    rows, channels = [prow, *hrows, zrow], [pchannel, *hchannels, zchannel]
    for index, (row, channel) in enumerate(zip(rows, channels, strict=True)):
        if index > 0 and response_handling(channel)[1] != "m^2/s^4/Hz":
            raise ValueError(
                "Horizontal/vertical response must have seismic input units"
            )
        if 0 < index < len(channels) - 1 and (
            channel.dip is None or abs(channel.dip) > 0.1
        ):
            raise ValueError("Pilot requires horizontal component orientations")
        standalone_snapshot = row.get("rayleigh_response_snapshot_file", "")
        if not channel.response or (
            not standalone_snapshot and row.ppsd_status not in {"complete", "partial"}
        ):
            raise ValueError("Every pilot channel needs a saved response snapshot")
        if channel.sample_rate != pchannel.sample_rate and event_resample_rate is None:
            raise ValueError("Unequal sample rates; no implicit event resampling")
        if not response_snapshot_path(row).is_file():
            raise ValueError("Missing saved response snapshot")
    return rows, channels


def review_assumptions(prow, tilt_method):
    """Working implementation; these flags identify unvalidated field assumptions."""
    flags = {
        "pressure_polarity_assumed": "Single supplied sign; verify hardware convention",
        "vertical_absolute_response_unverified": "Gain depends on Z calibration",
        "depth_from_inventory": "Verify bathymetry or explicitly supplied depth",
        "earth_model_needs_station_validation": "Validate the selected model locally",
        "water_properties_assumed": "Homogeneous water density and sound speed",
        "sensor_axes_assumed_orthogonal": "Verify sensor-axis geometry",
        "quiet_screen_and_coda_heuristic": "Training/coda screen needs testing",
        "tilt_stationarity_needs_testing": "Tilt may change during events",
        "rayleigh_window_mode_needs_review": "Verify isolated fundamental-mode packets",
        "qc_rules_need_field_validation": "QC/window/count rules are pilot choices",
        "magnitude_fallback_needs_testing": "Check weaker/mixed wavefields",
        "scalar_gain_frequency_independence": "No time-constant calibration",
        "uncertainty_incomplete": "Scatter omits systematic uncertainties",
    }
    sensor = str(prow.get("sensor_description", ""))
    if "DPG" not in sensor.upper() and "DIFFEREN" not in sensor.upper():
        flags["pressure_hardware_needs_confirmation"] = (
            f"Declared label {sensor!r} does not establish DPG hardware"
        )
    if tilt_method == "local_eigen":
        flags["local_tilt_method_needs_to_be_tested"] = "Experimental opt-in method"
    if tilt_method == "none":
        for key in (
            "sensor_axes_assumed_orthogonal",
            "quiet_screen_and_coda_heuristic",
            "tilt_stationarity_needs_testing",
        ):
            flags.pop(key)
    return [
        {"id": key, "status": "needs_review_or_testing", "note": note}
        for key, note in flags.items()
    ]


def run_calibration(
    *,
    catalog,
    dispersion_file=None,
    earth_model_file=None,
    master_dir="exports",
    data_dir=None,
    network="YV",
    station="RR52",
    band=(0.02, 0.06),
    depth=None,
    density=1028,
    sound_speed=1500,
    threshold=0.95,
    max_events=8,
    min_events=3,
    magnitude_fallback=True,
    event_sampling="all",
    event_id=None,
    show=False,
    tilt_method="tiskitpy",
    event_resample_rate=None,
    pressure_sign=None,
    native_rate_tolerance_ppm=0,
    window_overrides=None,
    pressure_epoch=None,
    master_dataset=None,
    waveform_records=None,
    _resolved_dispersion=None,
):
    """Save auditable candidate comparisons; never update inventory or masters."""
    if (
        not 0 < band[0] < band[1]
        or not 0 < threshold < 1
        or max_events < 1
        or min_events < 3
    ):
        raise ValueError("Invalid band, coherence threshold or event count")
    if event_sampling not in {"all", "bounded"}:
        raise ValueError("Unknown event sampling policy")
    if (
        not np.isfinite(native_rate_tolerance_ppm)
        or not 0 <= native_rate_tolerance_ppm <= 0.1
    ):
        raise ValueError("Native rate tolerance must be between 0 and 0.1 ppm")
    if native_rate_tolerance_ppm and event_resample_rate is None:
        raise ValueError("Native rate tolerance requires explicit event resampling")
    sign_origin = (
        "RR52_provisional_default" if pressure_sign is None else "explicit_argument"
    )
    if pressure_sign is None:
        if (network, station) != ("YV", "RR52"):
            raise ValueError("Supply an explicit pressure sign outside the RR52 pilot")
        pressure_sign = -1
    if pressure_sign not in (-1, 1):
        raise ValueError("Pressure sign must be -1 or +1")
    if tilt_method == "tiskitpy":
        require_clean_rotator()  # No silent fallback to the experimental method.
    elif tilt_method not in {"none", "local_eigen"}:
        raise ValueError("Unknown tilt method")
    dispersion, dispersion_model = _resolved_dispersion or load_dispersion(
        dispersion_file=dispersion_file, earth_model_file=earth_model_file, band=band
    )
    dispersion_model = dict(dispersion_model)
    model_input = Path(dispersion_model["source_file"])
    if (
        hashlib.sha256(model_input.read_bytes()).hexdigest()
        != dispersion_model["source_sha256"]
    ):
        raise ValueError("Model input changed during dispersion calculation")
    master = (
        load_master_dataset(master_dir) if master_dataset is None else master_dataset
    )
    rows, channels = select_calibration_channels(
        master, network, station, pressure_epoch, tilt_method, event_resample_rate
    )
    prow = rows[0]
    if depth is None:
        depth = float(prow.Profondeur)
    analysis_rate = min(float(c.sample_rate) for c in channels)
    if event_resample_rate is not None:
        if (
            not np.isfinite(event_resample_rate)
            or not 0 < event_resample_rate <= analysis_rate
        ):
            raise ValueError("Event resample rate must be positive and not upsample")
        analysis_rate = event_resample_rate
        if band[1] > analysis_rate / 10:
            raise ValueError(
                "Resampled fit band must be below one tenth of working rate"
            )
    if band[1] >= analysis_rate / 2:
        raise ValueError("Calibration band must be strictly below Nyquist")
    water_ratio([band[1]], dispersion[:, 1].max(), depth, density, sound_speed)
    records = waveform_records
    if records is None:
        records = discover_waveforms(
            data_dir or Path.home() / "Desktop/Data_Bruit_FM/Data",
            {r.Seed_id for r in rows},
        )
    output = (
        response_snapshot_path(prow).parent
        / "rayleigh_calibration"
        / f"run_{uuid4().hex[:12]}"
    )
    output.mkdir(parents=True, exist_ok=False)
    generated_dispersion = None
    if earth_model_file:
        generated_dispersion = output / "model_dispersion.csv"
        np.savetxt(
            generated_dispersion,
            dispersion,
            delimiter=",",
            header=DISPERSION_COLUMNS,
            comments="",
        )
        dispersion_model["generated_dispersion_sha256"] = hashlib.sha256(
            generated_dispersion.read_bytes()
        ).hexdigest()
    override_doc = (
        json.loads(Path(window_overrides).read_text()) if window_overrides else {}
    )
    if not isinstance(override_doc, dict):
        raise ValueError("Window overrides must be keyed by event ID")
    all_events = catalog_events(catalog)
    candidates = sorted(
        [
            e
            for e in all_events
            if e["magnitude"] >= 6
            and e["depth_km"] < 150
            and str(e["magnitude_type"]).lower().startswith("mw")
            and all(
                (not c.start_date or UTCDateTime(e["origin_utc"]) >= c.start_date)
                and (not c.end_date or UTCDateTime(e["origin_utc"]) < c.end_date)
                for c in channels
            )
        ],
        key=lambda e: e["origin_utc"],
    )
    if event_id:
        candidates = [e for e in candidates if e["id"] == event_id]
    summary = {
        "status": "unresolved",
        "implementation_status": "working_method_with_provisional_assumptions",
        "review_flags": review_assumptions(prow, tilt_method),
        "dispersion_model": dispersion_model,
        "station": f"{network}.{station}",
        "epoch_keys": [{k: str(r[k]) for k in KEYS} for r in rows],
        "settings": {
            "band_hz": list(band),
            "depth_m": depth,
            "density_kg_m3": density,
            "sound_speed_m_s": sound_speed,
            "coherence_squared_threshold": threshold,
            "phase_limit_deg": 20,
            "window_seconds": 1800,
            "segment_seconds": 200,
            "tilt_method": tilt_method,
            "tilt_correction_applied": tilt_method != "none",
            "event_resample_rate_hz": event_resample_rate,
            "native_rate_tolerance_ppm": native_rate_tolerance_ppm,
            "event_resampling": {
                "method": "polyphase Kaiser beta=12; UTC Lanczos a=12",
                "padding_seconds": 32 / event_resample_rate,
                "native_rates_hz": [float(c.sample_rate) for c in channels],
                "response_correction": (
                    "native complex response after resampling counts"
                ),
                "clock_correction_applied": False,
                "gaps_filled": False,
                "native_gridding": "linear at native rate using actual UTC timestamps",
                "native_rate_tolerance_ppm": native_rate_tolerance_ppm,
            }
            if event_resample_rate is not None
            else None,
            "pressure_sign": pressure_sign,
            "minimum_candidate_events": min_events,
            "event_sampling": event_sampling,
            "max_attempts_per_magnitude_level": max_events
            if event_sampling == "bounded"
            else None,
            "magnitude_fallback": bool(magnitude_fallback and not event_id),
            # Pilot assumptions, not published instrument acceptance limits.
            "source_depth_limit_km": 150,
            "mainshock_separation_hours": 24,
            "packet_center_frequency_hz": 0.04,
            "quiet_window_seconds": 7200 if tilt_method != "none" else None,
            "quiet_search_offsets_hours": [6, 12, 24, 48]
            if tilt_method != "none"
            else [],
            "quiet_coda_rule": "max(3600, (M-5.5)*43200) seconds"
            if tilt_method != "none"
            else None,
            "tilt_angle_limit_deg": 10 if tilt_method != "none" else None,
            "sample_start_tolerance_seconds": 1e-6,
            "gain_grid_bounds": [0.1, 3.0],
            "gain_grid_step": 0.001,
        },
        "pressure_sign_assumption": {
            "multiplier": pressure_sign,
            "origin": sign_origin,
            "status": "assumed_needs_confirmation",
            "reason": (
                "User recalls opposite signs; old RR52 plot inverted pressure"
                if sign_origin == "RR52_provisional_default"
                else "Caller-specified sign; hardware convention not confirmed here"
            ),
            "source_responses_changed": False,
        },
        "inputs": {
            str(p): hashlib.sha256(Path(p).read_bytes()).hexdigest()
            for p in [
                catalog,
                model_input,
                *([window_overrides] if window_overrides else []),
            ]
        },
        "source_inventory_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in {Path(r.inventory_file) for r in rows}
        },
        "response_snapshot_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in {response_snapshot_path(r) for r in rows}
        },
        "code_sha256": CODE_SHA256,
        "limitations": [
            "Candidate multiplier only; no response/master changes",
            "Selected Earth model and homogeneous water properties are assumptions",
            "Unrotated vertical component; no quiet training or tilt correction"
            if tilt_method == "none"
            else (
                "Quiet rotation assumes orthogonal axes; "
                "azimuth metadata are incomplete"
            ),
            "Vertical absolute calibration, clock and bathymetry are unverified",
            "No time-constant calibration or extrapolation outside the fit band",
            "Manual event isolation and uncertainty validation remain required",
        ],
        "events": [],
        "files": {"dispersion": str(generated_dispersion)}
        if generated_dispersion
        else {},
        "source_files": [],
        "magnitude_search": [],
    }
    if event_resample_rate is not None:
        summary["review_flags"].extend(
            [
                {
                    "id": "native_resampling_needs_field_validation",
                    "status": "needs_review_or_testing",
                    "note": "Verify field timing and filter sensitivity",
                },
                {
                    "id": "clock_alignment_needs_review",
                    "status": "needs_review_or_testing",
                    "note": "UTC gridding does not correct logger clock errors",
                },
            ]
        )
    for event in magnitude_schedule(
        candidates,
        summary,
        max_events,
        min_events,
        magnitude_fallback,
        event_id,
        sampling=event_sampling,
    ):
        origin = UTCDateTime(event["origin_utc"])
        event_sources = []
        item = event | {"status": "unresolved"}
        summary["events"].append(item)
        try:
            start, end, distance = event_geometry(
                event, float(prow.Lat), float(prow.Lon), dispersion
            )
            override = override_doc.get(event["id"], {})
            if not isinstance(override, dict):
                raise ValueError("Each window override must be an object")
            decision = override.get("decision", "pending")
            if decision not in {"pending", "use", "exclude"}:
                raise ValueError("Unknown manual window decision")
            if decision == "exclude":
                item.update(
                    status="excluded_manual",
                    reason=override.get("note", "Manual exclusion"),
                )
                continue
            if decision == "use":
                if not str(override.get("note", "")).strip():
                    raise ValueError("A reviewed window needs an explanation")
                if not {"start_utc", "end_utc"}.issubset(override):
                    raise ValueError("Reviewed window needs both UTC bounds")
                start, end = (
                    UTCDateTime(override["start_utc"]),
                    UTCDateTime(override["end_utc"]),
                )
                if end - start < 1600:
                    raise ValueError("Reviewed window needs eight 200 s segments")
            if start < origin or end - origin > 43200 or end - start > 7200:
                raise ValueError(
                    "Reviewed window must follow origin and span at most 2 h"
                )
            window_seconds = float(end - start)
            item.update(
                start_utc=str(start),
                end_utc=str(end),
                distance_km=distance,
                window_review="reviewed_override"
                if decision == "use"
                else "needs_visual_review",
                window_note=override.get("note", ""),
            )
            overlap_ids, other_windows = [], []
            for other in all_events:
                if other["id"] == event["id"] or other["magnitude"] < 6:
                    continue
                other_start, other_end, _ = event_geometry(
                    other, float(prow.Lat), float(prow.Lon), dispersion
                )
                if other_start < end and other_end > start:
                    overlap_ids.append(other["id"])
                    other_windows.append((other_start, other_end))
            item["other_predicted_arrival_ids"] = overlap_ids

            arrivals = TauPyModel(model="iasp91").get_travel_times(
                event["depth_km"],
                distance / 111.195,
                phase_list=["P", "Pdiff", "PKP", "S", "Sdiff", "SKS"],
            )
            p_arrivals = [a for a in arrivals if a.name.startswith("P")]
            s_arrivals = [a for a in arrivals if a.name.startswith("S")]
            item["predicted_first_p_utc"] = (
                str(origin + p_arrivals[0].time) if p_arrivals else None
            )
            item["predicted_first_s_utc"] = (
                str(origin + s_arrivals[0].time) if s_arrivals else None
            )
            group_band = (dispersion[:, 0] >= band[0]) & (dispersion[:, 0] <= band[1])
            group_f = np.r_[band[0], dispersion[group_band, 0], band[1]]
            speeds = np.interp(group_f, dispersion[:, 0], dispersion[:, 2])
            group_range = (
                distance * 1000 / speeds.max(),
                distance * 1000 / speeds.min(),
            )
            overview = output / f"{event['id']}_arrival_review.png"
            item["group_arrival_range_seconds"] = list(map(float, group_range))
            overview_seconds = np.ceil(float(end - origin + 900) / 60) * 60
            event_sources = source_metadata(
                records, origin, origin + overview_seconds + 1
            )
            try:
                plot_arrival_review(
                    records,
                    rows,
                    channels,
                    origin,
                    start,
                    end,
                    band,
                    group_range,
                    item["predicted_first_p_utc"],
                    item["predicted_first_s_utc"],
                    other_windows,
                    pressure_sign,
                    overview,
                )
                item.setdefault("files", {})["arrival_review"] = str(overview)
            except ValueError as exc:
                item["overview_reason"] = str(exc)
            quiet_start = None
            if tilt_method == "none":
                direction = np.array([1.0])
                angle, reduction = None, None
                tilt_info = {
                    "method": "none",
                    "applied": False,
                    "vertical_component": rows[-1].Seed_id,
                    "quiet_training": "not_requested",
                }
                meta = source_metadata(
                    records,
                    start - (32 / event_resample_rate if event_resample_rate else 0),
                    end
                    + (32 / event_resample_rate if event_resample_rate else 0)
                    + 2 / min(float(c.sample_rate) for c in channels),
                )
                event_sources.extend(meta)
            else:
                quiet_start = None
                for hours in (6, 12, 24, 48):
                    left, right = origin - hours * 3600, origin - (hours - 2) * 3600
                    # Conservative catalogue screen; quietness still needs review.
                    overlaps = False
                    for other in all_events:
                        os, oe, _ = event_geometry(
                            other, float(prow.Lat), float(prow.Lon), dispersion
                        )
                        coda = max(3600, (other["magnitude"] - 5.5) * 43200)
                        if os < right and oe + coda > left:
                            overlaps = True
                            break
                    if not overlaps:
                        quiet_start = left
                        break
                if quiet_start is None:
                    raise ValueError("No catalogue-screened quiet training interval")
                meta = source_metadata(
                    records,
                    quiet_start,
                    quiet_start + 7200 + 2 / float(channels[0].sample_rate),
                ) + source_metadata(
                    records, start, end + 2 / float(channels[0].sample_rate)
                )
                event_sources.extend(meta)
                timing_context = "quiet_training"
                quiet, quiet_actual, rate = load_aligned(
                    records, rows, channels, quiet_start, 7200
                )
                direction, angle, reduction, tilt_info = estimate_tilt(
                    quiet, rows, channels, rate, quiet_actual, tilt_method
                )
            item["tilt"] = tilt_info
            timing_context = "event_window"
            arrays, actual, rate = load_aligned(
                records,
                rows,
                channels,
                start,
                window_seconds,
                event_resample_rate,
                native_rate_tolerance_ppm,
            )
            if event_resample_rate is not None:
                item["resampling"] = summary["settings"]["event_resampling"]
            spectra = event_spectra(
                arrays, rows, channels, rate, direction, 200, band, pressure_sign
            )
            # Separate display spectra; the original fit/phase/QC band is unchanged.
            display_band = (0.005, min(0.1, rate / 2 * (1 - 1e-8)))
            try:
                paper_spectra = event_spectra(
                    arrays,
                    rows,
                    channels,
                    rate,
                    direction,
                    200,
                    display_band,
                    pressure_sign,
                )
            except ValueError as exc:
                # A display-only response problem must not change fit-band QC.
                item["paper_display_reason"] = str(exc)
                paper_spectra = spectra.copy()
            spectra.update(
                {f"paper_{key}": value for key, value in paper_spectra.items()}
            )
            f = spectra["frequency_hz"]
            phase_velocity = np.interp(f, dispersion[:, 0], dispersion[:, 1])
            theory = water_ratio(f, phase_velocity, depth, density, sound_speed)
            coh, phase = spectra["coherence_squared"], spectra["phase_deg"]
            valid = (coh >= threshold) & (abs(phase) <= 20)
            item.update(
                quiet_start_utc=str(quiet_start) if quiet_start is not None else None,
                rotation_coefficients=direction.tolist(),
                tilt_angle_deg=angle,
                quiet_vertical_power_ratio=reduction,
                median_coherence_squared=float(np.median(coh)),
                mean_coherence_squared=float(np.mean(coh)),
                accepted_bins=int(valid.sum()),
                raw_median_coherence_squared=float(
                    np.median(spectra["raw_coherence_squared"])
                ),
                actual_start_utc=str(actual),
                segment_count=spectra["segment_count"],
            )
            coherent = coh >= threshold
            coherent_event = (
                np.median(coh) >= threshold
                and np.mean(coh) >= threshold
                and coherent.sum() >= 3
            )
            accepted = coherent_event and valid.sum() >= 3 and not overlap_ids
            plot_gain = None
            plot_label = None
            fit = None
            if not coherent_event:
                item.update(status="rejected_qc", reason="Insufficient band coherence")
            elif not accepted:
                item.update(
                    status="polarity_review"
                    if valid.sum() < 3
                    else "arrival_overlap_review",
                    reason="Phase convention or arrival overlap needs review",
                )
            if coherent_event:
                fit_mask = valid if accepted else coherent
                try:
                    fit, grid, loss = fit_multiplier(
                        spectra["psd_ratio_pa"][fit_mask], theory[fit_mask]
                    )
                except ValueError as exc:
                    # No gain is recorded, but spectra and the QC outcome are kept.
                    item["fit_reason"] = str(exc)
                    if accepted:
                        item.update(status="unresolved", reason=str(exc))
            if fit is not None:
                plot_gain = fit["pressure_multiplier"]
                item.update(
                    diagnostic_pressure_amplitude_multiplier=plot_gain,
                    grid_pressure_amplitude_multiplier=fit["grid_pressure_multiplier"],
                    diagnostic_relative_rms=fit["relative_rms"],
                )
                for name, observed in (
                    ("h1", abs(spectra["h1_pa"][fit_mask])),
                    ("raw", spectra["raw_psd_ratio_pa"][fit_mask]),
                ):
                    try:
                        diagnostic, _, _ = fit_multiplier(observed, theory[fit_mask])
                        item[f"diagnostic_{name}_amplitude_multiplier"] = diagnostic[
                            "pressure_multiplier"
                        ]
                    except ValueError as exc:
                        item[f"diagnostic_{name}_fit_reason"] = str(exc)
                if accepted:
                    item.update(
                        fit,
                        status="candidate",
                        signed_pressure_multiplier=pressure_sign
                        * fit["pressure_multiplier"],
                        sign_status="assumed_needs_confirmation",
                    )
                    plot_label = f"Candidate amplitude ×{plot_gain:.3f}"
                else:
                    plot_label = f"Amplitude-only ×{plot_gain:.3f}; QC unresolved"
                write_csv(
                    output / f"{event['id']}_grid.csv",
                    [
                        {"pressure_multiplier": float(g), "normalized_mse": float(v)}
                        for g, v in zip(grid, loss, strict=True)
                    ],
                )
            spectra["expected_pa"] = theory
            spectra["accepted_bins"] = valid
            npz = output / f"{event['id']}_spectra.npz"
            np.savez_compressed(npz, **spectra)
            table = output / f"{event['id']}_spectra.csv"
            write_csv(
                table,
                [
                    {
                        "frequency_hz": float(f[i]),
                        "observed_pa": float(spectra["psd_ratio_pa"][i]),
                        "expected_pa": float(theory[i]),
                        "h1_amplitude_pa": float(abs(spectra["h1_pa"][i])),
                        "coherence_squared": float(coh[i]),
                        "phase_deg": float(phase[i]),
                        "nominal_phase_deg": float(spectra["nominal_phase_deg"][i]),
                        "pressure_sign": pressure_sign,
                        "accepted": bool(valid[i]),
                    }
                    for i in range(len(f))
                ],
            )
            waveform = output / f"{event['id']}_waveform.png"
            wf, coefficients, _ = corrected_fft(arrays, rows, channels, rate, *band)
            all_f = np.fft.rfftfreq(arrays.shape[1], 1 / rate)
            wave_band = (all_f >= band[0]) & (all_f <= band[1])
            reconstructed = []
            for values in (
                pressure_sign * coefficients[0],
                2j * np.pi * wf * (direction @ coefficients[1:]),
            ):
                fourier = np.zeros(len(all_f), complex)
                fourier[wave_band] = values
                reconstructed.append(np.fft.irfft(fourier, arrays.shape[1]))
            step = max(
                1, round(rate)
            )  # Display thinning only; spectra use every sample.
            elapsed = (
                float(actual - origin) + np.arange(0, arrays.shape[1], step) / rate
            )
            wavefig, waveaxes = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
            waveaxes[0].plot(elapsed, reconstructed[0][::step])
            waveaxes[0].set_ylabel(f"Pressure ×{pressure_sign:+d} (Pa; sign assumed)")
            waveaxes[1].plot(elapsed, reconstructed[1][::step], color="tab:orange")
            waveaxes[1].set_ylabel("Upward acceleration (m/s²)")
            waveaxes[1].set_xlabel("Seconds after event origin")
            wavefig.suptitle(
                f"{event['id']} — group-velocity candidate window\n"
                "Assumed sign; band-limited display; no amplitude gain applied"
            )
            wavefig.tight_layout()
            wavefig.savefig(waveform, dpi=140)
            plt.close(wavefig)
            fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
            axes[0].plot(
                f,
                spectra["psd_ratio_pa"] / (density * depth),
                "o-",
                label="Nominal P/a",
            )
            axes[0].plot(f, theory / (density * depth), label="Water-column prediction")
            if plot_gain is not None:
                axes[0].plot(
                    f,
                    plot_gain * spectra["psd_ratio_pa"] / (density * depth),
                    "--",
                    label=plot_label,
                )
            axes[0].set_ylabel("(P/a) / (ρH)")
            axes[0].legend()
            axes[1].plot(f, coh, "o-", label=tilt_info["method"])
            axes[1].plot(
                f, spectra["raw_coherence_squared"], "--", label="Original vertical"
            )
            axes[1].axhline(threshold, color="gray", ls=":")
            axes[1].set_ylim(0, 1.02)
            axes[1].set_ylabel("Squared coherence")
            axes[1].legend()
            axes[2].plot(f, spectra["nominal_phase_deg"], "--", label="Metadata sign")
            axes[2].plot(
                f, phase, "o-", label=f"Assumed pressure sign {pressure_sign:+d}"
            )
            axes[2].legend()
            axes[2].axhline(0, color="gray", ls=":")
            axes[2].axhspan(-20, 20, color="gray", alpha=0.1)
            axes[2].set_ylabel("P / upward-a phase (°)")
            axes[2].set_xlabel("Frequency (Hz)")
            for ax in axes:
                ax.grid(alpha=0.25)
            fig.suptitle(
                f"{network}.{station} — {event['id']} Mw {event['magnitude']}\n"
                f"{item['status']}; assumptions require review; no correction applied"
            )
            fig.tight_layout()
            plot = output / f"{event['id']}_comparison.png"
            fig.savefig(plot, dpi=140)
            if show:
                plt.show()
            plt.close(fig)
            item.setdefault("files", {}).update(
                {
                    "spectra": str(npz),
                    "table": str(table),
                    "plot": str(plot),
                    "waveform": str(waveform),
                }
            )
            summary["source_files"].extend(event_sources)
            for source in meta:
                stat = Path(source["path"]).stat()
                if (stat.st_size, stat.st_mtime_ns) != (
                    source["size"],
                    source["mtime_ns"],
                ):
                    raise ValueError("Source waveform changed during analysis")
            print(
                f"{event['id']}: {item['status']}; MSC={np.median(coh):.3f}", flush=True
            )
        except (ValueError, FileNotFoundError) as exc:
            if isinstance(exc, SampleAlignmentError):
                item["timing_review"] = {
                    "context": timing_context,
                    "channels": exc.timing,
                    "phase_bound_deg_at_band_max": 360 * band[1] * exc.span_seconds,
                    "decision": "excluded; investigate offsets before correcting",
                }
            elif isinstance(exc, SampleGridError):
                item["timing_review"] = {
                    "context": timing_context,
                    "seed_id": exc.seed_id,
                    "within_channel_offset_seconds": exc.residual_seconds,
                    "within_channel_offset_samples": (
                        exc.residual_seconds * exc.sample_rate
                    ),
                    "phase_bound_deg_at_band_max": (
                        360 * band[1] * exc.residual_seconds
                    ),
                    "decision": "excluded; native record grids were not merged",
                }
            item.update(status="unresolved", reason=str(exc))
            summary["source_files"].extend(event_sources)
            print(f"{event['id']}: unresolved ({exc})", flush=True)
    review_template = {
        e["id"]: {
            "decision": "pending",
            "start_utc": e["start_utc"],
            "end_utc": e["end_utc"],
            "note": "",
        }
        for e in summary["events"]
        if "start_utc" in e
    }
    template_path = output / "window_review_template.json"
    atomic_json(template_path, review_template)
    summary["files"]["window_review_template"] = str(template_path)
    amplitude_items = [
        e for e in summary["events"] if "diagnostic_pressure_amplitude_multiplier" in e
    ]
    summary["amplitude_review_events"] = len(amplitude_items)
    if amplitude_items:
        values = np.array(
            [e["diagnostic_pressure_amplitude_multiplier"] for e in amplitude_items]
        )
        summary["diagnostic_median_pressure_amplitude_multiplier"] = float(
            np.median(values)
        )
        summary["diagnostic_event_multiplier_std"] = (
            float(np.std(values, ddof=1)) if len(values) > 1 else None
        )
    accepted_items = [e for e in summary["events"] if e["status"] == "candidate"]
    summary["candidate_events"] = len(accepted_items)
    if len(accepted_items) >= min_events:
        values = np.array([e["pressure_multiplier"] for e in accepted_items])
        summary.update(
            status="candidate_comparison_under_assumed_sign",
            median_pressure_multiplier=float(np.median(values)),
            event_multiplier_std=float(np.std(values, ddof=1)),
            uncertainty_kind="event scatter, not a confidence interval",
        )
    else:
        summary["reason"] = (
            f"Fewer than {min_events} QC-accepted independent events; "
            "no station-level estimate"
        )
    plot_paper_comparison(
        summary, dispersion, output / "paper_figure5_comparison.png", show=show
    )
    for source in summary["source_files"]:
        stat = Path(source["path"]).stat()
        if (stat.st_size, stat.st_mtime_ns) != (source["size"], source["mtime_ns"]):
            raise ValueError("Source waveform changed during review")
    for path, digest in (
        summary["source_inventory_sha256"]
        | summary["response_snapshot_sha256"]
        | summary["inputs"]
    ).items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError("Calibration input changed during analysis")
    summary["files"]["result"] = str(output / "comparison.json")
    atomic_json(output / "comparison.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    model = parser.add_mutually_exclusive_group(required=True)
    model.add_argument("--dispersion-file", type=Path)
    model.add_argument(
        "--earth-model-file",
        type=Path,
        help="Layer CSV (km, km/s, km/s, g/cm3); compute dispersion with disba",
    )
    parser.add_argument("--master-dir", type=Path, default=Path("exports"))
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument("--network", default="YV")
    parser.add_argument("--station", default="RR52")
    parser.add_argument(
        "--depth", type=float, help="Default: inventory elevation estimate"
    )
    parser.add_argument("--density", type=float, default=1028)
    parser.add_argument("--sound-speed", type=float, default=1500)
    parser.add_argument("--frequency-min", type=float, default=0.02)
    parser.add_argument("--frequency-max", type=float, default=0.06)
    parser.add_argument("--threshold", type=float, default=0.95)
    parser.add_argument(
        "--max-events",
        type=int,
        default=8,
        help="Attempt cap per magnitude level, only with --event-sampling bounded",
    )
    parser.add_argument(
        "--min-events",
        type=int,
        default=3,
        help="Minimum accepted events for station estimate; not a cap in all mode",
    )
    parser.add_argument(
        "--event-sampling",
        choices=("all", "bounded"),
        default="all",
        help="all: exhaust eligible catalogue; bounded: historical pilot limits",
    )
    parser.add_argument(
        "--no-magnitude-fallback",
        dest="magnitude_fallback",
        action="store_false",
        help="Only search Mw>7; explicit --event-id also disables fallback",
    )
    parser.add_argument("--event-id")
    parser.add_argument(
        "--window-overrides",
        type=Path,
        help="Reviewed UTC windows keyed by event ID; no automatic picking",
    )
    parser.add_argument("--show", action="store_true")
    parser.add_argument("--native-rate-tolerance-ppm", type=float, default=0)
    parser.add_argument(
        "--event-resample-rate",
        type=float,
        help="Explicit working-copy anti-alias rate in Hz; requires --tilt-method none",
    )
    parser.add_argument(
        "--tilt-method",
        choices=("none", "tiskitpy", "local_eigen"),
        default="tiskitpy",
        help="none: pressure/BHZ only, without quiet training or tilt correction",
    )
    parser.add_argument(
        "--pressure-sign",
        type=int,
        choices=(-1, 1),
        help="RR52 default -1 (assumed); supply an explicit sign for other stations",
    )
    options = vars(parser.parse_args())
    options["band"] = (options.pop("frequency_min"), options.pop("frequency_max"))
    result = run_calibration(**options)
    print(result["files"]["result"])


if __name__ == "__main__":
    main()
