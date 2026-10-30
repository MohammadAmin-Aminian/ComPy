"""Local-earthquake P-wave pressure/vertical check for one channel pair.

Method: docs/local_p_method.md. Reads waveforms, StationXML and a catalogue in
place and writes only into a fresh results directory. It never applies a gain,
assigns a score or updates a master. Steady-state (Deng Eq. 5) windows give
conditional gain candidates; pre-multiple windows (p = rho c0 v / eta, time
gated) are diagnostics only. Every rejected event keeps its reason.
"""

import argparse
import csv
import hashlib
import json
import os
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from math import gcd
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from obspy import UTCDateTime, read, read_inventory  # noqa: E402
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees  # noqa: E402
from obspy.taup import TauPyModel  # noqa: E402
from scipy.signal import butter, sosfilt  # noqa: E402
from scipy.signal.windows import tukey  # noqa: E402

from compy_pressure import local_p_physics as lp  # noqa: E402
from compy_pressure.archive import WaveformFile  # noqa: E402
from compy_pressure.common import PROJECT_DIRECTORY  # noqa: E402
from compy_pressure.local_p_io import local_p_output  # noqa: E402
from compy_pressure.ppsd import atomic_json, response_handling  # noqa: E402
from compy_pressure.rayleigh_calibration import (  # noqa: E402
    SampleAlignmentError,
    catalog_events,
    load_aligned,
    load_resampled,
)

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

# Provenance: the git repository holding this code, and the paths in it
# whose uncommitted changes are recorded.
ROOT = Path(__file__).resolve().parent
CODE_PATHS = ["."]
PROTOCOL_VERSION = "local-p-v2-noise-guard-20261016"
DEFAULTS = {
    "max_distance_km": 300.0,
    "min_hypocentral_km": 10.0,
    "max_events": 80,
    "taup_model": "iasp91",
    "search_min_s": 3.0,
    "search_fraction": 0.10,
    "pick_band_hz": [2.0, 10.0],
    "pick_snr_min": 3.0,
    "onset_tolerance_s": 0.10,
    "lead_s": 0.2,
    "noise_gap_s": 1.0,
    "s_fraction": 0.85,
    "s_margin_s": 1.0,
    "min_window_s": 1.0,
    "max_window_s": 20.0,
    "steady_round_trips": lp.STEADY_STATE_ROUND_TRIPS,
    "pre_multiple_safety": 0.95,
    "nw": 2.5,
    "tapers": 4,
    "coherence_min": 0.8,
    "snr_power_min": 10.0,
    "response_floor": 0.1,
    "nyquist_fraction": 0.8,
    "resonance_margin": lp.RESONANCE_MARGIN,
    "bands_hz": [[0.5, 1.0], [1.0, 2.0], [2.0, 4.0], [4.0, 8.0], [8.0, 16.0]],
    "min_bins_per_band": 3,
    "min_events": 3,
    "depth_uncertainty": 0.05,
    "density": lp.DENSITY,
    "sound_speed": lp.SOUND_SPEED,
    "density_range": list(lp.DENSITY_RANGE),
    "sound_speed_range": list(lp.SOUND_SPEED_RANGE),
    "p_velocity_min": lp.P_VELOCITY_MIN,
    "polarity_tolerance_deg": 45.0,
    "timing_lag_search_s": 0.0,
    "lag_band_hz": [0.5, 2.0],
    "lag_refine_max_s": 3.0,
    # Optional clock-drift model for the pressure logger (timing only):
    # {"reference_utc", "intercept_s", "slope_s_per_day", "valid_from_utc",
    #  "valid_to_utc", "source"}; lag = intercept + slope * days since reference.
    "timing_drift": None,
    "native_rate_tolerance_ppm": 0.1,
    "overlap_before_s": 30.0,
    "screening_hours_m6": 3.0,
    "screening_hours_m55": 1.0,
    "event_plots": 6,
}
PROTECTED_DEFAULTS = (
    Path.home() / "Desktop/Data_Bruit_FM",
    Path.home() / "Projects/bruit-fm/exports",
    Path.home() / "Projects/bruit-fm/results",
    PROJECT_DIRECTORY / "exports",
)
VERTICAL_UNITS = {"M/S"}


class EventRejected(Exception):
    """An event is unusable for a recorded, explicit reason."""

    def __init__(self, reason, detail=""):
        super().__init__(f"{reason}: {detail}" if detail else reason)
        self.reason = reason
        self.detail = detail


# --------------------------------------------------------------------------
# Inputs and provenance


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def file_state(path):
    stat = Path(path).stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def git_state():
    try:
        commit = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(ROOT), "status", "--porcelain", "--", *CODE_PATHS],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        return {"commit": commit, "dirty": bool(dirty), "dirty_paths": CODE_PATHS}
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None, "dirty": None, "dirty_paths": CODE_PATHS}


def load_local_catalog(path):
    """ComCat CSV (with uncertainty columns) or repository GeoJSON."""
    path = Path(path)
    if path.suffix.lower() == ".json":
        return [
            {
                "id": e["id"],
                "time": UTCDateTime(e["origin_utc"]),
                "latitude": float(e["latitude"]),
                "longitude": float(e["longitude"]),
                "depth_km": float(e["depth_km"]),
                "magnitude": float(e["magnitude"]),
                "magnitude_type": e["magnitude_type"],
                "type": "earthquake",
                "horizontal_error_km": None,
                "depth_error_km": None,
                "network": "",
            }
            for e in catalog_events(path)
        ]

    def number(value):
        return float(value) if value not in ("", None) else None

    events = []
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            events.append(
                {
                    "id": row["id"],
                    "time": UTCDateTime(row["time"]),
                    "latitude": float(row["latitude"]),
                    "longitude": float(row["longitude"]),
                    "depth_km": float(row["depth"]),
                    "magnitude": number(row["mag"]),
                    "magnitude_type": row.get("magType", ""),
                    "type": row.get("type", ""),
                    "horizontal_error_km": number(row.get("horizontalError")),
                    "depth_error_km": number(row.get("depthError")),
                    "network": row.get("net", ""),
                }
            )
    return events


def channel_epochs(inventory_path, seed_id):
    net, sta, loc, cha = seed_id.split(".")
    inventory = read_inventory(str(inventory_path))
    found = []
    for n in inventory.select(network=net, station=sta, location=loc, channel=cha):
        for s in n:
            for c in s:
                found.append(SimpleNamespace(seed_id=seed_id, station=s, channel=c))
    if not found:
        raise ValueError(f"No StationXML epoch for {seed_id}")
    return found


def epoch_for(epochs, start, end):
    covering = [
        e
        for e in epochs
        if (e.channel.start_date is None or e.channel.start_date <= start)
        and (e.channel.end_date is None or end <= e.channel.end_date)
    ]
    if len(covering) != 1:
        raise EventRejected(
            "no_unique_response_epoch", f"{len(covering)} epochs cover the window"
        )
    return covering[0]


def water_depth(epoch):
    for item in (epoch.channel, epoch.station):
        elevation = getattr(item, "elevation", None)
        if elevation is not None and elevation < 0:
            return -float(elevation)
    return None


def duration_quantum(rates):
    """Smallest duration (s) holding an integer number of samples at every rate."""
    numerators, denominators = [], []
    for rate in rates:
        period = 1 / Fraction(str(rate))
        numerators.append(period.numerator)
        denominators.append(period.denominator)
    lcm_num = 1
    for value in numerators:
        lcm_num = lcm_num * value // gcd(lcm_num, value)
    gcd_den = 0
    for value in denominators:
        gcd_den = gcd(gcd_den, value)
    return float(Fraction(lcm_num, gcd_den))


def quantize(seconds, quantum, down=True):
    steps = np.floor(seconds / quantum + 1e-9) if down else np.ceil(seconds / quantum)
    return float(steps * quantum)


# --------------------------------------------------------------------------
# Waveform files: header-identified daily files from SDS or daily layouts


class DayFileIndex:
    """Map days to candidate files by filename, then trust only the headers.

    Filenames may carry location 00/-- while headers are blank (XO, 7D); the
    header SEED ID is used, and only that location discrepancy is accepted.
    """

    def __init__(self, roots, seed_id):
        self.seed_id = seed_id
        self.net, self.sta, self.loc, self.cha = seed_id.split(".")
        self.by_day = {}
        self.headers = {}
        self.states = {}
        self.duplicate_layout_days = set()
        self.location_mismatch = False
        for root in roots:
            for directory, _, names in os.walk(Path(root).expanduser()):
                for name in names:
                    day = self._day(name)
                    if day:
                        self.by_day.setdefault(day, []).append(Path(directory, name))
        if not self.by_day:
            raise FileNotFoundError(f"No waveform files for {seed_id} under {roots}")

    def _day(self, name):
        parts = name.split(".")
        if parts and parts[-1].lower() == "mseed" and len(parts) in (6, 8):
            body = parts[:-1]
        else:
            body = parts
        if len(body) < 5 or body[0] != self.net or body[1] != self.sta:
            return None
        if body[3] != self.cha or body[2] not in {self.loc, "", "--", "00"}:
            return None
        if (
            len(body) == 7
            and body[4] == "D"
            and body[5].isdigit()
            and body[6].isdigit()
        ):
            return f"{int(body[5]):04d}.{int(body[6]):03d}"
        if len(body) == 5 and len(body[4]) == 8 and body[4].isdigit():
            date = UTCDateTime(f"{body[4][:4]}-{body[4][4:6]}-{body[4][6:]}")
            return f"{date.year:04d}.{date.julday:03d}"
        return None

    def records(self, start, end):
        days, day = [], UTCDateTime(start.date) - 86400
        while day <= end:
            days.append(f"{day.year:04d}.{day.julday:03d}")
            day += 86400
        out = []
        for key in days:
            paths = sorted(self.by_day.get(key, []))
            # The same day may exist in SDS (CHA.D/) and daily (waveforms/) layouts.
            sds = [p for p in paths if p.parent.name.endswith(".D")]
            if sds and len(sds) < len(paths):
                self.duplicate_layout_days.add(key)
                paths = sds
            for path in paths:
                if path not in self.headers:
                    self.states[str(path)] = file_state(path)
                    self.headers[path] = self._read_headers(path)
                out.extend(
                    r for r in self.headers[path] if r.start < end and r.end >= start
                )
        return out

    def _read_headers(self, path):
        found = []
        for trace in read(str(path), format="MSEED", headonly=True):
            net, sta, loc, cha = trace.id.split(".")
            if (net, sta, cha) != (self.net, self.sta, self.cha):
                continue
            if loc != self.loc:
                if {loc, self.loc} - {"", "00", "--"}:
                    continue
                self.location_mismatch = True
            found.append(
                WaveformFile(
                    path,
                    self.seed_id,
                    trace.stats.starttime,
                    trace.stats.endtime,
                    float(trace.stats.sampling_rate),
                )
            )
        return found


# --------------------------------------------------------------------------
# Event geometry


@dataclass
class Geometry:
    distance_km: float
    hypocentral_km: float
    back_azimuth: float
    p_travel_s: float
    s_travel_s: float
    p_phase: str
    q_s_per_m: float
    eta: float


def event_geometry(model, event, latitude, longitude):
    distance_m, _, back_azimuth = gps2dist_azimuth(
        event["latitude"], event["longitude"], latitude, longitude
    )
    distance_km = distance_m / 1000
    depth_km = max(float(event["depth_km"]), 0.0)
    arrivals = model.get_travel_times(
        source_depth_in_km=depth_km,
        distance_in_degree=kilometers2degrees(distance_km),
        phase_list=["p", "P", "Pn", "Pg", "s", "S", "Sn", "Sg"],
    )
    p = [a for a in arrivals if a.name.upper().startswith("P")]
    s = [a for a in arrivals if a.name.upper().startswith("S")]
    if not p or not s:
        raise EventRejected("no_predicted_p_or_s")
    first_p = min(p, key=lambda a: a.time)
    q = lp.slowness_from_ray_parameter(first_p.ray_param)
    return Geometry(
        distance_km=distance_km,
        hypocentral_km=float(np.hypot(distance_km, depth_km)),
        back_azimuth=back_azimuth,
        p_travel_s=float(first_p.time),
        s_travel_s=float(min(a.time for a in s)),
        p_phase=first_p.name,
        q_s_per_m=q,
        eta=float(lp.vertical_cosine(q)),
    )


def search_half_width(geometry, event, options):
    """Pre-registered onset search: model and catalogue location uncertainty."""
    width = max(
        options["search_min_s"], options["search_fraction"] * geometry.p_travel_s
    )
    for key in ("horizontal_error_km", "depth_error_km"):
        if event.get(key):
            width += event[key] / 6.0  # km at ~6 km/s
    return float(width)


# --------------------------------------------------------------------------
# Signal helpers


def causal_band(data, rate, band):
    high = min(band[1], 0.4 * rate)
    if band[0] >= high:
        raise EventRejected("pick_band_above_nyquist")
    sos = butter(4, [band[0], high], btype="bandpass", fs=rate, output="sos")
    return sosfilt(sos, data - np.mean(data))


def measure_onset(data, rate, start_index, stop_index, band, snr_min):
    """AIC onset index within [start, stop) and its amplitude SNR."""
    filtered = causal_band(data, rate, band)
    segment = filtered[start_index:stop_index]
    index = lp.aic_onset(segment, minimum=max(10, int(0.5 * rate)))
    if index is None:
        return None, 0.0
    onset = start_index + index
    pre = filtered[max(0, onset - int(3 * rate)) : onset - int(0.5 * rate)]
    post = filtered[onset : onset + int(rate)]
    if len(pre) < rate or len(post) < rate / 2:
        return onset, 0.0
    snr = float(np.sqrt(np.mean(post**2) / max(np.mean(pre**2), 1e-300)))
    return onset, snr


def flat_topped(data, run=3):
    """Samples pinned at the window extreme for >= run consecutive samples."""
    for extreme in (np.max(data), np.min(data)):
        hits = np.flatnonzero(data == extreme)
        if len(hits) >= run and np.any(
            np.convolve(np.diff(hits) == 1, [1] * (run - 1), "valid") == run - 1
        ):
            return True
    return False


def envelope_lag(p, v, rate, band, max_lag):
    """Coarse lag (s) maximizing envelope correlation; positive = pressure later."""
    from scipy.signal import hilbert

    a = np.abs(hilbert(causal_band(p, rate, band)))
    b = np.abs(hilbert(causal_band(v, rate, band)))
    a, b = a - a.mean(), b - b.mean()
    n = len(a)
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    xc = np.fft.irfft(np.fft.rfft(a, nfft) * np.conj(np.fft.rfft(b, nfft)), nfft)
    lags = np.concatenate([np.arange(0, n), np.arange(-n + 1, 0)])
    values = np.concatenate([xc[:n], xc[nfft - n + 1 :]])
    keep = np.abs(lags) <= max_lag * rate
    best = lags[keep][np.argmax(values[keep])]
    norm = np.sqrt(np.sum(a**2) * np.sum(b**2))
    return best / rate, float(np.max(values[keep]) / norm) if norm else 0.0


# --------------------------------------------------------------------------
# Spectral analysis of one event


def model_transfer(f, regime, depth, eta, options, density=None, sound_speed=None):
    """Complex p/v model: steady state i2pi f rho h F, or rho c0 / eta.

    NumPy's forward FFT makes d/dt -> +i 2 pi f, so p = rho h F a gives
    p/v = i 2 pi f rho h F. Steady state is nan beyond the resonance margin.
    """
    density = options["density"] if density is None else density
    c0 = options["sound_speed"] if sound_speed is None else sound_speed
    q = np.sqrt(max(1 - eta**2, 0.0)) / c0
    if regime == "steady_state":
        theta = lp.water_phase(f, depth, q, c0)
        ok = theta <= options["resonance_margin"] * np.pi / 2
        out = np.full(f.shape, np.nan + 0j)
        out[ok] = 1j * 2 * np.pi * f[ok] * density * depth * lp.tan_ratio(theta[ok])
        return out
    return np.full(f.shape, density * c0 / eta + 0j)


def model_ratio_interval(f, regime, depth, eta, eta_bounds, options):
    """Min/max of |H_corner/H_nominal| over depth, eta, c0 and rho corners."""
    nominal = np.abs(model_transfer(f, regime, depth, eta, options))
    span = options["depth_uncertainty"]
    ratios = []
    for hh in (depth * (1 - span), depth * (1 + span)):
        for ee in eta_bounds:
            for c0 in options["sound_speed_range"]:
                for rho in options["density_range"]:
                    with np.errstate(invalid="ignore", divide="ignore"):
                        corner = model_transfer(f, regime, hh, ee, options, rho, c0)
                        ratios.append(np.abs(corner) / nominal)
    ratios = np.array(ratios)
    bad = ~np.isfinite(ratios).all(axis=0)
    ratios[:, bad] = np.nan
    low = np.full(f.shape, np.nan)
    high = np.full(f.shape, np.nan)
    low[~bad] = ratios[:, ~bad].min(axis=0)
    high[~bad] = ratios[:, ~bad].max(axis=0)
    return low, high


def passband_weight(f, edges):
    """Cosine-tapered weight: 0 below e0, 1 in [e1, e2], 0 above e3."""
    e0, e1, e2, e3 = edges
    weight = np.zeros_like(f)
    rise = (f > e0) & (f < e1)
    weight[rise] = 0.5 - 0.5 * np.cos(np.pi * (f[rise] - e0) / (e1 - e0))
    weight[(f >= e1) & (f <= e2)] = 1.0
    fall = (f > e2) & (f < e3)
    weight[fall] = 0.5 + 0.5 * np.cos(np.pi * (f[fall] - e2) / (e3 - e2))
    return weight


def model_passband(regime, depth, rate, options):
    """Edges where the forward model is defined (lowest possible first pole)."""
    if regime == "steady_state":
        pole = options["sound_speed_range"][0] / (
            4 * depth * (1 + options["depth_uncertainty"])
        )
        top = options["resonance_margin"] * pole
        return (0.1, 0.2, 0.85 * top, top)
    top = options["nyquist_fraction"] * rate / 2
    return (0.1, 0.2, 0.9 * top, top)


def forward_model(p_counts, v_counts, rate, responses, regime, depth, eta, options):
    """Predict pressure counts from vertical counts through model and responses.

    Observed pressure and vertical are band-limited by the same zero-phase
    weight, so observed and predicted pressure share spectral support. The
    identical multitaper estimator is then applied to both: windowing and
    leakage of a frequency-varying model affect both alike.
    """
    n = len(v_counts)
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    taper = tukey(n, 0.05)
    f = np.fft.rfftfreq(nfft, 1 / rate)
    edges = model_passband(regime, depth, rate, options)
    weight = passband_weight(f, edges)
    rp, rv, _, _ = responses(f)
    with np.errstate(invalid="ignore", divide="ignore"):
        model = model_transfer(f, regime, depth, eta, options)
        scale = np.where(weight > 0, model * rp / rv * weight, 0)
    scale[~np.isfinite(scale)] = 0

    def spectrum(x):
        x = np.asarray(x, dtype=float)
        trend = np.polyval(np.polyfit(np.arange(n), x, 1), np.arange(n))
        return np.fft.rfft((x - trend) * taper, nfft)

    v_spec = spectrum(v_counts)
    predicted = np.fft.irfft(v_spec * scale, nfft)[:n]
    observed = np.fft.irfft(spectrum(p_counts) * weight, nfft)[:n]
    vertical = np.fft.irfft(v_spec * weight, nfft)[:n]
    return observed, predicted, vertical, edges


def analyse_event(
    p_seg,
    v_seg,
    signal_slice,
    noise_slice,
    rate,
    responses,
    regime,
    depth,
    eta,
    options,
):
    """Forward-model comparison, QC mask and band summaries for one event."""
    observed, predicted, vertical, edges = forward_model(
        p_seg, v_seg, rate, responses, regime, depth, eta, options
    )
    n = signal_slice.stop - signal_slice.start
    nfft = 1 << int(np.ceil(np.log2(4 * n)))
    nw, k = options["nw"], options["tapers"]

    def products(sl):
        f, ok, vk = lp.multitaper_spectra(observed[sl], vertical[sl], rate, nw, k, nfft)
        _, pk, _ = lp.multitaper_spectra(predicted[sl], vertical[sl], rate, nw, k, nfft)
        obs = lp.cross_spectral_products(ok, vk)
        pred = lp.cross_spectral_products(pk, vk)
        return f, obs, pred

    f, signal, model_sig = products(signal_slice)
    _, noise, model_noise = products(noise_slice)
    _, _, rp_rel, rv_rel = responses(f)
    half_bandwidth = nw * rate / n
    eta_bounds = (
        lp.eta_min_for_p(options["sound_speed"], options["p_velocity_min"]),
        1.0,
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        snr_p = signal["spp"] / noise["spp"]
        snr_v = signal["svv"] / noise["svv"]
        # g * p_nominal = p_model; both pressures are counts through the same response.
        g_h1 = np.abs(model_sig["spv"]) / np.abs(signal["spv"])
        g_h2 = np.abs(model_sig["h1"]) / np.abs(signal["h2"])
        g_nc = np.sqrt(
            (model_sig["spp"] - model_noise["spp"]) / (signal["spp"] - noise["spp"])
        )
        g_nc[
            (signal["spp"] <= noise["spp"]) | (model_sig["spp"] <= model_noise["spp"])
        ] = np.nan
        low, high = model_ratio_interval(f, regime, depth, eta, eta_bounds, options)
        phase = np.angle(signal["spv"] / model_sig["spv"], deg=True)
    checks = {
        "coherence": signal["coherence"] >= options["coherence_min"],
        "snr_pressure": snr_p >= options["snr_power_min"],
        "snr_vertical": snr_v >= options["snr_power_min"],
        "response_floor": (rp_rel >= options["response_floor"])
        & (rv_rel >= options["response_floor"]),
        "nyquist": f <= options["nyquist_fraction"] * rate / 2,
        "resolution": f >= 2 * half_bandwidth,
        "model_passband": (f >= edges[1] + 2 * half_bandwidth)
        & (f <= edges[2] - 2 * half_bandwidth),
        "model_interval_defined": np.isfinite(low) & np.isfinite(high),
    }
    valid = np.logical_and.reduce(list(checks.values())) & np.isfinite(g_nc)
    bands = []
    for lo, hi in options["bands_hz"]:
        sel = valid & (f >= lo) & (f < hi)
        count = int(sel.sum())
        independent = (hi - lo) / (2 * half_bandwidth)
        item = {
            "band_hz": [lo, hi],
            "valid_bins": count,
            "independent_estimates": independent,
        }
        if count >= options["min_bins_per_band"] and independent >= 0.5:
            item.update(
                {
                    "gain_noise_corrected": float(np.median(g_nc[sel])),
                    "gain_h1_upper": float(np.median(g_h1[sel])),
                    "gain_h2_lower": float(np.median(g_h2[sel])),
                    "gain_model_low": float(np.median(g_nc[sel] * low[sel])),
                    "gain_model_high": float(np.median(g_nc[sel] * high[sel])),
                    "phase_deg": lp.circular_mean_deg(phase[sel]),
                    "coherence_median": float(np.median(signal["coherence"][sel])),
                    "snr_pressure_median": float(np.median(snr_p[sel])),
                    "snr_vertical_median": float(np.median(snr_v[sel])),
                }
            )
            item["polarity"] = int(
                lp.polarity_class(
                    [item["phase_deg"]], options["polarity_tolerance_deg"]
                )[0]
            )
        bands.append(item)
    arrays = {
        "frequency_hz": f,
        "s_obs_v": signal["spv"],
        "s_pred_v": model_sig["spv"],
        "s_obs_obs": signal["spp"],
        "s_pred_pred": model_sig["spp"],
        "s_vv": signal["svv"],
        "n_obs_obs": noise["spp"],
        "n_pred_pred": model_noise["spp"],
        "n_vv": noise["svv"],
        "coherence": signal["coherence"],
        "gain_noise_corrected": g_nc,
        "gain_h1_upper": g_h1,
        "gain_h2_lower": g_h2,
        "model_ratio_low": low,
        "model_ratio_high": high,
        "phase_obs_minus_pred_deg": phase,
        "valid": valid,
        **{f"check_{name}": mask for name, mask in checks.items()},
    }
    failing = {
        name: int(
            (
                ~mask
                & (f >= options["bands_hz"][0][0])
                & (f < options["bands_hz"][-1][1])
            ).sum()
        )
        for name, mask in checks.items()
    }
    return (
        bands,
        arrays,
        {
            "half_bandwidth_hz": half_bandwidth,
            "model_passband_hz": list(edges),
            "bins_failing_in_bands": failing,
        },
    )


# --------------------------------------------------------------------------
# One station


def responses_for(p_epoch, v_epoch):
    _, p_kind = response_handling(p_epoch.channel)
    _, v_kind = response_handling(v_epoch.channel)
    if p_kind != "Pa^2/Hz":
        raise ValueError("Pressure channel response is not in Pa")
    units = v_epoch.channel.response.instrument_sensitivity.input_units
    if units.upper().replace(" ", "") not in VERTICAL_UNITS:
        raise ValueError(f"Vertical response must be velocity (m/s), not {units}")
    if v_epoch.channel.dip is None or not np.isclose(
        abs(v_epoch.channel.dip), 90, atol=0.1
    ):
        raise ValueError("Vertical channel dip is not +-90 degrees")
    # SEED dip -90 means positive output for upward motion.
    v_sign = -1.0 if v_epoch.channel.dip > 0 else 1.0
    cache = {}

    def evaluate(f):
        key = (len(f), float(f[1]) if len(f) > 1 else 0.0)
        if key not in cache:
            g = f.copy()
            g[0] = g[1] if len(g) > 1 else 1.0
            rp = p_epoch.channel.response.get_evalresp_response_for_frequencies(
                g, output="DEF"
            )
            rv = (
                v_sign
                * v_epoch.channel.response.get_evalresp_response_for_frequencies(
                    g, output="VEL"
                )
            )
            sp = abs(p_epoch.channel.response.instrument_sensitivity.value)
            sv = abs(v_epoch.channel.response.instrument_sensitivity.value)
            cache[key] = (rp, rv, np.abs(rp) / sp, np.abs(rv) / sv)
        return cache[key]

    return evaluate


def load_pair(p_index, v_index, p_epoch, v_epoch, start, seconds, rate, lag, options):
    """Common-grid counts for pressure (shifted by an integer-sample lag) and Z."""
    rows = [
        SimpleNamespace(Seed_id=p_index.seed_id),
        SimpleNamespace(Seed_id=v_index.seed_id),
    ]
    native = [float(p_epoch.channel.sample_rate), float(v_epoch.channel.sample_rate)]
    tol = options["native_rate_tolerance_ppm"]
    p_records = p_index.records(start + lag - 120, start + lag + seconds + 120)
    v_records = v_index.records(start - 120, start + seconds + 120)
    for records, nominal in ((p_records, native[0]), (v_records, native[1])):
        if not records:
            raise EventRejected("missing_waveform_file")
        for record in records:
            if abs(record.sample_rate / nominal - 1) * 1e6 > max(tol, 1e-6):
                raise EventRejected(
                    "response_sample_rate_mismatch",
                    f"header {record.sample_rate} Hz vs response {nominal} Hz",
                )
    exact = all(
        r.sample_rate == native[i]
        for i, rr in enumerate((p_records, v_records))
        for r in rr
    )
    try:
        if lag == 0 and native[0] == native[1] == rate and exact:
            try:
                arrays, actual, _ = load_aligned(
                    p_records + v_records,
                    rows,
                    [p_epoch.channel, v_epoch.channel],
                    start,
                    seconds,
                )
                return arrays[0], arrays[1], "native_aligned", actual
            except SampleAlignmentError:
                pass
        p_arr, _, _ = load_resampled(
            p_records, rows[:1], [p_epoch.channel], start + lag, seconds, rate, tol
        )
        v_arr, _, _ = load_resampled(
            v_records, rows[1:], [v_epoch.channel], start, seconds, rate, tol
        )
        return p_arr[0], v_arr[0], "resampled_common_grid", start
    except ValueError as exc:
        reason = "gap_or_missing_data"
        if "epoch" in str(exc):
            reason = "window_crosses_response_epoch"
        raise EventRejected(reason, str(exc)) from exc


def choose_window(geometry, onset_offset_s, depth, options):
    """Regime and window length after the measured onset (method note 3.3)."""
    sp = geometry.s_travel_s - geometry.p_travel_s
    s_limit = options["s_fraction"] * sp - options["s_margin_s"] + options["lead_s"]
    s_limit -= max(onset_offset_s, 0.0)  # a late onset leaves less time before S
    h = (
        depth * (1 - options["depth_uncertainty"]),
        depth * (1 + options["depth_uncertainty"]),
    )
    eta_bounds = (
        lp.eta_min_for_p(options["sound_speed"], options["p_velocity_min"]),
        1.0,
    )
    c_bounds = tuple(options["sound_speed_range"])
    tau_min = 2 * h[0] * eta_bounds[0] / c_bounds[1]
    tau_max = 2 * h[1] * eta_bounds[1] / c_bounds[0]
    steady_length = min(s_limit, options["max_window_s"])
    if steady_length >= options["steady_round_trips"] * tau_max + options["lead_s"]:
        return "steady_state", steady_length, tau_min, tau_max
    gated = min(options["pre_multiple_safety"] * tau_min + options["lead_s"], s_limit)
    if gated >= options["min_window_s"]:
        return "pre_multiple", gated, tau_min, tau_max
    reason = (
        "window_regime_mixed" if s_limit >= options["min_window_s"] else "s_too_close"
    )
    raise EventRejected(
        reason, f"S-limited length {s_limit:.2f} s, tau_min {tau_min:.2f} s"
    )


def noise_window_start(signal_start, length, rate, earliest_p_index, gap_s):
    """Equal-length noise window before both the signal gap and earliest P.

    Move the noise window earlier, not the signal or acceptance thresholds.
    Indices are on the common vertical UTC grid; stop is exclusive.
    """
    stop = min(signal_start - int(round(gap_s * rate)), int(np.floor(earliest_p_index)))
    return stop - length


def screened(event_time, screening, options):
    for item in screening:
        delay = event_time - item["time"]
        hours = (
            options["screening_hours_m6"]
            if (item["magnitude"] or 0) >= 6
            else options["screening_hours_m55"]
        )
        if 0 <= delay <= hours * 3600:
            return item["id"]
    return None


def overlapping_event(event, catalog, model, latitude, longitude, start, end):
    for other in catalog:
        if other["id"] == event["id"] or abs(other["time"] - event["time"]) > 900:
            continue
        try:
            geometry = event_geometry(model, other, latitude, longitude)
        except EventRejected:
            continue
        if start <= other["time"] + geometry.p_travel_s <= end:
            return other["id"]
    return None


def process_event(event, ctx):
    o = ctx.options
    p_epoch = epoch_for(ctx.p_epochs, event["time"], event["time"] + 600)
    v_epoch = epoch_for(ctx.v_epochs, event["time"], event["time"] + 600)
    station = v_epoch.station
    geometry = event_geometry(ctx.model, event, station.latitude, station.longitude)
    record = {
        "distance_km": round(geometry.distance_km, 3),
        "hypocentral_km": round(geometry.hypocentral_km, 3),
        "predicted_p_phase": geometry.p_phase,
        "predicted_p_utc": str(event["time"] + geometry.p_travel_s),
        "predicted_s_minus_p_s": round(geometry.s_travel_s - geometry.p_travel_s, 3),
        "q_s_per_km": round(geometry.q_s_per_m * 1000, 5),
        "eta_taup": round(geometry.eta, 5),
        "q_source": f"TauP {o['taup_model']} first P ray parameter",
    }
    ctx.current.update(record)
    if geometry.distance_km > o["max_distance_km"]:
        raise EventRejected("beyond_max_distance")
    if geometry.hypocentral_km < o["min_hypocentral_km"]:
        raise EventRejected("too_close_for_plane_wave")
    screening_hit = screened(event["time"], ctx.screening, o)
    if screening_hit:
        raise EventRejected("screening_catalogue_event", screening_hit)
    rate = ctx.rate
    search = search_half_width(geometry, event, o)
    lag_search = o["timing_lag_search_s"]
    pre = search + o["max_window_s"] + o["noise_gap_s"] + o["lead_s"] + 8.0
    post = search + o["max_window_s"] + 4.0
    quantum = ctx.quantum
    t_pred = event["time"] + geometry.p_travel_s
    seg_start = t_pred - quantize(pre, quantum, down=False)
    seg_len = quantize(pre + post, quantum, down=False)
    lag = 0.0
    lag_corr = None
    drift = o["timing_drift"]
    if drift and lag_search > 0:
        raise ValueError("Use either a lag search or a drift model, not both")
    if drift:
        if not (
            UTCDateTime(drift["valid_from_utc"])
            <= event["time"]
            <= UTCDateTime(drift["valid_to_utc"])
        ):
            raise EventRejected("outside_drift_model_range")
        days = (event["time"] - UTCDateTime(drift["reference_utc"])) / 86400
        lag = drift["intercept_s"] + drift["slope_s_per_day"] * days
        lag = round(lag * rate) / rate
    if lag_search > 0:
        wide = quantize(seg_len + 2 * lag_search, quantum, down=False)
        p_wide, v_wide, _, _ = load_pair(
            ctx.p_index,
            ctx.v_index,
            p_epoch,
            v_epoch,
            seg_start - lag_search,
            wide,
            rate,
            0.0,
            o,
        )
        lag, lag_corr = envelope_lag(p_wide, v_wide, rate, o["lag_band_hz"], lag_search)
        lag = round(lag * rate) / rate
    p_seg, v_seg, grid, _ = load_pair(
        ctx.p_index, ctx.v_index, p_epoch, v_epoch, seg_start, seg_len, rate, lag, o
    )
    i0 = int(round((t_pred - search - seg_start) * rate))
    # Stop the AIC search before the predicted S can dominate the variance change.
    s_minus_p = geometry.s_travel_s - geometry.p_travel_s
    i1 = int(round((t_pred + min(search, 0.8 * s_minus_p) - seg_start) * rate))
    band, snr_min = o["pick_band_hz"], o["pick_snr_min"]
    on_v, snr_v = measure_onset(v_seg, rate, i0 - int(4 * rate), i1, band, snr_min)
    on_p, snr_p = measure_onset(p_seg, rate, i0 - int(4 * rate), i1, band, snr_min)
    lag_coarse = lag
    if (lag_search > 0 or drift) and on_v is not None and on_p is not None:
        # Refine the coarse lag by the onset difference (bounded). The
        # sign of the pressure is untouched, so polarity stays testable; the
        # phase slope (timing) no longer is.
        shift = (on_p - on_v) / rate
        if abs(shift) > o["lag_refine_max_s"] or min(snr_v, snr_p) < snr_min:
            raise EventRejected(
                "lag_refinement_failed",
                f"onset shift {shift:+.2f} s; onset SNR Z {snr_v:.1f}, P {snr_p:.1f}",
            )
        lag = round((lag + shift) * rate) / rate
        p_seg, v_seg, grid, _ = load_pair(
            ctx.p_index, ctx.v_index, p_epoch, v_epoch, seg_start, seg_len, rate, lag, o
        )
        on_p, snr_p = measure_onset(p_seg, rate, i0 - int(4 * rate), i1, band, snr_min)
    ctx.current.update(
        {
            "grid": grid,
            "timing_lag_coarse_s": lag_coarse,
            "timing_lag_s": lag,
            "timing_lag_envelope_corr": lag_corr,
            "timing_mode": "drift_model"
            if drift
            else ("search" if lag_search else "none"),
        }
    )
    ctx.current.update(
        {
            "onset_vertical_utc": str(seg_start + on_v / rate)
            if on_v is not None
            else "",
            "onset_pressure_utc": str(seg_start + lag + on_p / rate)
            if on_p is not None
            else "",
            "onset_method": "auto_aic_causal_bandpass",
            "onset_snr_vertical": round(snr_v, 3),
            "onset_snr_pressure": round(snr_p, 3),
        }
    )
    if on_v is None or snr_v < o["pick_snr_min"]:
        raise EventRejected("vertical_onset_not_measured", f"SNR {snr_v:.2f}")
    if on_p is None or snr_p < o["pick_snr_min"]:
        raise EventRejected("pressure_onset_not_measured", f"SNR {snr_p:.2f}")
    delta = (on_p - on_v) / rate
    ctx.current["onset_difference_s"] = round(delta, 4)
    if abs(delta) > o["onset_tolerance_s"]:
        raise EventRejected("pressure_vertical_onset_mismatch", f"{delta:+.3f} s")
    onset_time = seg_start + on_v / rate
    regime, length, tau_min, tau_max = choose_window(
        geometry, float(onset_time - t_pred), ctx.depth, o
    )
    length = quantize(length, 1 / rate)
    n = int(round(length * rate))
    s0 = on_v - int(round(o["lead_s"] * rate))
    earliest_p = (t_pred - search - seg_start) * rate
    n0 = noise_window_start(s0, n, rate, earliest_p, o["noise_gap_s"])
    if n0 < 0 or s0 + n > len(v_seg):
        raise EventRejected("segment_too_short_for_windows")
    other = overlapping_event(
        event,
        ctx.catalog,
        ctx.model,
        v_epoch.station.latitude,
        v_epoch.station.longitude,
        seg_start + n0 / rate - o["overlap_before_s"],
        seg_start + (s0 + n) / rate,
    )
    if other:
        raise EventRejected("overlapping_catalogue_event", other)
    for name, data in (("pressure", p_seg), ("vertical", v_seg)):
        if flat_topped(data[s0 : s0 + n]):
            raise EventRejected("possible_clipping", name)
    ctx.current.update(
        {
            "regime": regime,
            "window_start_utc": str(seg_start + s0 / rate),
            "window_seconds": round(n / rate, 4),
            "noise_start_utc": str(seg_start + n0 / rate),
            "noise_end_utc": str(seg_start + (n0 + n) / rate),
            "earliest_predicted_p_utc": str(t_pred - search),
            "noise_guard_passed": n0 + n <= earliest_p,
            "round_trip_min_s": round(tau_min, 4),
            "round_trip_max_s": round(tau_max, 4),
            "max_abs_counts_pressure": float(np.max(np.abs(p_seg[s0 : s0 + n]))),
            "max_abs_counts_vertical": float(np.max(np.abs(v_seg[s0 : s0 + n]))),
        }
    )
    bands, arrays, diagnostics = analyse_event(
        p_seg,
        v_seg,
        slice(s0, s0 + n),
        slice(n0, n0 + n),
        rate,
        ctx.responses_for(p_epoch, v_epoch),
        regime,
        ctx.depth,
        geometry.eta,
        o,
    )
    ctx.current.update(diagnostics)
    ctx.current["bands"] = bands
    # Same string form as the master's channel_start_time / channel_end_time.
    ctx.current["p_epoch"] = [
        str(p_epoch.channel.start_date or ""),
        str(p_epoch.channel.end_date or ""),
    ]
    ctx.current["v_epoch"] = [
        str(v_epoch.channel.start_date or ""),
        str(v_epoch.channel.end_date or ""),
    ]
    usable = [b for b in bands if "gain_noise_corrected" in b]
    if not usable:
        raise EventRejected(
            "no_valid_band", json.dumps(diagnostics["bins_failing_in_bands"])
        )
    plot = {
        "rate": rate,
        "t0": float(seg_start - onset_time),
        "p": p_seg,
        "v": v_seg,
        "s0": s0,
        "n0": n0,
        "n": n,
        "tau": (tau_min, tau_max),
        "s_minus_p": geometry.s_travel_s - geometry.p_travel_s,
    }
    return arrays, plot


def station_summary(events, options, processed):
    accepted = [e for e in events if e["status"] == "accepted"]
    bands = []
    for lo, hi in options["bands_hz"]:
        values = [
            (e, b)
            for e in accepted
            for b in e["bands"]
            if b["band_hz"] == [lo, hi] and "gain_noise_corrected" in b
        ]
        item = {"band_hz": [lo, hi], "events": len(values)}
        if values:
            gains = np.array([b["gain_noise_corrected"] for _, b in values])
            median = float(np.median(gains))
            mad = float(np.median(np.abs(gains - median)))
            regimes = sorted({e["regime"] for e, _ in values})
            polarity = [b["polarity"] for _, b in values]
            item.update(
                {
                    "regimes": regimes,
                    "gain_median": median,
                    "gain_mad": mad,
                    "gain_relative_mad": mad / median if median else None,
                    "gain_h1_upper_median": float(
                        np.median([b["gain_h1_upper"] for _, b in values])
                    ),
                    "gain_h2_lower_median": float(
                        np.median([b["gain_h2_lower"] for _, b in values])
                    ),
                    "gain_model_low_median": float(
                        np.median([b["gain_model_low"] for _, b in values])
                    ),
                    "gain_model_high_median": float(
                        np.median([b["gain_model_high"] for _, b in values])
                    ),
                    "phase_deg_circular_mean": lp.circular_mean_deg(
                        [b["phase_deg"] for _, b in values]
                    ),
                    "polarity_counts": {str(k): polarity.count(k) for k in (1, -1, 0)},
                    "event_ids": [e["event_id"] for e, _ in values],
                }
            )
        bands.append(item)
    steady = [
        b
        for b in bands
        if b["events"] >= options["min_events"] and b.get("regimes") == ["steady_state"]
    ]
    any_band = [b for b in bands if b["events"]]
    if steady:
        outcome = "conditional_gain_candidate"
    elif any_band:
        outcome = "diagnostic_only"
    elif processed:
        outcome = "unresolved"
    else:
        outcome = "no_suitable_events"
    usable = [b for b in (steady or any_band) if b["events"] >= options["min_events"]]
    frequency = {"assessment": "insufficient_events"}
    if len(usable) == 1:
        frequency = {"assessment": "single_band"}
    elif len(usable) > 1:
        medians = np.array([b["gain_median"] for b in usable])
        spread = float(np.median([b["gain_relative_mad"] for b in usable]))
        ratio = float(medians.max() / medians.min())
        consistent = ratio <= 1 + 2 * spread
        frequency = {
            "assessment": "consistent_within_spread"
            if consistent
            else "frequency_dependent",
            "max_over_min_band_median": ratio,
            "median_relative_mad": spread,
        }
        if consistent:
            all_gains = [
                b2["gain_noise_corrected"]
                for e in accepted
                for b2 in e["bands"]
                if "gain_noise_corrected" in b2
                and any(b2["band_hz"] == b["band_hz"] for b in usable)
            ]
            frequency["scalar_median_all_bands"] = float(np.median(all_gains))
    polarity = {}
    for b in any_band:
        for key, value in b["polarity_counts"].items():
            polarity[key] = polarity.get(key, 0) + value
    votes = {k: v for k, v in polarity.items() if k != "0"}
    dominant = None
    if votes and max(votes.values()) > 0:
        winners = [k for k, v in votes.items() if v == max(votes.values())]
        if len(winners) == 1:
            dominant = winners[0]
    return {
        "outcome": outcome,
        "bands": bands,
        "frequency_dependence": frequency,
        "polarity_band_counts": polarity,
        "polarity_majority": int(dominant) if dominant else None,
        "polarity_summary_basis": "dominant classified band sign; not confidence",
        "polarity_classified_fraction": (
            sum(votes.values()) / sum(polarity.values())
            if sum(polarity.values())
            else None
        ),
        "gain_convention": "pressure multiplier: g * p_nominal = p_model",
        "score": None,
        "score_note": "No score assigned; scoring policy not agreed.",
    }


def plot_station(path, events, spectra, summary, title):
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    for key, arrays in spectra.items():
        m = arrays["valid"]
        if not m.any():
            continue
        f = arrays["frequency_hz"]
        axes[0].semilogx(f[m], arrays["gain_noise_corrected"][m], ".", ms=2, alpha=0.5)
        axes[1].semilogx(
            f[m], arrays["phase_obs_minus_pred_deg"][m], ".", ms=2, alpha=0.5
        )
        axes[2].semilogx(f[m], arrays["coherence"][m], ".", ms=2, alpha=0.5)
    for b in summary["bands"]:
        if b["events"]:
            lo, hi = b["band_hz"]
            axes[0].hlines(b["gain_median"], lo, hi, color="k", lw=2)
            axes[0].fill_between(
                [lo, hi],
                b["gain_median"] - b["gain_mad"],
                b["gain_median"] + b["gain_mad"],
                color="k",
                alpha=0.15,
            )
            axes[0].text(
                np.sqrt(lo * hi),
                b["gain_median"],
                f"n={b['events']}",
                ha="center",
                va="bottom",
            )
    axes[0].set_ylabel("Pressure multiplier g (valid bins)")
    axes[1].set_ylabel("Phase observed - predicted (deg)")
    axes[1].set_ylim(-180, 180)
    axes[2].set_ylabel("Squared coherence")
    axes[2].set_ylim(0, 1.02)
    axes[2].set_xlabel("Frequency (Hz)")
    for axis in axes:
        axis.grid(which="both", alpha=0.3)
    fig.suptitle(
        f"{title}\noutcome: {summary['outcome']} (no score; diagnostic figure)"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def plot_event(path, event, plot, depth, options):
    rate = plot["rate"]
    t = plot["t0"] + np.arange(len(plot["v"])) / rate
    band = options["pick_band_hz"]
    p = causal_band(plot["p"], rate, band)
    v = causal_band(plot["v"], rate, band)
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    axes[0].plot(t, p / np.max(np.abs(p)), lw=0.6, label="pressure counts (norm.)")
    axes[0].plot(
        t, v / np.max(np.abs(v)), lw=0.6, alpha=0.7, label="vertical counts (norm.)"
    )
    axes[0].legend(loc="upper right", fontsize=7)
    for axis in axes:
        s0, n0, n = plot["s0"], plot["n0"], plot["n"]
        axis.axvspan(t[s0], t[s0 + n - 1], color="g", alpha=0.15)
        axis.axvspan(t[n0], t[n0 + n - 1], color="grey", alpha=0.15)
        for k in range(1, 4):
            axis.axvline(k * plot["tau"][0], color="r", ls=":", lw=0.8)
        axis.axvline(plot["s_minus_p"], color="m", ls="--", lw=0.8)
    zoom = (t > -3) & (t < plot["s_minus_p"] + 5)
    axes[1].plot(t[zoom], p[zoom], lw=0.6)
    axes[1].set_xlabel(
        "Time after vertical onset (s); red: surface multiples; magenta: predicted S"
    )
    axes[0].set_xlim(t[plot["n0"]] - 2, plot["s_minus_p"] + 8)
    fig.suptitle(
        f"{event['event_id']} M{event['magnitude']} {event['distance_km']:.0f} km; "
        f"regime {event['regime']}; green signal, grey noise (filtered counts)"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def flatten(event):
    row = {
        k: v for k, v in event.items() if k not in {"bands", "bins_failing_in_bands"}
    }
    row["bins_failing_in_bands"] = json.dumps(event.get("bins_failing_in_bands", {}))
    for b in event.get("bands", []):
        tag = f"{b['band_hz'][0]:g}_{b['band_hz'][1]:g}hz"
        row[f"{tag}_valid_bins"] = b["valid_bins"]
        for key in (
            "gain_noise_corrected",
            "gain_h1_upper",
            "gain_h2_lower",
            "phase_deg",
            "polarity",
        ):
            row[f"{tag}_{key}"] = b.get(key, "")
    return row


def assert_fresh_output(path, protected):
    """Outputs go only to an empty folder under some results/LOCAL_P/ tree."""
    path = local_p_output(path)
    for root in protected:
        root = Path(root).expanduser().resolve()
        inside = path == root or root in path.parents
        local_p = any(
            parent.name == "LOCAL_P" and parent.parent.name == "results"
            for parent in [path, *path.parents]
            if root in parent.parents
        )
        if inside and (root.name != "results" or not local_p):
            raise ValueError(f"Output {path} lies inside protected input root {root}")
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {path}")
    return path


def run_station(
    *,
    pressure,
    vertical,
    inventory,
    waveform_dirs,
    catalog,
    results_dir,
    water_depth_m=None,
    water_depth_source=None,
    screening_catalog=None,
    start=None,
    end=None,
    options=None,
    protected=PROTECTED_DEFAULTS,
):
    o = {**DEFAULTS, **(options or {})}
    out = assert_fresh_output(results_dir, protected)
    out.mkdir(parents=True, exist_ok=True)
    started = datetime.now(UTC).isoformat()
    inputs = [Path(inventory), Path(catalog)] + (
        [Path(screening_catalog)] if screening_catalog else []
    )
    before = {str(p): {**file_state(p), "sha256": sha256(p)} for p in inputs}
    p_epochs = channel_epochs(inventory, pressure)
    v_epochs = channel_epochs(inventory, vertical)
    rates = sorted({float(e.channel.sample_rate) for e in p_epochs + v_epochs})
    rate = min(rates)
    summary = {
        "check": "calibration_check_pressure_vertical_local_p_waves",
        "protocol_version": PROTOCOL_VERSION,
        "pressure_seed_id": pressure,
        "vertical_seed_id": vertical,
        "started_utc": started,
        "code": git_state(),
        "options": o,
        "inputs": before,
        "method_note": "docs/local_p_method.md",
    }
    events_out, spectra, plots = [], {}, []
    regime_counts = {}
    try:
        if rate / 2 * o["nyquist_fraction"] < o["bands_hz"][0][1]:
            summary.update(
                {
                    "execution_status": "skipped",
                    "outcome": "not_applicable",
                    "reason": f"sample rate {rate} Hz below the local-P band",
                }
            )
            return summary
        epoch0 = v_epochs[0]
        depth = water_depth_m or water_depth(epoch0)
        if not depth:
            raise ValueError("No water depth from inventory; pass --water-depth")
        source = water_depth_source or (
            "explicit" if water_depth_m else "inventory elevation estimate"
        )
        catalog_events_all = load_local_catalog(catalog)
        screening = load_local_catalog(screening_catalog) if screening_catalog else []
        ctx = SimpleNamespace(
            options=o,
            p_epochs=p_epochs,
            v_epochs=v_epochs,
            p_index=DayFileIndex(waveform_dirs, pressure),
            v_index=DayFileIndex(waveform_dirs, vertical),
            model=TauPyModel(o["taup_model"]),
            catalog=catalog_events_all,
            screening=screening,
            rate=rate,
            quantum=duration_quantum(
                [rate] + [float(e.channel.sample_rate) for e in p_epochs + v_epochs]
            ),
            depth=float(depth),
            responses_for=responses_for,
            current={},
        )
        lat, lon = epoch0.station.latitude, epoch0.station.longitude
        days = sorted(set(ctx.p_index.by_day) & set(ctx.v_index.by_day))
        first = UTCDateTime(f"{days[0][:4]}-01-01") + (int(days[0][5:]) - 1) * 86400
        last = UTCDateTime(f"{days[-1][:4]}-01-01") + int(days[-1][5:]) * 86400
        start = max(first, UTCDateTime(start)) if start else first
        end = min(last, UTCDateTime(end)) if end else last
        candidates = []
        for event in catalog_events_all:
            if not start <= event["time"] < end - 600:
                continue
            if event["type"] and event["type"] != "earthquake":
                continue
            distance = (
                gps2dist_azimuth(event["latitude"], event["longitude"], lat, lon)[0]
                / 1000
            )
            if distance > o["max_distance_km"] or event["magnitude"] is None:
                continue
            # Deterministic priority: magnitude corrected for geometrical spreading.
            score = event["magnitude"] - 1.5 * np.log10(max(distance, 10.0) / 10.0)
            candidates.append((-score, str(event["time"]), event, distance))
        candidates.sort(key=lambda item: (item[0], item[1]))
        summary.update(
            {
                "station": {"latitude": lat, "longitude": lon},
                "water_depth_m": float(depth),
                "water_depth_source": source,
                "analysis_rate_hz": rate,
                "native_rates_hz": rates,
                "eta_min_p": lp.eta_min_for_p(o["sound_speed"], o["p_velocity_min"]),
                "first_resonance_hz_vertical": float(
                    lp.resonance_frequencies(depth)[0]
                ),
                "round_trip_vertical_s": lp.round_trip_time(depth),
                "coverage_utc": [str(start), str(end)],
                "catalogue_events": len(catalog_events_all),
                "candidates_within_distance": len(candidates),
                "candidates_processed": min(len(candidates), o["max_events"]),
            }
        )
        for _, _, event, distance in candidates[: o["max_events"]]:
            ctx.current = {
                "event_id": event["id"],
                "origin_utc": str(event["time"]),
                "latitude": event["latitude"],
                "longitude": event["longitude"],
                "depth_km": event["depth_km"],
                "magnitude": event["magnitude"],
                "magnitude_type": event["magnitude_type"],
                "catalogue_network": event["network"],
                "horizontal_error_km": event["horizontal_error_km"],
                "depth_error_km": event["depth_error_km"],
            }
            try:
                arrays, plot = process_event(event, ctx)
                ctx.current["status"] = "accepted"
                ctx.current["reason"] = ""
                spectra[event["id"]] = arrays
                plots.append((dict(ctx.current), plot))
                regime_counts[ctx.current["regime"]] = (
                    regime_counts.get(ctx.current["regime"], 0) + 1
                )
            except EventRejected as exc:
                ctx.current.update(
                    {"status": "rejected", "reason": exc.reason, "detail": exc.detail}
                )
            events_out.append(dict(ctx.current))
        for _, _, event, _ in candidates[o["max_events"] :]:
            events_out.append(
                {
                    "event_id": event["id"],
                    "origin_utc": str(event["time"]),
                    "magnitude": event["magnitude"],
                    "status": "not_processed",
                    "reason": "beyond_max_events_priority",
                }
            )
        processed = sum(e["status"] in {"accepted", "rejected"} for e in events_out)
        result = station_summary(events_out, o, processed)
        summary.update(result)
        reasons = {}
        for e in events_out:
            reasons[e["reason"] or "accepted"] = (
                reasons.get(e["reason"] or "accepted", 0) + 1
            )
        summary["event_reason_counts"] = reasons
        summary["regime_counts"] = regime_counts
        # Master keys: every pressure epoch overlapping the analysed coverage,
        # and the epochs that supplied accepted events.
        summary["pressure_epochs_in_coverage"] = [
            [str(e.channel.start_date or ""), str(e.channel.end_date or "")]
            for e in p_epochs
            if (e.channel.start_date is None or e.channel.start_date < end)
            and (e.channel.end_date is None or e.channel.end_date > start)
        ]
        summary["pressure_epochs_accepted"] = sorted(
            {tuple(e["p_epoch"]) for e in events_out if e.get("status") == "accepted"}
        )
        summary["epoch_results"] = []
        for epoch in summary["pressure_epochs_accepted"]:
            group = [
                e
                for e in events_out
                if tuple(e.get("p_epoch", [])) == epoch and e["status"] == "accepted"
            ]
            summary["epoch_results"].append(
                {
                    "pressure_epoch": list(epoch),
                    "vertical_epochs": sorted({tuple(e["v_epoch"]) for e in group}),
                    "regime_counts": {
                        r: sum(e["regime"] == r for e in group)
                        for r in {e["regime"] for e in group}
                    },
                    **station_summary(group, o, len(group)),
                }
            )
        summary["waveform_files_read"] = len(ctx.p_index.states) + len(
            ctx.v_index.states
        )
        summary["duplicate_layout_days_sds_preferred"] = sorted(
            ctx.p_index.duplicate_layout_days | ctx.v_index.duplicate_layout_days
        )
        changed = [
            path
            for index in (ctx.p_index, ctx.v_index)
            for path, state in index.states.items()
            if file_state(path) != state
        ]
        if changed:
            raise RuntimeError(f"Waveform files changed during the run: {changed[:3]}")
        summary["location_code_mismatch_in_filenames"] = (
            ctx.p_index.location_mismatch or ctx.v_index.location_mismatch
        )
        summary["timing_mode"] = (
            "drift_model"
            if o["timing_drift"]
            else ("search" if o["timing_lag_search_s"] > 0 else "none")
        )
        summary["timing_note"] = (
            "no lag correction; onsets measured independently on both channels"
            if summary["timing_mode"] == "none"
            else "lag estimated per event; phase slope is not an independent check"
        )
        summary["execution_status"] = "completed"
    except Exception as exc:  # recorded, then re-raised for the caller's log
        summary.update(
            {"execution_status": "failed", "error": f"{type(exc).__name__}: {exc}"}
        )
        raise
    finally:
        rows = [flatten(e) for e in events_out]
        if rows:
            fields = []
            for row in rows:
                fields.extend(k for k in row if k not in fields)
            with (out / "events.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
        if spectra:
            flat = {
                f"{key}__{name}": value
                for key, arrays in spectra.items()
                for name, value in arrays.items()
            }
            np.savez_compressed(out / "spectra.npz", **flat)
        if "bands" in summary:
            figure = out / "station_summary.png"
            plot_station(
                figure, events_out, spectra, summary, f"{pressure} / {vertical}"
            )
            summary["figure"] = str(figure)
            for event, plot in plots[: o["event_plots"]]:
                plot_event(
                    out / f"event_{event['event_id']}.png",
                    event,
                    plot,
                    summary["water_depth_m"],
                    o,
                )
        after = {str(p): {**file_state(p), "sha256": sha256(p)} for p in inputs}
        summary["inputs_unchanged"] = after == before
        summary["finished_utc"] = datetime.now(UTC).isoformat()
        summary["events_file"] = str(out / "events.csv")
        atomic_json(out / "summary.json", json.loads(json.dumps(summary, default=str)))
        if not summary["inputs_unchanged"]:
            raise RuntimeError("Protected input changed during the run")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pressure", required=True, help="Header SEED ID, e.g. XO.LT02..HDH"
    )
    parser.add_argument(
        "--vertical", required=True, help="Header SEED ID, e.g. XO.LT02..HHZ"
    )
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--waveform-dir", required=True, type=Path, nargs="+")
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--screening-catalog", type=Path)
    parser.add_argument("--water-depth", type=float)
    parser.add_argument("--water-depth-source")
    parser.add_argument("--start")
    parser.add_argument("--end")
    parser.add_argument("--max-events", type=int)
    parser.add_argument("--max-distance-km", type=float)
    parser.add_argument("--timing-lag-search", type=float, help="seconds; 0 disables")
    parser.add_argument("--options-json", type=Path, help="explicit option overrides")
    args = parser.parse_args(argv)
    options = json.loads(args.options_json.read_text()) if args.options_json else {}
    for key, value in (
        ("max_events", args.max_events),
        ("max_distance_km", args.max_distance_km),
        ("timing_lag_search_s", args.timing_lag_search),
    ):
        if value is not None:
            options[key] = value
    unknown = set(options) - set(DEFAULTS)
    if unknown:
        raise SystemExit(f"Unknown options: {sorted(unknown)}")
    summary = run_station(
        pressure=args.pressure,
        vertical=args.vertical,
        inventory=args.inventory,
        waveform_dirs=args.waveform_dir,
        catalog=args.catalog,
        results_dir=args.results_dir,
        water_depth_m=args.water_depth,
        water_depth_source=args.water_depth_source,
        screening_catalog=args.screening_catalog,
        start=args.start,
        end=args.end,
        options=options,
    )
    keys = (
        "execution_status",
        "outcome",
        "event_reason_counts",
        "regime_counts",
        "frequency_dependence",
        "polarity_majority",
    )
    print(json.dumps({k: summary.get(k) for k in keys}, indent=2, default=str))
    for band in summary.get("bands", []):
        if band["events"]:
            low, high = band["gain_h2_lower_median"], band["gain_h1_upper_median"]
            m_low, m_high = (
                band["gain_model_low_median"],
                band["gain_model_high_median"],
            )
            print(
                f"  {band['band_hz']} Hz: n={band['events']} "
                f"g={band['gain_median']:.3f} MAD={band['gain_mad']:.3f} "
                f"noise [{low:.3f}, {high:.3f}] model [{m_low:.3f}, {m_high:.3f}] "
                f"phase={band['phase_deg_circular_mean']:.0f} deg {band['regimes']}"
            )


if __name__ == "__main__":
    main()
