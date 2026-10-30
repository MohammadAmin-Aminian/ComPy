"""Microseism pressure–vertical consistency check on saved paired spectra.

Reads one saved paired-spectra run and one fixed Rayleigh reference result. It
never refits, applies or stores a pressure gain, and never writes beside inputs.

Theory: a homogeneous acoustic water layer with a free surface carries one
plane wave of horizontal slowness 1/c. The gauge (Lagrangian) pressure at the
seafloor is then ``P = rho_w * L(f, c) * a_up``, where the admittance height is
``L = H tan(y)/y`` (c > alpha) or ``H tanh(y)/y`` (c < alpha) with
``y = 2 pi f H sqrt(|1/alpha^2 - 1/c^2|)``; this is Zha & Webb (2016) A5, also
used by ``rayleigh_calibration.water_ratio``. Below the first water-column
resonance ``f_r = alpha/(4H)``, ``0 < L <= L(c=inf)`` for every real slowness,
so pressure polarity and an upper bound on the pressure multiplier follow
without a structure model. Above it, L spans (0, inf) and depends on the
unknown phase velocity, i.e. sediments and mode mixture (Ruan et al. 2014).
"""

import argparse
import hashlib
import json
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from compy_pressure.common import write_csv
from compy_pressure.phase import (
    PHASE_CODE_SHA256,
    PHASE_SUMMARY_METHOD,
    circular_median_degrees,
)
from compy_pressure.ppsd import atomic_json

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

PAIRED_FILES = (
    "comparison.json",
    "spectra.npz",
    "spectra.csv",
    "windows.csv",
    "spectra_coherence_transfer.png",
)
PAIRED_CONVENTIONS = {
    "cross_convention": "S_dp = mean(D * conjugate(P))",
    "transfer_direction": "pressure_to_upward_displacement",
    "transfer_units": "m/Pa",
    "pressure_psd_units": "Pa^2/Hz",
    "vertical_psd_units": "m^2/Hz",
}
VERTICAL_UNITS = {"M", "M/S", "M/S**2", "M/S^2", "M/S/S"}
STANDARD_GRAVITY = 9.80665
LITERATURE_BAND_HZ = (0.1, 0.2)  # Ruan et al. (2014) sediment-sensitive band
MODES = (0, 1, 2)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _water(depth, sound_speed):
    if not all(np.isfinite(v) and v > 0 for v in (depth, sound_speed)):
        raise ValueError("Water depth and sound speed must be positive and finite")


def admittance_height(f, phase_velocity, depth, sound_speed, gravity=0.0):
    """Seafloor gauge P / (rho_w * upward acceleration), in metres.

    ``phase_velocity`` may be ``inf`` (vertical incidence). ``gravity`` adds the
    first-order free-surface term for the Lagrangian gauge pressure; interior
    compressibility-gravity coupling is not modelled. No branch is rejected.
    """
    _water(depth, sound_speed)
    f, c = np.broadcast_arrays(
        np.asarray(f, dtype=float), np.asarray(phase_velocity, dtype=float)
    )
    if np.any(~np.isfinite(f) | (f <= 0)) or np.any(np.isnan(c) | (c <= 0)):
        raise ValueError("Frequency and phase velocity must be positive")
    if not np.isfinite(gravity) or gravity < 0:
        raise ValueError("Gravity must be finite and non-negative")
    omega = 2 * np.pi * f
    q = 1 / c**2 - 1 / sound_speed**2  # nu^2 / omega^2
    y = omega * depth * np.sqrt(abs(q))
    e = gravity * np.sqrt(abs(q)) / omega  # |g nu / omega^2|
    evanescent = q > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(evanescent, np.tanh(y), np.tan(y))
        shape = np.where(y < 1e-6, 1 + np.sign(-q) * y**2 / 3, t / y)
        surface = np.where(
            evanescent, (1 - e**2) / (1 - e * t), (1 + e**2) / (1 + e * t)
        )
    return depth * shape * surface


def first_resonance_hz(depth, sound_speed):
    """Lowest frequency at which vertical incidence reaches tan(y) = inf."""
    _water(depth, sound_speed)
    return sound_speed / (4 * depth)


def vertical_incidence_limit(f, depth, sound_speed):
    """Largest first-branch height; inf at or above the first resonance."""
    f = np.asarray(f, dtype=float)
    below = f < first_resonance_hz(depth, sound_speed)
    limit = np.full(f.shape, np.inf)
    limit[below] = admittance_height(f[below], np.inf, depth, sound_speed)
    return limit


def phase_velocity_from_height(f, height, depth, sound_speed):
    """Invert positive heights on the first (lowest-velocity) branch, no gravity.

    Returns NaN where the height is non-positive, implies c < 1 m/s, or exceeds
    the vertical-incidence limit below the first water-column resonance.
    """
    _water(depth, sound_speed)
    f, height = np.broadcast_arrays(
        np.asarray(f, dtype=float), np.asarray(height, dtype=float)
    )
    result = np.full(f.shape, np.nan)
    for index in np.ndindex(f.shape):
        frequency, target = f[index], height[index] / depth
        if not (np.isfinite(target) and target > 0 and frequency > 0):
            continue
        scale = 2 * np.pi * frequency * depth / sound_speed

        def residual(u):  # u = alpha^2/c^2 - 1; first branch only
            y = scale * np.sqrt(abs(u))
            if y < 1e-6:
                return 1 + np.sign(-u) * y**2 / 3 - target
            return (np.tanh(y) if u > 0 else np.tan(y)) / y - target

        resonance = -((np.pi / (2 * scale)) ** 2)
        low = resonance * (1 - 1e-12) if resonance > -1 else -1.0
        high = (sound_speed / 1.0) ** 2 - 1  # c = 1 m/s
        at_low = residual(low)
        if low == -1.0 and abs(at_low) <= 1e-12 * target:
            result[index] = np.inf  # exactly vertical incidence
            continue
        if not at_low > 0 > residual(high):
            continue
        u = brentq(residual, low, high, xtol=1e-15, rtol=1e-13, maxiter=500)
        result[index] = np.inf if u <= -1 else sound_speed / np.sqrt(1 + u)
    return result


def load_layers(path):
    layers = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    if layers.shape[1] != 4 or not np.isfinite(layers).all() or layers[-1, 0] != 0:
        raise ValueError("Expected km, km/s, km/s, g/cm3 layers ending in half-space")
    return layers


def modal_phase_velocity(layers, frequency, modes=MODES, water=None):
    """disba fundamental/overtone Rayleigh phase velocity in m/s; NaN if absent.

    ``water=(depth_m, sound_speed_m_s, density_kg_m3)`` prepends a fluid layer.
    """
    from disba import PhaseDispersion

    model = np.asarray(layers, dtype=float)
    if water is not None:
        depth, speed, density = water
        _water(depth, speed)
        model = np.vstack([[depth / 1e3, speed / 1e3, 0.0, density / 1e3], model])
    frequency = np.asarray(frequency, dtype=float)
    periods = np.sort(1 / frequency)
    solver = PhaseDispersion(*model.T)
    result = np.full((len(modes), len(frequency)), np.nan)
    for row, mode in enumerate(modes):
        solved = solver(periods, mode=mode, wave="rayleigh")
        found = 1 / np.asarray(solved.period)
        speed = np.asarray(solved.velocity) * 1e3
        for column, value in enumerate(frequency):
            match = np.isclose(found, value, rtol=1e-10, atol=0)
            if match.any():
                result[row, column] = speed[np.argmax(match)]
    return result


def check_paired_run(run_directory):
    """Return validated saved arrays; refuse convention or arithmetic drift."""
    run_directory = Path(run_directory)
    summary = json.loads((run_directory / "comparison.json").read_text())
    for key, expected in PAIRED_CONVENTIONS.items():
        if summary.get(key) != expected:
            raise ValueError(f"Unexpected paired convention {key}: {summary.get(key)}")
    units = {k: str(v).upper() for k, v in summary.get("input_units", {}).items()}
    if units.get("pressure") != "PA" or units.get("vertical") not in VERTICAL_UNITS:
        raise ValueError("Paired run must convert Pa pressure and a seismic vertical")
    with np.load(run_directory / "spectra.npz") as saved:
        data = {key: saved[key] for key in saved.files}
    f = data["frequency_hz"]
    if np.any(np.diff(f) <= 0) or not data["valid"].all():
        raise ValueError("Saved frequencies must be increasing and all bins valid")
    count = len(data["retained_window_starts"])
    if count != summary.get("retained_segments") or count < 2:
        raise ValueError("Retained window count disagrees with saved summary")
    for mean, windows in (
        ("pressure_psd", "pressure_window_psd"),
        ("displacement_psd", "displacement_window_psd"),
        ("cross_dp", "cross_window_dp"),
    ):
        if not np.allclose(data[mean], data[windows].mean(axis=0), rtol=1e-12, atol=0):
            raise ValueError(f"Saved {mean} is not the mean of its windows")
    pp, dd, dp = data["pressure_psd"], data["displacement_psd"], data["cross_dp"]
    if not np.allclose(data["transfer_dp"], dp / pp, rtol=1e-12, atol=0):
        raise ValueError("Saved transfer is not S_dp / S_pp")
    if not np.allclose(
        data["coherence_squared"], abs(dp) ** 2 / (pp * dd), rtol=1e-10, atol=0
    ):
        raise ValueError("Saved coherence is not |S_dp|^2 / (S_pp S_dd)")
    return summary, data


def transfer_estimates(pp, dd, dp):
    """H1 = S_dp/S_pp (pressure noise biases low); H2 = S_dd/S_pd (biased high)."""
    return dp / pp, dd / np.conj(dp)


def implied_multiplier(f, height, transfer, density):
    """Pressure multiplier g for which theory height reproduces |transfer|."""
    return (2 * np.pi * f) ** 2 * density * height * abs(transfer)


def _median(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.median(values)) if len(values) else None


def _ratio_pct(value, reference):
    return None if value is None or reference is None else 100 * (value / reference - 1)


def daily_transfer(data, mask):
    """Per-UTC-day H1/H2 from saved hourly windows (descriptive stability)."""
    days = np.array([start[:10] for start in data["retained_window_starts"]])
    result = {}
    for day in sorted(set(days)):
        rows = days == day
        pp = data["pressure_window_psd"][rows][:, mask].mean(axis=0)
        dd = data["displacement_window_psd"][rows][:, mask].mean(axis=0)
        dp = data["cross_window_dp"][rows][:, mask].mean(axis=0)
        result[day] = (int(rows.sum()), *transfer_estimates(pp, dd, dp))
    return result


def run_check(
    paired_run,
    reference_result,
    layers_file,
    output_root,
    *,
    water=None,
    reference_sha256=None,
    frequency_min=0.02,
    frequency_max=0.35,
    illustrative_density=1035.0,
    illustrative_depth_fraction=0.02,
    sound_speed_range=(1480.0, 1540.0),
    show=False,
):
    """Diagnose one saved pair against a fixed reference result, or, without
    one, against ``water=(depth_m, density_kg_m3, sound_speed_m_s)`` only."""
    paired_run, layers_file = Path(paired_run), Path(layers_file)
    output_root = Path(output_root)
    if (reference_result is None) == (water is None):
        raise ValueError("Supply exactly one reference result or water settings")
    if output_root.resolve().is_relative_to(paired_run.resolve().parent):
        raise ValueError("Output must not be written beside the saved paired runs")
    inputs = [
        paired_run / name for name in PAIRED_FILES if (paired_run / name).exists()
    ]
    if reference_result is not None:
        reference_result = Path(reference_result)
        inputs.append(reference_result)
    elif reference_sha256:
        raise ValueError("A reference hash needs a reference result")
    inputs.append(layers_file)
    before = {
        str(p): (sha256(p), p.stat().st_size, p.stat().st_mtime_ns) for p in inputs
    }
    if reference_sha256 and before[str(reference_result)][0] != reference_sha256:
        raise ValueError("Reference result hash differs from the recorded hash")
    summary, data = check_paired_run(paired_run)
    station = summary["pressure_epoch"]["Seed_id"].rsplit(".", 2)[0]
    if reference_result is None:
        reference, g_ref, sign = None, None, None
        depth, density, alpha = map(float, water)
        _water(depth, alpha)
        water_source = "caller (master Profondeur); unverified bathymetry"
    else:
        reference = json.loads(reference_result.read_text())
        settings = reference["settings"]
        g_ref = float(reference["median_pressure_multiplier"])
        sign = int(settings["pressure_sign"])
        depth = float(settings["depth_m"])
        density = float(settings["density_kg_m3"])
        alpha = float(settings["sound_speed_m_s"])
        water_source = "reference settings (inventory elevation); unverified bathymetry"
        if reference.get("station") != station:
            raise ValueError("Reference and paired run belong to different stations")
        epochs = reference.get("epoch_keys", [])
        if any(summary[k] not in epochs for k in ("pressure_epoch", "vertical_epoch")):
            raise ValueError("Reference and paired run use different channel epochs")

    f_all = data["frequency_hz"]
    keep = (f_all >= frequency_min) & (f_all <= frequency_max)
    f = f_all[keep]
    pp, dd, dp = (
        data[k][keep] for k in ("pressure_psd", "displacement_psd", "cross_dp")
    )
    coherence = data["coherence_squared"][keep]
    h1, h2 = transfer_estimates(pp, dd, dp)
    phase = np.angle(h1, deg=True)
    # First branch (L > 0): sign -1 gives 0 deg, +1 gives 180 deg below f_r.
    expected_phase = None if sign is None else (0.0 if sign == -1 else 180.0)
    omega2 = (2 * np.pi * f) ** 2

    f_res = first_resonance_hz(depth, alpha)
    f_ig = np.sqrt(STANDARD_GRAVITY / (2 * np.pi * depth))
    bands = {
        "long_wave": (1.5 * f_ig, 0.9 * f_res),
        "literature_0p1_0p2": LITERATURE_BAND_HZ,
        "above_resonance": (f_res, 2 * f_res),
    }

    layers = load_layers(layers_file)
    water_model = modal_phase_velocity(layers, f, water=(depth, alpha, density))
    dry_model = modal_phase_velocity(layers, f, modes=(0,))[0]
    curves = {f"mode{m}_water": water_model[i] for i, m in enumerate(MODES)}
    curves["mode0_no_water"] = dry_model
    heights = {
        name: np.where(
            np.isfinite(c),
            admittance_height(f, np.where(np.isfinite(c), c, 1.0), depth, alpha),
            np.nan,
        )
        for name, c in curves.items()
    }
    limit = vertical_incidence_limit(f, depth, alpha)
    observed_height = (
        {}
        if g_ref is None
        else {
            k: g_ref / (omega2 * density * abs(h)) for k, h in (("h1", h1), ("h2", h2))
        }
    )
    apparent_c = {
        k: phase_velocity_from_height(f, v, depth, alpha)
        for k, v in observed_height.items()
    }
    g_upper = {
        k: implied_multiplier(f, limit, h, density) for k, h in (("h1", h1), ("h2", h2))
    }
    g_lower_alpha = implied_multiplier(f, depth, h1, density)
    g_model = {
        f"{name}_{k}": implied_multiplier(f, heights[name], h, density)
        for name in heights
        for k, h in (("h1", h1), ("h2", h2))
    }

    def in_band(name):
        low, high = bands[name]
        return (f >= low) & (f <= high)

    results = {}
    for name in bands:
        m = in_band(name)
        if not m.any():
            continue
        entry = {
            "frequency_hz": list(bands[name]),
            "bins": int(m.sum()),
            "median_coherence_squared": _median(coherence[m]),
            "median_phase_deg": circular_median_degrees(phase[m]),
            "phase_summary_method": PHASE_SUMMARY_METHOD,
            "fraction_phase_within_90deg_of_zero": float(
                np.mean(np.cos(np.radians(phase[m])) > 0)
            ),
            "median_reference_height_over_depth": {
                k: _median(v[m] / depth) for k, v in observed_height.items()
            },
            "median_apparent_phase_velocity_m_s": {
                k: _median(v[m]) for k, v in apparent_c.items()
            },
            "apparent_phase_velocity_defined_fraction": {
                k: float(np.mean(np.isfinite(v[m]))) for k, v in apparent_c.items()
            },
            "median_model_phase_velocity_m_s": {
                k: _median(v[m]) for k, v in curves.items()
            },
            "model_conditional_multiplier": {
                k: {
                    "median": _median(v[m]),
                    "vs_reference_pct": _ratio_pct(_median(v[m]), g_ref),
                }
                for k, v in g_model.items()
            },
        }
        if expected_phase is not None:
            entry["fraction_phase_within_90deg_of_expected"] = float(
                np.mean(np.cos(np.radians(phase[m] - expected_phase)) > 0)
            )
        below = m & np.isfinite(limit)
        if below.any() and g_ref is None:
            entry["structure_free_upper_bound"] = {
                "bins": int(below.sum()),
                "median_upper_bound_h1": _median(g_upper["h1"][below]),
                "median_upper_bound_h2": _median(g_upper["h2"][below]),
            }
            entry["lower_bound_if_c_at_least_water_speed"] = {
                "median_h1": _median(g_lower_alpha[below]),
            }
        elif below.any():
            lo_ok = g_ref <= g_upper["h1"][below]
            hi_bad = g_ref > g_upper["h2"][below]
            entry["structure_free_upper_bound"] = {
                "bins": int(below.sum()),
                "median_upper_bound_h1": _median(g_upper["h1"][below]),
                "median_upper_bound_h2": _median(g_upper["h2"][below]),
                "fraction_reference_consistent": float(np.mean(lo_ok)),
                "fraction_reference_violates": float(np.mean(hi_bad)),
                "fraction_inconclusive": float(np.mean(~lo_ok & ~hi_bad)),
            }
            entry["lower_bound_if_c_at_least_water_speed"] = {
                "median_h1": _median(g_lower_alpha[below]),
                "fraction_reference_above": float(
                    np.mean(g_ref >= g_lower_alpha[below])
                ),
            }
        results[name] = entry

    long_wave = in_band("long_wave")
    sensitivity = {}
    perturbed = {
        "depth_minus": dict(depth=depth * (1 - illustrative_depth_fraction)),
        "depth_plus": dict(depth=depth * (1 + illustrative_depth_fraction)),
        "sound_speed_low": dict(alpha=sound_speed_range[0]),
        "sound_speed_high": dict(alpha=sound_speed_range[1]),
        "column_density": dict(density=illustrative_density),
        "free_surface_gravity": dict(gravity=STANDARD_GRAVITY),
    }
    c0 = curves["mode0_water"][long_wave]
    for label, change in perturbed.items():
        d_, a_ = change.get("depth", depth), change.get("alpha", alpha)
        rho_ = change.get("density", density)
        height = admittance_height(f[long_wave], c0, d_, a_, change.get("gravity", 0.0))
        value = _median(implied_multiplier(f[long_wave], height, h1[long_wave], rho_))
        sensitivity[label] = {
            "setting": change,
            "median_h1": value,
            "vs_reference_pct": _ratio_pct(value, g_ref),
        }

    # Fixed equal log-width split of the long-wave band; not chosen from the data.
    edges = np.geomspace(*bands["long_wave"], 5)
    days = daily_transfer(data, keep)
    sub_bands = []
    for low, high in zip(edges[:-1], edges[1:]):
        m = (f >= low) & (f < high)
        height = heights["mode0_water"][m]
        daily = {
            day: [
                _median(implied_multiplier(f[m], height, d[m], density))
                for d in (d1, d2)
            ]
            for day, (_, d1, d2) in days.items()
        }
        sub_bands.append(
            {
                "frequency_hz": [float(low), float(high)],
                "bins": int(m.sum()),
                "median_coherence_squared": _median(coherence[m]),
                "median_phase_deg": circular_median_degrees(phase[m]),
                "phase_summary_method": PHASE_SUMMARY_METHOD,
                "mode0_water_h1": _median(g_model["mode0_water_h1"][m]),
                "mode0_water_h2": _median(g_model["mode0_water_h2"][m]),
                "upper_bound_h1": _median(g_upper["h1"][m]),
                "upper_bound_h2": _median(g_upper["h2"][m]),
                "daily_mode0_water_h1_h2": daily,
            }
        )
    window_counts = {day: value[0] for day, value in days.items()}

    response = data["pressure_response"][keep][long_wave]
    rows = [
        {
            "frequency_hz": f"{f[i]:.9g}",
            "coherence_squared": f"{coherence[i]:.6g}",
            "h1_m_per_pa": f"{abs(h1[i]):.9g}",
            "h2_m_per_pa": f"{abs(h2[i]):.9g}",
            "phase_deg": f"{phase[i]:.4f}",
            **{
                f"reference_height_over_depth_{k}": f"{v[i] / depth:.6g}"
                for k, v in observed_height.items()
            },
            **{
                f"apparent_phase_velocity_{k}_m_s": f"{v[i]:.6g}"
                for k, v in apparent_c.items()
            },
            **{f"c_{k}_m_s": f"{v[i]:.6g}" for k, v in curves.items()},
            "upper_bound_h1": f"{g_upper['h1'][i]:.6g}",
            "upper_bound_h2": f"{g_upper['h2'][i]:.6g}",
            "lower_bound_c_ge_alpha_h1": f"{g_lower_alpha[i]:.6g}",
            **{f"g_{k}": f"{v[i]:.6g}" for k, v in g_model.items()},
        }
        for i in range(len(f))
    ]

    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    output = output_root / station / f"run_{stamp}_{uuid4().hex[:8]}"
    output.mkdir(parents=True, exist_ok=False)
    files = {
        "summary": output / "microseism_check.json",
        "bins": output / "microseism_check_bins.csv",
        "plot": output / "microseism_check.png",
    }
    write_csv(files["bins"], rows)
    plot_check(
        f,
        coherence,
        phase,
        expected_phase,
        observed_height,
        heights,
        limit,
        g_upper,
        g_lower_alpha,
        g_model,
        apparent_c,
        curves,
        g_ref,
        depth,
        bands,
        f_res,
        f_ig,
        f"{station} microseism pressure–vertical check (diagnostic)",
        files["plot"],
        show,
    )

    after = {
        str(p): (sha256(p), p.stat().st_size, p.stat().st_mtime_ns) for p in inputs
    }
    if after != before:
        raise RuntimeError("An input changed during the microseism check")
    result = {
        "status": "diagnostic_only",
        "station": station,
        "pressure_seed_id": summary["pressure_epoch"]["Seed_id"],
        "vertical_seed_id": summary["vertical_seed_id"],
        "interval_utc": [summary["start_utc"], summary["end_utc"]],
        "retained_windows": summary["retained_segments"],
        "identifiability": (
            "The P/Z ratio constrains g/g_z jointly with the horizontal slowness of "
            "the coherent wavefield. Below alpha/(4H) it fixes polarity and an upper "
            "bound on g without a structure model; a two-sided interval needs a "
            "minimum phase velocity, and a value needs a dispersion/mode model. "
            "Above alpha/(4H) the admissible range is unbounded."
        ),
        "reference": None
        if reference is None
        else {
            "path": str(reference_result),
            "sha256": before[str(reference_result)][0],
            "pressure_multiplier": g_ref,
            "pressure_sign": sign,
            "pressure_sign_status": reference["pressure_sign_assumption"]["status"],
            "role": "fixed comparison reference; not refitted or applied",
        },
        "water": {
            "depth_m": depth,
            "density_kg_m3": density,
            "sound_speed_m_s": alpha,
            "first_resonance_hz": f_res,
            "infragravity_cutoff_hz": float(f_ig),
            "source": water_source,
        },
        "model": {
            "layers_file": str(layers_file),
            "layers_sha256": before[str(layers_file)][0],
            "water_layer_added": True,
            "disba_version": version("disba"),
            "status": "historical model; not a validated local Earth model",
        },
        "estimators": {
            "h1": "S_dp/S_pp; uncorrelated pressure noise biases |H| low",
            "h2": "S_dd/conj(S_dp); uncorrelated vertical noise biases |H| high",
            "multiplier": "g = (2 pi f)^2 rho_w L |H|; phase checked separately",
        },
        "expected_phase_deg": expected_phase,
        "bands": results,
        "long_wave_sensitivity_mode0_water": sensitivity,
        "long_wave_sub_bands": sub_bands,
        "daily_window_counts": window_counts,
        "nominal_pressure_response_long_wave_amplitude_range": float(
            abs(response).max() / abs(response).min()
        ),
        "inputs": {
            k: {"sha256": v[0], "size": v[1], "mtime_ns": v[2]}
            for k, v in before.items()
        },
        "inputs_unchanged": True,
        "code_sha256": sha256(__file__),
        "phase_summary_method": PHASE_SUMMARY_METHOD,
        "phase_summary_code_sha256": PHASE_CODE_SHA256,
        "numpy_version": np.__version__,
        "files": {k: str(v) for k, v in files.items()},
        "gain_applied": False,
        "master_or_scores_changed": False,
        "limitations": [
            "Vertical response assumed correct: results constrain g/g_z",
            "Depth from inventory elevation; density/sound speed homogeneous",
            "Single horizontally stratified water layer; plane waves, real slowness",
            "Model dispersion is the historical RR52 layer model, not validated",
            "Hourly windows over one week are correlated; no confidence interval",
            "No tilt correction; incoherent tilt noise is bracketed by H1/H2 only",
            "Polarity is inferred from phase; hardware polarity is not verified",
        ],
    }
    atomic_json(files["summary"], result)
    return result


def plot_check(
    f,
    coherence,
    phase,
    expected_phase,
    observed_height,
    heights,
    limit,
    g_upper,
    g_lower_alpha,
    g_model,
    apparent_c,
    curves,
    g_ref,
    depth,
    bands,
    f_res,
    f_ig,
    title,
    destination,
    show,
):
    figure, axes = plt.subplots(3, 2, figsize=(13, 12), sharex=True)
    axes = axes.ravel()
    for axis in axes[:5]:
        axis.set_xscale("log")
        low, high = bands["long_wave"]
        axis.axvspan(low, high, color="tab:green", alpha=0.08)
        axis.axvspan(*LITERATURE_BAND_HZ, color="tab:orange", alpha=0.08)
        axis.axvline(f_res, color="k", ls=":", lw=0.8)
        axis.grid(True, which="both", alpha=0.25)
    axes[0].plot(f, coherence, lw=0.6)
    axes[0].set(
        ylabel="coherence²", ylim=(0, 1), title="Saved coherence (no threshold)"
    )
    axes[1].plot(f, phase, ".", ms=1.5)
    axes[1].axhline(0, color="r", lw=0.8)  # sign -1 below f_r
    for edge in (-180, 180):
        axes[1].axhline(edge, color="gray", ls="--", lw=0.8)  # sign +1
    axes[1].set(
        ylabel="arg H_dp (deg)", ylim=(-180, 180), title="Phase vs expected polarity"
    )
    if observed_height:
        axes[2].fill_between(
            f,
            observed_height["h2"] / depth,
            observed_height["h1"] / depth,
            color="k",
            alpha=0.35,
            label=f"data × g_ref={g_ref:.4f} (H2–H1)",
        )
    for name, style in (
        ("mode0_water", "b-"),
        ("mode1_water", "m-"),
        ("mode0_no_water", "c--"),
    ):
        axes[2].plot(f, heights[name] / depth, style, lw=1, label=name)
    axes[2].plot(f, limit / depth, "r:", lw=1, label="vertical incidence (bound)")
    axes[2].axhline(1, color="g", lw=0.8, label="c = water sound speed")
    axes[2].set(ylabel="L / H", ylim=(0, 6), title="Admittance height P/(ρ a_z H)")
    axes[2].legend(fontsize=7)
    norm = 1.0 if g_ref is None else g_ref
    shown = f >= bands["long_wave"][0]  # below it, compliance/tilt/noise dominate
    g_model = {k: np.where(shown, v, np.nan) for k, v in g_model.items()}
    g_upper = {k: np.where(shown, v, np.nan) for k, v in g_upper.items()}
    g_lower_alpha = np.where(shown, g_lower_alpha, np.nan)
    axes[3].fill_between(
        f,
        g_model["mode0_water_h1"] / norm,
        g_model["mode0_water_h2"] / norm,
        color="b",
        alpha=0.3,
        label="mode0 water (H1–H2)",
    )
    axes[3].plot(
        f,
        g_model["mode0_no_water_h1"] / norm,
        "c-",
        lw=0.6,
        label="mode0 no water (H1)",
    )
    axes[3].plot(f, g_upper["h2"] / norm, "r-", lw=0.6, label="upper bound (H2, c=∞)")
    axes[3].plot(f, g_lower_alpha / norm, "g-", lw=0.6, label="lower bound if c≥α (H1)")
    axes[3].axhline(1, color="k", lw=0.8)
    axes[3].set(
        ylabel="g" if g_ref is None else "g / g_ref",
        ylim=(0.3, 2.0) if g_ref is None else (0.5, 1.5),
        title="Implied multiplier (diagnostic)",
    )
    axes[3].legend(fontsize=7)
    for name, style in (
        ("mode0_water", "b-"),
        ("mode1_water", "m-"),
        ("mode2_water", "y-"),
        ("mode0_no_water", "c--"),
    ):
        axes[4].plot(f, curves[name], style, lw=1, label=name)
    if apparent_c:
        axes[4].plot(f, apparent_c["h1"], "k.", ms=1.5, label="apparent (g_ref, H1)")
        axes[4].plot(
            f, apparent_c["h2"], ".", color="gray", ms=1.5, label="apparent (g_ref, H2)"
        )
    axes[4].set(
        ylabel="phase velocity (m/s)",
        ylim=(0, 6000),
        xlabel="frequency (Hz)",
        title="First-branch apparent vs model velocity",
    )
    axes[4].legend(fontsize=7)
    axes[5].axis("off")
    axes[5].text(
        0.02,
        0.98,
        "Green: long-wave band 1.5 f_IG–0.9 f_r\n"
        "Grey band: data with g_ref between H1 and H2\n"
        "Phase: red 0° = sign −1, grey ±180° = sign +1\n"
        "Orange: Ruan et al. 0.1–0.2 Hz band\n"
        f"Dotted: first water resonance f_r={f_res:.4f} Hz\n"
        f"f_IG={f_ig:.4f} Hz\n\n"
        "No gain fitted or applied; model curves use the\n"
        "historical RR52 layers, not a validated model.",
        va="top",
        fontsize=9,
        transform=axes[5].transAxes,
    )
    axes[3].set_xlabel("frequency (Hz)")
    axes[3].tick_params(labelbottom=True)
    figure.suptitle(title)
    figure.tight_layout()
    figure.savefig(destination, dpi=130)
    if show:
        plt.show()
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--paired-run", type=Path, required=True)
    parser.add_argument("--reference-result", type=Path, required=True)
    parser.add_argument("--reference-sha256")
    parser.add_argument(
        "--layers-file",
        type=Path,
        default=Path(__file__).with_name("reference_data")
        / "rr52_historical_rayleigh_layers.csv",
    )
    parser.add_argument(
        "--output-root", type=Path, default=Path("results/EXPERIMENTS/MICROSEISM")
    )
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    result = run_check(
        args.paired_run,
        args.reference_result,
        args.layers_file,
        args.output_root,
        reference_sha256=args.reference_sha256,
        show=args.show,
    )
    print(f"{result['status']}: {result['files']['summary']}")


if __name__ == "__main__":
    main()
