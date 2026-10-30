"""Opt-in scientific trial helpers; never change production responses or defaults."""

import numpy as np
from scipy.optimize import least_squares

from compy_pressure.phase import PHASE_SUMMARY_METHOD, circular_median_degrees
from compy_pressure.rayleigh_calibration import fit_multiplier


def real_pressure_poles(response):
    """Select positive real poles only in Laplace radians/second stages."""
    poles = []
    for stage in response.response_stages:
        positive = [p for p in getattr(stage, "poles", []) if p.real > 0]
        if not positive:
            continue
        if (
            getattr(stage, "pz_transfer_function_type", "")
            != "LAPLACE (RADIANS/SECOND)"
        ):
            raise ValueError("Positive pressure poles require explicit rad/s units")
        if any(abs(p.imag) > 1e-12 for p in positive):
            raise ValueError("Complex positive poles are outside this RC experiment")
        poles.extend(float(p.real) for p in positive)
    return poles


def pole_factor(frequency, positive_poles):
    """R_nominal/R_alternative, so P_alternative=P_nominal*factor.

    Reflect +a to -a without changing gain or normalization. This all-pass
    factor changes pressure phase only, not PSD or coherence. No XML write.
    """
    f = np.asarray(frequency, dtype=float)
    if not np.isfinite(f).all() or np.any(f <= 0):
        raise ValueError("Positive finite frequencies required")
    factor = np.ones(f.shape, dtype=complex)
    for a in positive_poles:
        if not np.isfinite(a) or a <= 0:
            raise ValueError("Positive finite real poles required")
        factor *= (2j * np.pi * f + a) / (2j * np.pi * f - a)
    return factor


def apply_response_hypothesis(spectra, pressure_ratio, vertical_ratio):
    """Saved event spectra as they would read under alternative responses.

    Ratios are R_nominal/R_alternative on the saved grid, so P_alt = P*rp and
    A_alt = A*rz. Deterministic per-bin factors leave coherence unchanged; the
    cross-spectrum gains rp*conj(rz) and the |P|/|A| ratio |rp|/|rz|. Returns
    only the fields `evaluate` reads, as new arrays; saved arrays and response
    files are untouched. Experiment-only: no production default uses this.
    """
    f = np.asarray(spectra["frequency_hz"], dtype=float)
    rp = np.broadcast_to(np.asarray(pressure_ratio, dtype=complex), f.shape)
    rz = np.broadcast_to(np.asarray(vertical_ratio, dtype=complex), f.shape)
    if not (np.isfinite(rp).all() and np.isfinite(rz).all()) or np.any(
        (rp == 0) | (rz == 0)
    ):
        raise ValueError("Response hypothesis ratios must be finite and nonzero")
    cross = np.asarray(spectra["cross_pa"]) * rp * np.conj(rz)
    return {
        "frequency_hz": f.copy(),
        "coherence_squared": np.array(spectra["coherence_squared"], dtype=float),
        "expected_pa": np.array(spectra["expected_pa"], dtype=float),
        "segment_count": spectra["segment_count"],
        "cross_pa": cross,
        "phase_deg": np.angle(cross, deg=True),
        "psd_ratio_pa": np.asarray(spectra["psd_ratio_pa"]) * abs(rp) / abs(rz),
    }


def log_gain(observed, expected, coherence, segments, coherence_cap=0.999):
    """SOC-004 sensitivity fit to ln(T/R), with approximate H1 variance weights.

    The cap prevents singular weights at unity coherence. These weights are
    an experiment, not a calibrated uncertainty model for the PSD-ratio estimator.
    Hann-bin covariance and unknown channel-noise partition remain unmodeled.
    """
    r, t, c = [np.asarray(x, dtype=float) for x in (observed, expected, coherence)]
    if r.ndim != 1 or r.shape != t.shape or r.shape != c.shape or len(r) < 3:
        raise ValueError("At least three corresponding frequency bins required")
    if (
        not all(np.isfinite(x).all() for x in (r, t, c))
        or np.any(r <= 0)
        or np.any(t <= 0)
        or np.any(c <= 0)
        or np.any(c > 1)
        or segments < 1
        or not 0 < coherence_cap < 1
    ):
        raise ValueError("Invalid ratios, coherence, segments or weight cap")
    c = np.minimum(c, coherence_cap)
    weights = 2 * segments * c / (1 - c)
    target = np.log(t) - np.log(r)
    value = float(np.exp(np.average(target, weights=weights)))
    return value


def evaluate(spectra, settings, overlap, factor, method):
    """Same coherence/phase/overlap rules as production, optional fit objective."""
    c = np.asarray(spectra["coherence_squared"])
    phase = np.angle(spectra["cross_pa"] * factor, deg=True)
    coherent = c >= settings["coherence_squared_threshold"]
    valid = coherent & (abs(phase) <= settings["phase_limit_deg"])
    out = {
        "phase_valid_bins": int(valid.sum()),
        "median_phase_deg": circular_median_degrees(phase),
        "phase_summary_method": PHASE_SUMMARY_METHOD,
        "pressure_multiplier": None,
        "status": "rejected_qc",
    }
    threshold = settings["coherence_squared_threshold"]
    if np.median(c) < threshold or np.mean(c) < threshold or coherent.sum() < 3:
        return out
    if valid.sum() < 3:
        out["status"] = "polarity_review"
        return out
    if overlap:
        out["status"] = "arrival_overlap_review"
        return out
    r, t = spectra["psd_ratio_pa"][valid], spectra["expected_pa"][valid]
    try:
        if method == "absolute":
            g = fit_multiplier(r, t)[0]["pressure_multiplier"]
        elif method == "weighted_log":
            g = log_gain(r, t, c[valid], int(spectra["segment_count"]))
            low, high = settings["gain_grid_bounds"]
            if not low < g < high:
                raise ValueError("Weighted log gain outside unchanged trial bounds")
        else:
            raise ValueError("Unknown experimental fit method")
        out.update(
            status="candidate",
            pressure_multiplier=g,
            relative_rms=float(np.linalg.norm(g * r - t) / np.linalg.norm(t)),
        )
        if valid.sum() >= 4:
            # Remove the highest accepted bin; do not silently change the band.
            if method == "absolute":
                top_removed = fit_multiplier(r[:-1], t[:-1])[0]["pressure_multiplier"]
            else:
                top_removed = log_gain(
                    r[:-1], t[:-1], c[valid][:-1], int(spectra["segment_count"])
                )
            out["top_bin_removed_gain"] = top_removed
    except ValueError as exc:
        if out["pressure_multiplier"] is None:
            out.update(status="unresolved", fit_reason=str(exc))
        else:
            out["top_bin_removed_reason"] = str(exc)
    return out


def highpass_diagnostic(spectra, factor, threshold):
    """Joint log-amplitude/phase diagnostic; never creates an accepted gain.

    Model observed/theory = G*(i*f/fc)/(1+i*f/fc). G is the transfer's
    high-frequency limit, not a scalar pressure calibration for the fit band.
    """
    good = spectra["coherence_squared"] >= threshold
    if good.sum() < 4:
        return None
    f = spectra["frequency_hz"][good]
    ratio = spectra["psd_ratio_pa"][good] / spectra["expected_pa"][good]
    phase = np.angle(spectra["cross_pa"][good] * factor[good])

    def residual(parameters):
        gain, corner = np.exp(parameters)
        model = gain * (1j * f / corner) / (1 + 1j * f / corner)
        return np.r_[
            np.log(abs(model)) - np.log(ratio),
            np.angle(model * np.exp(-1j * phase)),
        ]

    fit = least_squares(
        residual,
        np.log([np.median(ratio), 0.03]),
        bounds=(np.log([1e-4, 1e-6]), np.log([100, 1])),
    )
    gain, corner = np.exp(fit.x)
    n = len(f)
    return {
        "diagnostic_only": True,
        "converged": bool(fit.success),
        "transfer_gain_G": float(gain),
        "corner_hz": float(corner),
        "time_constant_s": float(1 / (2 * np.pi * corner)),
        "log_amplitude_rms": float(np.sqrt(np.mean(fit.fun[:n] ** 2))),
        "phase_rms_deg": float(np.degrees(np.sqrt(np.mean(fit.fun[n:] ** 2)))),
        "fit_band_hz": [float(f.min()), float(f.max())],
        "corner_at_lower_bound": bool(corner < 1.01e-6),
    }
