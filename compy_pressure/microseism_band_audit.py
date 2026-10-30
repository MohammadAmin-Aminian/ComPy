"""Descriptive phase audit of the exact bands selected by the frozen v1 rule.

This does not change v1 selection, estimate a gain, or establish that tilt is
absent. In particular a real, pressure-coherent contribution to the vertical
can bias amplitude without failing a phase test.
"""

import numpy as np

from compy_pressure.microseism_validation import EXACT_PHASE_LABELS

COHERENT = 0.5
MINIMUM_BINS = 10
FITS_DEG = 10.0
UNEXPLAINED_DEG = 20.0
AUDIT_ID = "gain-band-phase-audit-v1"


def wrap(degrees):
    return (np.asarray(degrees) + 180.0) % 360.0 - 180.0


def hypotheses(f, pressure_poles=(), vertical_poles=()):
    """Same four hypotheses as the legacy screen, without fitting parameters.

    Positive-real poles are in rad/s; pressure and vertical enter with opposite
    signs in H = D/P. Reflection and polarity affect phase only, not |H|.
    """
    f = np.asarray(f, dtype=float)
    if np.any(~np.isfinite(f) | (f <= 0)):
        raise ValueError("Frequencies must be positive and finite")
    pole = np.zeros_like(f)
    for poles, sign in ((pressure_poles, -1), (vertical_poles, 1)):
        for p in poles:
            if not np.isfinite(p) or p <= 0:
                raise ValueError("Positive-real poles must be positive rad/s")
            pole += sign * np.degrees(2 * np.arctan(p / (2 * np.pi * f)))
    result = {"nominal": np.zeros_like(f), "reversed_nominal": np.full_like(f, 180)}
    if pressure_poles or vertical_poles:
        result |= {"reflected_pole": pole, "reversed_reflected_pole": 180 + pole}
    return result


def phase_band(
    f, phase, coherence, edges, prescribed, pressure_poles=(), vertical_poles=()
):
    """Report RMS against the *prescribed* whole-band label, never reselect it.

    Band edges are [low, high), matching the saved gain estimator. The 0.5,
    10-bin and 10/20-degree descriptive cutoffs come from the legacy screen.
    A per-band best alternative is reported only as a diagnostic.
    """
    f, phase, coherence = [np.asarray(a, dtype=float) for a in (f, phase, coherence)]
    if f.ndim != 1 or not (f.shape == phase.shape == coherence.shape):
        raise ValueError("Frequency, phase and coherence must be equal 1D arrays")
    if np.any(np.diff(f) <= 0):
        raise ValueError("Frequencies must increase strictly")
    if len(edges) != 2 or not 0 < edges[0] < edges[1] or not np.isfinite(edges).all():
        raise ValueError("Expected finite positive half-open band edges")
    if prescribed not in EXACT_PHASE_LABELS:
        raise ValueError("Expected an exact frozen-rule phase label")
    models = hypotheses(f, pressure_poles, vertical_poles)
    if prescribed not in models:
        raise ValueError("Reflected label requires recorded positive-real poles")
    inside = (f >= edges[0]) & (f < edges[1])
    if np.any(~np.isfinite(phase[inside])) or np.any(
        ~np.isfinite(coherence[inside])
        | (coherence[inside] < 0)
        | (coherence[inside] > 1 + 1e-9)
    ):
        raise ValueError("Nonfinite phase or invalid coherence inside gain band")
    coherent = inside & (coherence >= COHERENT)

    def rms(mask, model):
        return (
            float(np.sqrt(np.mean(wrap(phase[mask] - model[mask]) ** 2)))
            if mask.any()
            else None
        )

    residuals = {name: rms(coherent, model) for name, model in models.items()}
    best = min(residuals, key=residuals.get) if coherent.any() else None
    prescribed_rms = residuals[prescribed]
    status = "insufficient_coherent_bins"
    if coherent.sum() >= MINIMUM_BINS:
        status = (
            "consistent_with_prescribed"
            if prescribed_rms <= FITS_DEG
            else "approximate_prescribed"
            if prescribed_rms <= UNEXPLAINED_DEG
            else "inconsistent_with_prescribed"
        )
    result = {
        "frequency_hz": list(edges),
        "edge_convention": "low <= f < high",
        "bins": int(inside.sum()),
        "coherent_bins": int(coherent.sum()),
        "coherent_fraction": float(coherent.sum() / inside.sum())
        if inside.any()
        else None,
        "prescribed_hypothesis": prescribed,
        "prescribed_rms_deg": prescribed_rms,
        "all_bin_prescribed_rms_deg": rms(inside, models[prescribed]),
        "best_alternative_diagnostic": best,
        "rms_by_hypothesis_deg": residuals,
        "status": status,
    }
    return result


def model_band_context(bins, mask, reference=None):
    """Compare cached model variants on identical bins, without a new inversion.

    Cached per-bin g values have six significant figures. This sensitivity is
    only water-loaded vs dry historical mode 0, NOT an Earth-model uncertainty
    budget. H1/H2 are noise-bias diagnostics, not confidence bounds.
    """
    result = {}
    for name, key in (
        ("water_h1", "g_mode0_water_h1"),
        ("water_h2", "g_mode0_water_h2"),
        ("dry_h1", "g_mode0_no_water_h1"),
        ("dry_h2", "g_mode0_no_water_h2"),
        ("upper_h2", "upper_bound_h2"),
    ):
        values = np.asarray(bins[key])[mask]
        if not len(values) or not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError(
                "Gain-band cached model values must be positive and finite"
            )
        result[name] = float(np.median(values))
    result["dry_vs_water_h1_pct"] = 100 * (result["dry_h1"] / result["water_h1"] - 1)
    if reference is not None:
        if not np.isfinite(reference) or reference <= 0:
            raise ValueError("Fixed comparison reference must be positive and finite")
        result |= {
            "fixed_rayleigh_g": reference,
            "reference_inside_model_h1_h2": result["water_h1"]
            <= reference
            <= result["water_h2"],
            "required_height_over_model_h1": reference / result["water_h1"],
            "required_height_over_model_h2": reference / result["water_h2"],
        }
    return result
