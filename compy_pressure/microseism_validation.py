"""Frozen frequency-consistency rule for conditional microseism estimates.

Version 1 is transcribed from Claude's 15 October handoff (b603948). YL was
the design set. This is a reproducibility screen, not absolute calibration or
a confidence interval; correlated tilt and an incorrect vertical/model remain
possible even when every condition passes.
"""

import math
import statistics

RULE_ID = "claude-frequency-v1-20261015"
EXACT_PHASE_LABELS = frozenset(
    {"nominal", "reflected_pole", "reversed_nominal", "reversed_reflected_pole"}
)
RULE = {
    "id": RULE_ID,
    "source_commit": "b603948",
    "design_network": "YL",
    "maximum_sub_band_edge_hz": 0.10,
    "maximum_conditioning_ratio": 1.5,
    "minimum_eligible_sub_bands": 2,
    "maximum_sub_band_spread_pct": 10.0,
    "minimum_passing_weeks": 3,
    "maximum_week_spread_pct": 10.0,
    "stable_week_spread_pct": 5.0,
    "exact_phase_labels": sorted(EXACT_PHASE_LABELS),
    "coherence_policy": "reported descriptively; no added coherence threshold",
    "phase_policy": "legacy whole-band phase screen; not a gain-band tilt test",
}


def _positive(value):
    return isinstance(value, (int, float)) and math.isfinite(value) and value > 0


def frequency_week(check, phase_label):
    """Apply the exact frozen per-week rule, without any reference gain input."""
    eligible = []
    for band in check.get("long_wave_sub_bands", []):
        edges = band.get("frequency_hz", [])
        if (
            len(edges) != 2
            or not all(_positive(v) for v in edges)
            or edges[0] >= edges[1]
            or edges[1] > RULE["maximum_sub_band_edge_hz"]
        ):
            continue
        g1, g2, bound = (
            band.get(key)
            for key in ("mode0_water_h1", "mode0_water_h2", "upper_bound_h2")
        )
        if not all(_positive(v) for v in (g1, g2, bound)) or g1 > g2:
            continue
        conditioning = bound / g1
        if conditioning > RULE["maximum_conditioning_ratio"]:
            continue
        eligible.append(
            {
                "frequency_hz": list(edges),
                "g_h1": g1,
                "g_h2": g2,
                "midpoint_g": (g1 + g2) / 2,
                "conditioning_ratio": conditioning,
                "median_coherence_squared": band.get("median_coherence_squared"),
            }
        )
    # A duplicate/overlapping band must never masquerade as corroborating data.
    eligible.sort(key=lambda b: b["frequency_hz"])
    for left, right in zip(eligible, eligible[1:]):
        if left["frequency_hz"][1] > right["frequency_hz"][0]:
            raise ValueError("Eligible microseism sub-bands overlap or repeat")
    mids = [band["midpoint_g"] for band in eligible]
    spread = 100 * (max(mids) / min(mids) - 1) if mids else None
    reasons = []
    if phase_label not in EXACT_PHASE_LABELS:
        reasons.append("phase_not_explained_exactly")
    if len(eligible) < RULE["minimum_eligible_sub_bands"]:
        reasons.append("fewer_than_two_eligible_sub_bands")
    if spread is not None and spread > RULE["maximum_sub_band_spread_pct"] + 1e-12:
        reasons.append("sub_bands_disagree")
    passed = not reasons
    return {
        "rule_id": RULE_ID,
        "phase_label": phase_label,
        "passes": passed,
        "reasons": reasons,
        "eligible_sub_bands": eligible,
        "sub_band_spread_pct": spread,
        "conditional_week_g": statistics.median(mids) if passed else None,
        "g_h1": min(b["g_h1"] for b in eligible) if passed else None,
        "g_h2": max(b["g_h2"] for b in eligible) if passed else None,
    }


def frequency_station(weeks):
    """Summarise passing weeks; v1 needs >=3 weeks and <=10% spread.

    The caller verifies run provenance and non-overlapping time intervals.
    Sign consistency is reported separately: passing this amplitude rule does
    not certify polarity, any channel response, or an absolute pressure gain.
    """
    passed = [w for w in weeks if w.get("frequency_check", {}).get("passes")]
    values = [w["frequency_check"]["conditional_week_g"] for w in passed]
    signs = {w["sign"] for w in passed if w.get("sign") is not None}
    result = {
        "rule_id": RULE_ID,
        "passing_weeks": len(passed),
        "conditional_reportable": False,
        "median_g": statistics.median(values) if values else None,
        "week_spread_pct": None,
        "sign_stable": len(signs) == 1 if signs else None,
        "status": "insufficient_passing_weeks",
    }
    if not values:
        result["status"] = "no_passing_weeks"
        return result
    spread = 100 * (max(values) - min(values)) / statistics.median(values)
    result.update(
        min_week_g=min(values),
        max_week_g=max(values),
        week_spread_pct=spread,
        pooled_h1=min(w["frequency_check"]["g_h1"] for w in passed),
        pooled_h2=max(w["frequency_check"]["g_h2"] for w in passed),
    )
    if len(passed) < RULE["minimum_passing_weeks"]:
        return result
    if spread > RULE["maximum_week_spread_pct"] + 1e-12:
        result["status"] = "variable"
        return result
    result["conditional_reportable"] = True
    result["status"] = (
        "stable" if spread <= RULE["stable_week_spread_pct"] + 1e-12 else "moderate"
    )
    return result
