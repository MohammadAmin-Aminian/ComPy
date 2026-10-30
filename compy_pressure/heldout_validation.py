"""Frozen-gain checks on unused events, with no calibration refitting."""

import numpy as np
from obspy import UTCDateTime

from compy_pressure.phase import PHASE_SUMMARY_METHOD, circular_median_degrees


def select_unused_events(catalogue, reviewed, minimum_magnitude=6.0):
    """Exclude all previously recorded IDs and nearby origins, not just fit events.

    Retain the existing Mw/depth/24 h rules. Select all remaining mainshocks,
    largest first; no success-based stopping or maximum event count.
    """
    if not np.isfinite(minimum_magnitude) or minimum_magnitude < 6:
        raise ValueError("Minimum magnitude must be at least 6")
    ids = {e["id"] for e in reviewed}
    origins = [
        UTCDateTime(e["origin_utc"]) for e in {e["id"]: e for e in reviewed}.values()
    ]
    eligible, excluded, selected = [], [], []
    for event in catalogue:
        origin = UTCDateTime(event["origin_utc"])
        if event["id"] in ids:
            reason = "previously_recorded_in_YV_trial"
        elif (
            event["magnitude"] < minimum_magnitude
            or event["depth_km"] >= 150
            or not str(event["magnitude_type"]).lower().startswith("mw")
        ):
            reason = "magnitude_depth_or_type"
        elif any(abs(origin - previous) < 86400 for previous in origins):
            reason = "within_24h_of_reviewed_event"
        else:
            eligible.append(event)
            continue
        excluded.append({"id": event["id"], "reason": reason})
    for event in sorted(
        eligible, key=lambda e: (-e["magnitude"], e["origin_utc"], e["id"])
    ):
        origin = UTCDateTime(event["origin_utc"])
        if any(abs(origin - UTCDateTime(e["origin_utc"])) < 86400 for e in selected):
            excluded.append({"id": event["id"], "reason": "heldout_mainshock_cluster"})
        else:
            selected.append(event)
    return selected, excluded


def frozen_gain_check(spectra, settings, gain, factor):
    """Measure prediction errors at fixed gain; never reject by amplitude error."""
    if not np.isfinite(gain) or gain <= 0:
        raise ValueError("Frozen gain must be positive and finite")
    c = np.asarray(spectra["coherence_squared"])
    r, t = spectra["psd_ratio_pa"], spectra["expected_pa"]
    phase = np.angle(spectra["cross_pa"] * factor, deg=True)
    if (
        not all(np.isfinite(a).all() for a in (c, r, t, phase))
        or np.any(r <= 0)
        or np.any(t <= 0)
    ):
        raise ValueError("Invalid held-out spectra")
    threshold = settings["coherence_squared_threshold"]
    coherent = c >= threshold
    valid = coherent & (abs(phase) <= settings["phase_limit_deg"])
    relative = gain * r / t - 1
    result = {
        "frozen_gain": float(gain),
        "gain_refitted": False,
        "phase_valid_bins": int(valid.sum()),
        "median_coherence": float(np.median(c)),
        "median_phase_deg": circular_median_degrees(phase),
        "phase_summary_method": PHASE_SUMMARY_METHOD,
        "all_band_median_error_percent": float(100 * np.median(relative)),
        "status": "coherence_review",
    }
    if np.median(c) < threshold or np.mean(c) < threshold or coherent.sum() < 3:
        return result
    if valid.sum() < 3:
        result["status"] = "phase_review"
        return result
    residual = relative[valid]
    result.update(
        status="qc_usable",
        median_error_percent=float(100 * np.median(residual)),
        rms_error_percent=float(100 * np.sqrt(np.mean(residual**2))),
        frequency_error_percent=(100 * relative).tolist(),
        phase_valid_mask=valid.tolist(),
    )
    if valid.sum() >= 4:
        result["top_bin_removed_median_error_percent"] = float(
            100 * np.median(residual[:-1])
        )
    return result


def summarize_errors(events, *, require_isolated=True):
    """Conditional event spread, not absolute calibration uncertainty or scoring."""
    errors = [
        e["reflected"]["median_error_percent"]
        for e in events
        if e.get("reflected", {}).get("status") == "qc_usable"
        and not e.get("other_packet_ids")
        and (not require_isolated or not e.get("body_arrivals_in_window"))
    ]
    result = {
        "usable_events": len(errors),
        "status": "heldout_review"
        if len(errors) >= 3
        else "insufficient_heldout_events",
        "uncertainty_kind": "conditional event spread; systematics omitted",
    }
    if errors:
        result.update(
            median_error_percent=float(np.median(errors)),
            event_error_percentiles_16_84=np.percentile(errors, [16, 84]).tolist(),
            event_error_range_percent=[float(min(errors)), float(max(errors))],
        )
    if len(errors) >= 5:
        values = np.asarray(errors)
        rng = np.random.default_rng(0)
        medians = np.median(
            rng.choice(values, size=(5000, len(values)), replace=True), axis=1
        )
        result["event_resampling_interval_of_median_error_percent"] = np.percentile(
            medians, [2.5, 97.5]
        ).tolist()
        result["interval_note"] = (
            "IID-event descriptive bootstrap, not calibration confidence"
        )
    return result
