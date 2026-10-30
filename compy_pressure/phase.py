"""Circular phase summaries; descriptive only, never a scientific QC threshold."""

import hashlib
from pathlib import Path

import numpy as np

PHASE_SUMMARY_METHOD = "circular_l1_median_v1"
PHASE_CODE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def circular_median_degrees(values):
    """Minimize total shortest angular distance, returning [-180, 180).

    A circular-mean-centred median is preferred when it is a global L1 minimizer;
    otherwise choose the minimizing observation nearest that centre. Prefix sums
    make the candidate search O(n log n), without a quadratic distance matrix.
    Empty/nonfinite-only or numerically directionless data return None, not a
    fabricated zero-degree phase. No coherence/polarity acceptance rule is implied.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 1:
        raise ValueError("Phase summary requires a one-dimensional array")
    values = values[np.isfinite(values)]
    if not values.size:
        return None
    wrapped = (values + 180) % 360 - 180
    resultant = np.mean(np.exp(1j * np.radians(wrapped)))
    epsilon = 64 * np.finfo(float).eps
    if abs(resultant) <= epsilon:
        return None
    centre = np.degrees(np.angle(resultant))
    candidate = float(
        (centre + np.median((wrapped - centre + 180) % 360 - 180) + 180) % 360 - 180
    )
    angles = np.sort(wrapped % 360)
    extended = np.concatenate((angles - 360, angles, angles + 360))
    prefix = np.concatenate(([0.0], np.cumsum(extended)))
    left = np.searchsorted(extended, angles - 180, side="left")
    middle = np.searchsorted(extended, angles, side="left")
    right = left + len(angles)
    losses = (
        angles * (middle - left)
        - (prefix[middle] - prefix[left])
        + (prefix[right] - prefix[middle])
        - angles * (right - middle)
    )
    minimum = losses.min()
    tolerance = epsilon * len(angles) * 360
    candidate_loss = np.abs((wrapped - candidate + 180) % 360 - 180).sum()
    if abs(candidate_loss - minimum) <= tolerance:
        return candidate
    winners = angles[np.abs(losses - minimum) <= tolerance]
    winner = winners[np.argmin(np.abs((winners - centre + 180) % 360 - 180))]
    return float((winner + 180) % 360 - 180)
