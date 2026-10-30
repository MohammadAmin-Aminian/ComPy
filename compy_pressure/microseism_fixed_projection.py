"""Out-of-training sensitivity to an already saved quiet-data projection.

No fitting, pressure reference or new physical rotation is performed here.
Weights are interpreted strictly in the recorded H1,H2,upward-Z basis.
"""

import numpy as np
from obspy import UTCDateTime

RULE_ID = "saved-quiet-projection-v1"
QUIET_SECONDS = 7200  # historical CleanRotator training documented in tilt_methods.md


def prior_projections(events, test_start):
    """Select solely by training time; never by event QC or target amplitudes.

    Primary = latest completed quiet interval before the earliest test sample.
    All earlier eligible projections are sensitivity cases, not independent
    calibration estimates. A future/overlapping fit cannot enter this test.
    """
    test_start = UTCDateTime(test_start)
    result = []
    for e in events:
        if not e.get("quiet_start_utc") or e.get("rotation_coefficients") is None:
            continue
        if e.get("tilt", {}).get("method") != "tiskitpy.CleanRotator":
            continue
        start = UTCDateTime(e["quiet_start_utc"])
        end = start + QUIET_SECONDS
        if end > test_start:
            continue
        if e["tilt"].get("response_output") != "VEL":
            raise ValueError(
                "Historical weights must come from physical velocity traces"
            )
        weights = np.asarray(e["rotation_coefficients"], dtype=float)
        validate_weights(weights)
        result.append(
            {
                "event_id": e["id"],
                "event_status_not_used_for_selection": e.get("status"),
                "training_start_utc": str(start),
                "training_end_utc": str(end),
                "weights_h1_h2_z": weights.tolist(),
                "tilt_source": e["tilt"],
            }
        )
    result.sort(key=lambda r: r["training_start_utc"])
    if not result:
        raise ValueError("No nonoverlapping prior quiet-data projection")
    if len({r["training_start_utc"] for r in result}) != len(result):
        raise ValueError("Ambiguous duplicate quiet-training intervals")
    for r in result:
        r["role"] = "primary_latest_prior" if r is result[-1] else "earlier_sensitivity"
    return result


def validate_weights(weights):
    w = np.asarray(weights, dtype=float)
    if (
        w.shape != (3,)
        or not np.isfinite(w).all()
        or not np.isclose(np.linalg.norm(w), 1, rtol=1e-8)
    ):
        raise ValueError("Expected three finite unit-norm H1,H2,Z weights")
    if w[2] < np.cos(np.radians(10)):
        raise ValueError("Projection must preserve upward near-vertical orientation")
    return w


def projection_products(coefficients, weights):
    """Apply fixed relative-basis weights to temporary Fourier copies only.

    Velocity-trained constant real weights commute with conversion to displacement
    on the same grid. Input order P,H1,H2,Z; transfer = mean(Z conj(P))/mean(|P|²).
    """
    x = np.asarray(coefficients, dtype=complex)
    if x.ndim != 3 or x.shape[0] < 6 or x.shape[1] != 4 or not np.isfinite(x).all():
        raise ValueError("Expected >=6 finite windows in P,H1,H2,Z order")
    w = validate_weights(weights)
    p, z = x[:, 0], x[:, 3]
    projected = np.einsum("i,wif->wf", w, x[:, 1:])
    pp, zz, rr = (np.mean(abs(v) ** 2, axis=0) for v in (p, z, projected))
    if np.any((pp <= 0) | (zz <= 0) | (rr <= 0)):
        raise ValueError("Positive spectral powers are required")
    zp, rp = (np.mean(v * p.conj(), axis=0) for v in (z, projected))
    valid = abs(zp) > 0
    ratio = np.full(zp.shape, np.nan + 1j * np.nan)
    ratio[valid] = rp[valid] / zp[valid]
    return {
        "nominal_transfer": zp / pp,
        "projected_transfer": rp / pp,
        "nominal_coherence_squared": np.clip(abs(zp) ** 2 / (pp * zz), 0, 1),
        "projected_coherence_squared": np.clip(abs(rp) ** 2 / (pp * rr), 0, 1),
        "vertical_power_ratio": rr / zz,
        "complex_transfer_ratio": ratio,
    }


def ideal_rayleigh_distortion(weights, ellipticities=(0.5, 1, 2, 5)):
    """Illustrative quadrature-polarized signal only; not an empirical bound.

    Sweeps every relative azimuth with a hypothetical orthogonal H pair. This
    neither verifies the native channel geometry nor covers arbitrary mixtures.
    """
    w = validate_weights(weights)
    theta = np.radians(np.arange(360))
    result = []
    for ratio in ellipticities:
        if not np.isfinite(ratio) or ratio <= 0:
            raise ValueError("Illustrative H/Z ellipticity must be positive")
        value = w[2] + 1j * ratio * (w[0] * np.cos(theta) + w[1] * np.sin(theta))
        result.append(
            {
                "horizontal_over_vertical_amplitude": float(ratio),
                "maximum_amplitude_change_pct": float(
                    np.max(abs(abs(value) - 1)) * 100
                ),
                "maximum_phase_change_deg": float(
                    np.max(abs(np.angle(value, deg=True)))
                ),
            }
        )
    return result
