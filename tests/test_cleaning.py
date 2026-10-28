"""Synthetic tests for unattended cleaning steps."""

import numpy as np
import pytest
from obspy import Stream, Trace, UTCDateTime

from compy_cleaning import (
    periodic_template,
    remove_periodic_transient,
    rotate_smoothly,
    rotation_coefficients,
    sta_lta_spans,
)

FS = 1.25


def trace(data, channel="BHZ", start=UTCDateTime(2013, 1, 1)):
    return Trace(
        np.asarray(data, float),
        header={"channel": channel, "sampling_rate": FS, "starttime": start,
                "network": "XX", "station": "T"},
    )  # fmt: skip


def glitch_series(period, days, amplitude, seed):
    rng = np.random.default_rng(seed)
    n = int(days * 86400 * FS)
    t = np.arange(n) / FS
    phase = np.mod(t, period)
    pulse = (
        amplitude
        * np.exp(-(((phase - 600) / 120) ** 2))
        * np.sin(2 * np.pi * phase / 400)
    )
    return t, pulse, rng.normal(0, 0.2 * amplitude, n)


def test_periodic_template_recovers_period_and_removes_transient():
    t, pulse, noise = glitch_series(3599.37, 20, 1.0, 1)
    tr = trace(pulse + noise)
    template = periodic_template(tr, period_step=0.02)
    assert template.period == pytest.approx(3599.37, abs=0.02)
    assert template.strength > 5
    cleaned = remove_periodic_transient(tr, template)
    residual = cleaned.data - noise
    # Floor: median-stack noise, ~1.25 * 0.2 / sqrt(480 cycles) = 0.011.
    assert np.std(residual) < 0.15 * np.std(pulse)


def test_periodic_template_reports_absence_of_transient():
    _, _, noise = glitch_series(3600, 10, 1.0, 2)
    template = periodic_template(trace(noise))
    assert template.strength < 3


@pytest.mark.parametrize("angle,azimuth", [(0.0, 0.0), (2.0, 45.0), (5.0, 250.0)])
def test_rotation_coefficients_match_tiskitpy(angle, azimuth):
    from tiskitpy.util_classes.seis_rotate import SeisRotate

    rng = np.random.default_rng(3)
    st = Stream([trace(rng.normal(size=100), c) for c in ("BHZ", "BH1", "BH2")])
    rotator = SeisRotate(st.copy())
    rotator.zrotate(angle, azimuth, horiz_also=False)
    cz, c1, c2 = rotation_coefficients(angle, azimuth)
    expected = cz * st[0].data + c1 * st[1].data + c2 * st[2].data
    np.testing.assert_allclose(rotator.Z.data, expected, rtol=1e-10, atol=1e-12)


def test_rotate_smoothly_removes_constant_tilt_without_steps():
    rng = np.random.default_rng(4)
    n = int(2 * 86400 * FS)
    h1 = np.cumsum(rng.normal(size=n))
    h2 = np.cumsum(rng.normal(size=n))
    z0 = rng.normal(size=n)
    angle, azimuth = 1.5, 30.0
    cz, c1, c2 = rotation_coefficients(angle, azimuth)
    # Data whose vertical is the inverse-rotated "true" vertical.
    leak1, leak2 = -c1 / cz, -c2 / cz
    st = Stream([trace(z0 + leak1 * h1 + leak2 * h2, "BHZ"), trace(h1, "BH1"),
                 trace(h2, "BH2")])  # fmt: skip
    hours = [st[0].stats.starttime + 1800 + 3600 * i for i in range(48)]
    out = rotate_smoothly(st, hours, [angle] * 48, [azimuth] * 48)
    np.testing.assert_allclose(out[0].data / cz, z0, atol=1e-6 * np.std(h1))


def test_sta_lta_finds_injected_event():
    rng = np.random.default_rng(5)
    n = int(86400 * FS)
    x = rng.normal(size=n)
    onset = int(43200 * FS)
    x[onset : onset + 600] += 40 * np.sin(np.arange(600) * 2 * np.pi * 0.1 / FS)
    spans = sta_lta_spans(trace(x))
    assert len(spans) == 1
    start, end = spans[0]
    t_event = trace(x).stats.starttime + 43200
    assert start <= t_event <= end
