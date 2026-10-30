"""Known complex count gains and sign, independent of the water reference."""

import numpy as np
import pytest
from obspy import Trace, UTCDateTime

from compy_pressure.rayleigh_validation import count_domain_products
from scripts.validate_yv_rayleigh import load_counts


def test_independent_welch_recovers_known_pressure_acceleration_and_rotation():
    rate, f0, pressure_per_acceleration = 4, 0.045, 2500000
    t = np.arange(4 * 1800) / rate
    weights = np.array([0.1, 0.2, np.sqrt(0.95)])
    amplitude = np.array([4, 8, 12, 16])
    phase = np.array([0.7, -0.3, 0.2, -0.5])
    # Velocity is cos; upward acceleration is -omega*sin. Pressure is reversed.
    physical_amp = np.r_[pressure_per_acceleration * 2 * np.pi * f0, weights]
    signal_phase = np.r_[-np.pi / 2, np.zeros(3)]
    counts = (
        amplitude[:, None]
        * physical_amp[:, None]
        * np.cos(2 * np.pi * f0 * t + (phase + signal_phase)[:, None])
    )
    before = counts.copy()
    complex_response = amplitude * np.exp(1j * phase)

    def response(f):
        return np.broadcast_to(complex_response[:, None], (4, len(f)))

    result = count_domain_products(counts, rate, response, weights, (0.03, 0.07))
    at = np.argmin(abs(result["frequency_hz"] - f0))
    assert result["psd_ratio_pa"][at] == pytest.approx(pressure_per_acceleration)
    assert result["phase_deg"][at] == pytest.approx(0, abs=1e-10)
    assert result["coherence_squared"][at] == pytest.approx(1)
    opposite = count_domain_products(
        counts, rate, response, weights, (0.03, 0.07), pressure_sign=1
    )
    assert abs(opposite["phase_deg"][at]) == pytest.approx(180)
    np.testing.assert_array_equal(counts, before)


def test_independent_welch_rejects_invalid_counts_and_response():
    x = np.ones((4, 7200))
    weights = [0, 0, 1]
    with pytest.raises(ValueError, match="band"):
        count_domain_products(x, 4, lambda f: np.ones((4, len(f))), weights, (0, 0.07))
    with pytest.raises(ValueError, match="response"):
        count_domain_products(
            x, 4, lambda f: np.zeros((4, len(f))), weights, (0.03, 0.07)
        )


@pytest.mark.parametrize("fractional_shift", [0, 0.0016])
def test_independent_reader_checks_duplicate_overlap_values_and_timing(
    tmp_path, fractional_shift
):
    start = UTCDateTime("2020-01-01")
    paths = []
    for i, offset in enumerate([0, 30 - fractional_shift]):
        tr = Trace(np.arange(i * 1500, i * 1500 + 2000, dtype=np.int32))
        tr.stats.update(
            dict(
                network="XX",
                station="TEST",
                channel="BHZ",
                sampling_rate=50,
                starttime=start + offset,
            )
        )
        path = tmp_path / f"XX.TEST..BHZ.D.2020.{i + 1:03}"
        tr.write(str(path), format="MSEED")
        paths.append(str(path))
    if fractional_shift:
        with pytest.raises(ValueError, match="Fractional sample-grid"):
            load_counts(paths, "XX.TEST..BHZ", start, 60, 50)
    else:
        np.testing.assert_array_equal(
            load_counts(paths, "XX.TEST..BHZ", start, 60, 50), np.arange(3000)
        )
