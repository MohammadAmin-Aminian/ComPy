from types import SimpleNamespace

import numpy as np
import pytest
from obspy import Stream, Trace, UTCDateTime

import Pressure_calibration as calibration


def test_theoretical_ratio_passes_model_and_handles_zero_vertical_wavenumber(
    monkeypatch,
):
    model = np.array([[1, 3, 2, 2], [0, 4, 3, 3]])

    def dispersion(velocity_model, plot_condition, t):
        np.testing.assert_equal(velocity_model, model)
        return np.array([0.01]), np.array([1.5])

    monkeypatch.setattr(calibration, "phase_dispersion", dispersion)
    f, ratio = calibration.theoretical_p_a_ratio(
        alpha=1500, plot_condition=False, velocity_model=model
    )
    np.testing.assert_equal(ratio, [1])


def test_synthetic_event_calibration_includes_single_event_and_preserves_input(
    monkeypatch,
):
    depth, rho, fs = 4000.0, 1028.0, 2.0
    p = np.random.default_rng(7).normal(size=14400)
    st = Stream(
        [
            Trace(
                data.copy(),
                header={
                    "channel": channel,
                    "sampling_rate": fs,
                    "starttime": UTCDateTime(0),
                    "network": "XX",
                    "station": "TEST",
                },
            )
            for channel, data in [("BDH", p), ("BHZ", p / (depth * rho))]
        ]
    )
    before = st.copy()
    inventory = SimpleNamespace(get_coordinates=lambda *args: {"elevation": -depth})
    spans = SimpleNamespace(start_times=[UTCDateTime(0)])
    # The synthetic traces are already in physical units; this isolates event
    # handling/spectra/gain fitting from instrument-response implementation.
    monkeypatch.setattr(Stream, "remove_response", lambda self, **kwargs: self)

    def theory(*, t, **kwargs):
        return 1 / t, np.full(len(t), 0.6637)

    monkeypatch.setattr(calibration, "theoretical_p_a_ratio", theory)
    gain = calibration.calculate_spectral_ratio(
        st, inventory=inventory, event_spans=spans
    )
    assert np.isclose(gain, 0.6637, rtol=1e-12)
    assert st == before


def test_real_disba_dispersion_smoke():
    model = np.array([[1, 5, 3, 2.5], [0, 7, 4, 3.0]])
    f, v = calibration.phase_dispersion(
        model, plot_condition=False, t=np.array([10.0, 20.0, 30.0])
    )
    assert len(f) == len(v) == 3
    assert np.all(v > 0)


def test_rayleigh_window_uses_peak_magnitude():
    x = np.zeros(7200)
    x[3600] = -10
    st = Stream(
        [
            Trace(
                x,
                header={
                    "channel": "BHZ",
                    "sampling_rate": 2.0,
                    "starttime": UTCDateTime(0),
                },
            )
        ]
    )
    _, t1, t2 = calibration.rayleigh_arrival(st, timelag=-2, window=20)
    assert t1 == 1680
    assert t2 == 2880
