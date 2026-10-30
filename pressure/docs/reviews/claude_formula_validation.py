"""Synthetic checks for docs/reviews/claude_formula_review.md (C1-C6).

Read-only: analytic/ODE solutions and synthetic signals only; no waveforms,
masters or saved results are read and nothing is written. Run from the root:

    PYTHONPATH=.. python docs/reviews/claude_formula_validation.py
"""

import warnings
from pathlib import Path

import numpy as np
from obspy import Trace, UTCDateTime
from obspy.core.inventory import Channel, Inventory, Network, Response, Station
from obspy.signal import PPSD
from scipy.integrate import solve_ivp

from compy_pressure.infragravity_test import dispersion_kh, godin_pressure_db
from compy_pressure.local_p_check import model_transfer
from compy_pressure.microseism_check import admittance_height
from compy_pressure.rayleigh_calibration import water_ratio

ROOT = Path(__file__).resolve().parents[2]
RHO, ALPHA, G = 1028.0, 1500.0, 9.81


def ode_height(f, c, depth, gravity):
    """Integrate the acoustic ODE floor->surface; Lagrangian P / (rho a_up)."""
    w = 2 * np.pi * f
    k = w / c

    def rhs(z, y):  # y = [Eulerian pressure, upward displacement]
        return [RHO * w**2 * y[1], y[0] * (k**2 / w**2 - 1 / ALPHA**2) / RHO]

    def surface_residual(p0):
        sol = solve_ivp(rhs, (-depth, 0), [p0, 1.0], rtol=1e-11, atol=1e-14)
        return sol.y[0, -1] - RHO * gravity * sol.y[1, -1]

    r0, r1 = surface_residual(0.0), surface_residual(1.0)  # linear in p0
    p0 = -r0 / (r1 - r0)
    return (p0 - RHO * gravity) / (-RHO * w**2)


def c1_water_column():
    worst = 0.0
    for depth in (900, 2918, 5000):
        for f in (0.02, 0.05, 0.08):
            for c in (900.0, 3900.0, 1e9):
                q = 1 / ALPHA**2 - 1 / c**2
                if q > 0 and 2 * np.pi * f * depth * np.sqrt(q) > 1.3:
                    continue
                for gravity in (0.0, G):
                    ode = ode_height(f, c, depth, gravity)
                    code = float(admittance_height(f, c, depth, ALPHA, gravity=gravity))
                    worst = max(worst, abs(code / ode - 1))
                    if gravity == 0 and ALPHA < c < 1e8:
                        rayleigh = water_ratio([f], [c], depth, RHO, ALPHA)[0] / RHO
                        worst = max(worst, abs(rayleigh / ode - 1))
    print(f"C1 water column vs ODE (gravity, evanescent): max rel. error {worst:.1e}")


def c2_godin():
    for rho in (1000.0, 1025.0, 1028.0):
        value = 2.1e-5 * (rho * G) ** 2  # Zabotin et al. (2015) B = 2.1e-5 m^2/Hz
        print(
            f"C2 Godin constant rho={rho:.0f}: B(rho g)^2 = {value:.0f} Pa^2/Hz, "
            f"{10 * np.log10(value / 2000):+.2f} dB vs code 2000"
        )
    x = dispersion_kh([0.005, 0.02], 4000.0)
    target = (2 * np.pi * np.array([0.005, 0.02])) ** 2 * 4000 / G
    print(
        f"C2 dispersion x tanh x residual {np.abs(x * np.tanh(x) - target).max():.0e}"
    )


def c3_depth_sensitivity():
    for depth in (2000.0, 3000.0, 4500.0):
        parts = []
        for f in (0.005, 0.01, 0.015, 0.02):
            base = godin_pressure_db([f], depth)[0]
            d5 = godin_pressure_db([f], depth * 1.05)[0] - base
            d10 = godin_pressure_db([f], depth * 1.10)[0] - base
            parts.append(f"{f * 1e3:.0f} mHz {d5:+.1f}/{d10:+.1f}")
        print(f"C3 Godin level, depth +5%/+10% at {depth:.0f} m: " + "; ".join(parts))


def c4_local_p_sign():
    rate, n, depth = 50.0, 2**14, 300.0
    t = np.arange(n) / rate
    v = np.exp(-(((t - 160) / 2.0) ** 2)) * np.sin(2 * np.pi * 0.5 * (t - 160))
    options = {"density": RHO, "sound_speed": ALPHA, "resonance_margin": 0.9}
    f = np.fft.rfftfreq(n, 1 / rate)
    model = model_transfer(f, "steady_state", depth, 1.0, options)
    low = f < 0.15
    p_model = np.fft.irfft(
        np.where(np.isfinite(model) & low, model, 0) * np.fft.rfft(v)
    )
    dvdt = np.gradient(v, 1 / rate)
    p_time = np.fft.irfft(np.where(low, np.fft.rfft(RHO * depth * dvdt), 0), n)
    print(
        "C4 local-P model vs rho*h*dv/dt: corr "
        f"{np.corrcoef(p_model, p_time)[0, 1]:.6f}, "
        f"scale {np.std(p_model) / np.std(p_time):.4f} (F=tan/theta ~1.01)"
    )


def c5_water_layer_dispersion():
    try:
        from disba import PhaseDispersion
    except ImportError:
        print("C5 skipped: disba not installed")
        return
    layers = np.loadtxt(
        ROOT / "../compy_pressure/reference_data/rr52_historical_rayleigh_layers.csv",
        delimiter=",",
        skiprows=1,
    )
    wet = np.vstack([[2.918, 1.5, 0.0, 1.028], layers])
    freqs = np.array([0.03, 0.05, 0.07])
    periods = np.sort(1 / freqs)
    dry = PhaseDispersion(*layers.T)(periods, mode=0, wave="rayleigh").velocity
    water = PhaseDispersion(*wet.T)(periods, mode=0, wave="rayleigh").velocity
    change = water_ratio(freqs, water[::-1] * 1e3, 2918.0, RHO, ALPHA) / water_ratio(
        freqs, dry[::-1] * 1e3, 2918.0, RHO, ALPHA
    )
    print(
        "C5 RR52 water layer in dispersion, T change at 0.03/0.05/0.07 Hz: "
        + ", ".join(f"{100 * (v - 1):+.2f}%" for v in change)
    )


def c6_ppsd_pressure_level():
    rate, sensitivity, sigma = 62.5, 2.0e4, 3.0
    response = Response.from_paz(
        zeros=[],
        poles=[],
        stage_gain=sensitivity,
        stage_gain_frequency=1.0,
        input_units="PA",
        output_units="COUNTS",
    )
    epoch = dict(start_date=UTCDateTime(2019, 1, 1), end_date=UTCDateTime(2021, 1, 1))
    channel = Channel(
        "BDH", "00", 0, 0, -3000, 0, sample_rate=rate, response=response, **epoch
    )
    inventory = Inventory(
        [Network("XX", [Station("TEST", 0, 0, -3000, channels=[channel])])]
    )
    data = np.random.default_rng(1).normal(0, sigma, int(rate * 86400))
    trace = Trace(
        data * sensitivity,
        header=dict(
            network="XX",
            station="TEST",
            location="00",
            channel="BDH",
            sampling_rate=rate,
            starttime=UTCDateTime(2020, 1, 1),
        ),
    )
    ppsd = PPSD(
        trace.stats,
        inventory,
        ppsd_length=10800,
        overlap=0.5,
        special_handling="hydrophone",
        db_bins=(-250, 100, 1.0),
        period_limits=(2 / rate, 1000),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ppsd.add(trace)
    values = np.asarray(ppsd.psd_values)
    expected = 10 * np.log10(2 * sigma**2 / rate)
    f = 1 / np.asarray(ppsd.period_bin_centers)
    parts = [
        f"{lo:g}-{hi:g} Hz {values[:, (f >= lo) & (f < hi)].mean() - expected:+.2f} dB"
        for lo, hi in ((0.001, 0.003), (0.01, 0.1), (0.1, 1), (1, 25))
    ]
    print(f"C6 PPSD white pressure noise minus 2*sigma^2/fs: {'; '.join(parts)}")


if __name__ == "__main__":
    c1_water_column()
    c2_godin()
    c3_depth_sensitivity()
    c4_local_p_sign()
    c5_water_layer_dispersion()
    c6_ppsd_pressure_level()
