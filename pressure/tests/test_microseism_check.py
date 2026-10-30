import hashlib
import json
from importlib.resources import files
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from compy_pressure.microseism_check import (
    STANDARD_GRAVITY,
    admittance_height,
    first_resonance_hz,
    modal_phase_velocity,
    phase_velocity_from_height,
    run_check,
    vertical_incidence_limit,
)
from compy_pressure.rayleigh_calibration import water_ratio

DEPTH, SPEED, DENSITY = 2918.0, 1500.0, 1028.0
REFERENCE = Path(str(files("compy_pressure").joinpath("reference_data")))


def numerical_height(f, c, depth, speed, gravity):
    """Integrate p'' = nu^2 p downward from a gravity free surface."""
    omega = 2 * np.pi * f
    nu2 = omega**2 * (1 / c**2 - 1 / speed**2)
    rho = 1.0  # cancels in the height
    p0, dp0 = -rho * gravity, rho * omega**2  # u_z(0) = 1 (downward)
    solution = solve_ivp(
        lambda z, s: [s[1], nu2 * s[0]],
        (0, depth),
        [p0, dp0],
        rtol=1e-12,
        atol=1e-12,
    )
    p, dp = solution.y[:, -1]
    u_down = dp / (rho * omega**2)
    gauge = p + rho * gravity * u_down  # Lagrangian pressure on moving seafloor
    return gauge / (rho * omega**2 * u_down)


@pytest.mark.parametrize(
    "f,c",
    [(0.05, 3900.0), (0.15, 2243.0), (0.15, 4070.0), (0.3, 800.0), (0.07, np.inf)],
)
@pytest.mark.parametrize("gravity", [0.0, STANDARD_GRAVITY])
def test_height_matches_numerical_water_column(f, c, gravity):
    expected = numerical_height(f, c, DEPTH, SPEED, gravity)
    actual = admittance_height(f, c, DEPTH, SPEED, gravity)
    assert actual == pytest.approx(expected, rel=1e-8)


def test_height_matches_existing_first_branch_implementation():
    f = np.linspace(0.01, 0.07, 25)
    c = np.linspace(2500, 4500, 25)
    expected = water_ratio(f, c, DEPTH, DENSITY, SPEED) / DENSITY
    np.testing.assert_allclose(admittance_height(f, c, DEPTH, SPEED), expected)


def test_long_wave_limit_is_water_depth_for_any_slowness():
    for c in (300.0, 1500.0, 4000.0, np.inf):
        assert admittance_height(1e-5, c, DEPTH, SPEED) == pytest.approx(DEPTH)


def test_vertical_incidence_bounds_heights_only_below_resonance():
    resonance = first_resonance_hz(DEPTH, SPEED)
    assert resonance == pytest.approx(SPEED / (4 * DEPTH))
    speeds = np.r_[np.geomspace(100, 1e6, 400), np.inf]
    below = np.linspace(0.01, 0.99 * resonance, 12)
    heights = admittance_height(below[:, None], speeds[None, :], DEPTH, SPEED)
    limit = vertical_incidence_limit(below, DEPTH, SPEED)
    assert np.all(heights > 0)
    assert np.all(heights <= limit[:, None] * (1 + 1e-12))
    above = admittance_height(1.2 * resonance, speeds, DEPTH, SPEED)
    assert above.min() < 0 < above.max()
    assert np.isinf(vertical_incidence_limit([1.2 * resonance], DEPTH, SPEED))


def test_phase_velocity_inversion_round_trip_and_refusals():
    f = np.array([0.03, 0.06, 0.1, 0.15, 0.2, 0.25])
    c = np.array([3900.0, 3500.0, 1200.0, 2200.0, 1700.0, 600.0])
    height = admittance_height(f, c, DEPTH, SPEED)
    recovered = phase_velocity_from_height(f, height, DEPTH, SPEED)
    np.testing.assert_allclose(recovered, c, rtol=1e-8)
    limit = vertical_incidence_limit([0.06], DEPTH, SPEED)[0]
    refused = phase_velocity_from_height(
        [0.06, 0.06, 0.15], [1.001 * limit, -DEPTH, 0.0], DEPTH, SPEED
    )
    assert np.isnan(refused).all()
    assert np.isinf(phase_velocity_from_height(0.06, limit, DEPTH, SPEED))


def test_disba_water_layer_matches_independent_secular_equation():
    pytest.importorskip("disba")
    depth, vp, vs, rho_s = 2000.0, 6000.0, 3500.0, 2700.0

    def secular(c, f):
        k = 2 * np.pi * f / c
        a, b = np.sqrt(1 - c**2 / vp**2), np.sqrt(1 - c**2 / vs**2)
        rayleigh = (2 - c**2 / vs**2) ** 2 - 4 * a * b
        loading = DENSITY / rho_s * (c / vs) ** 4 * a * k
        return rayleigh + loading * admittance_height(f, c, depth, SPEED)

    layers = np.array([[0.0, vp / 1e3, vs / 1e3, rho_s / 1e3]])
    for f in (0.01, 0.05, 0.1, 0.2):
        grid = np.linspace(100, vs * (1 - 1e-6), 20000)
        values = np.array([secular(c, f) for c in grid])
        roots = [
            brentq(secular, grid[i], grid[i + 1], args=(f,), xtol=1e-10)
            for i in np.flatnonzero(np.sign(values[:-1]) != np.sign(values[1:]))
        ]
        roots = [r for r in roots if abs(secular(r, f)) < 1e-6]
        solved = modal_phase_velocity(
            layers, [f], modes=(0,), water=(depth, SPEED, DENSITY)
        )
        assert solved[0, 0] == pytest.approx(min(roots), rel=1e-5)


def test_dry_model_reproduces_bundled_reference_dispersion():
    pytest.importorskip("disba")
    from compy_pressure.microseism_check import load_layers

    table = np.loadtxt(
        REFERENCE / "rr52_historical_rayleigh_dispersion.csv", delimiter=",", skiprows=1
    )[1:-1:7]  # interior rows; end frequencies were rounded to the band
    layers = load_layers(REFERENCE / "rr52_historical_rayleigh_layers.csv")
    solved = modal_phase_velocity(layers, table[:, 0], modes=(0,))
    np.testing.assert_allclose(solved[0], table[:, 1], rtol=1e-6)


def write_synthetic_run(root, g_true, *, g_ref=0.7, sign=-1, transfer_direction=None):
    """Noise-free saved pair whose wavefield is the water-model fundamental mode."""
    layers = np.loadtxt(
        REFERENCE / "rr52_historical_rayleigh_layers.csv", delimiter=",", skiprows=1
    )
    f = np.arange(1, 241) / 600.0 + 0.01
    c = modal_phase_velocity(layers, f, modes=(0,), water=(DEPTH, SPEED, DENSITY))[0]
    height = admittance_height(f, c, DEPTH, SPEED)
    rng = np.random.default_rng(4)
    windows = 6
    p_true = rng.normal(size=(windows, len(f))) + 1j * rng.normal(
        size=(windows, len(f))
    )
    d = -p_true / ((2 * np.pi * f) ** 2 * DENSITY * height)
    p_nominal = p_true / (g_true * sign)
    pp, dd, dp = abs(p_nominal) ** 2, abs(d) ** 2, d * np.conj(p_nominal)
    epoch = dict(
        Seed_id="XX.SYN.00.BDH",
        channel_start_time="2013-01-01T00:00:00.000000Z",
        channel_end_time="2014-01-01T00:00:00.000000Z",
        inventory_file="inv.xml",
    )
    vertical = {**epoch, "Seed_id": "XX.SYN.00.BHZ"}
    run = root / "paired_spectra" / "run_synthetic"
    run.mkdir(parents=True)
    np.savez(
        run / "spectra.npz",
        frequency_hz=f,
        pressure_psd=pp.mean(0),
        displacement_psd=dd.mean(0),
        cross_dp=dp.mean(0),
        transfer_dp=dp.mean(0) / pp.mean(0),
        coherence_squared=abs(dp.mean(0)) ** 2 / (pp.mean(0) * dd.mean(0)),
        valid=np.ones(len(f), bool),
        pressure_window_psd=pp,
        displacement_window_psd=dd,
        cross_window_dp=dp,
        pressure_response=np.ones(len(f), complex),
        displacement_response=np.ones(len(f), complex),
        retained_window_starts=np.array(
            [f"2013-06-1{i // 3}T0{i % 3}:00:00.000000Z" for i in range(windows)]
        ),
    )
    (run / "comparison.json").write_text(
        json.dumps(
            {
                "cross_convention": "S_dp = mean(D * conjugate(P))",
                "transfer_direction": transfer_direction
                or "pressure_to_upward_displacement",
                "transfer_units": "m/Pa",
                "pressure_psd_units": "Pa^2/Hz",
                "vertical_psd_units": "m^2/Hz",
                "input_units": {"pressure": "PA", "vertical": "M/S"},
                "retained_segments": windows,
                "pressure_epoch": epoch,
                "vertical_epoch": vertical,
                "vertical_seed_id": vertical["Seed_id"],
                "start_utc": "2013-06-10T00:00:00Z",
                "end_utc": "2013-06-12T00:00:00Z",
            }
        )
    )
    reference = root / "reference.json"
    reference.write_text(
        json.dumps(
            {
                "station": "XX.SYN",
                "epoch_keys": [epoch, vertical],
                "median_pressure_multiplier": g_ref,
                "pressure_sign_assumption": {"status": "assumed_needs_confirmation"},
                "settings": {
                    "depth_m": DEPTH,
                    "density_kg_m3": DENSITY,
                    "sound_speed_m_s": SPEED,
                    "pressure_sign": sign,
                },
            }
        )
    )
    return run, reference


def digest(paths):
    return {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def test_run_check_recovers_synthetic_gain_and_preserves_inputs(tmp_path):
    pytest.importorskip("disba")
    run, reference = write_synthetic_run(tmp_path / "saved", g_true=0.8, g_ref=0.8)
    layers = REFERENCE / "rr52_historical_rayleigh_layers.csv"
    protected = [*run.iterdir(), reference]
    before = digest(protected)
    result = run_check(run, reference, layers, tmp_path / "out")
    assert digest(protected) == before
    assert sorted(p.name for p in run.iterdir()) == ["comparison.json", "spectra.npz"]
    band = result["bands"]["long_wave"]
    model = band["model_conditional_multiplier"]
    assert model["mode0_water_h1"]["median"] == pytest.approx(0.8, rel=1e-6)
    assert model["mode0_water_h2"]["median"] == pytest.approx(0.8, rel=1e-6)
    assert band["fraction_phase_within_90deg_of_expected"] == 1.0
    bound = band["structure_free_upper_bound"]
    assert bound["fraction_reference_consistent"] == 1.0
    assert bound["median_upper_bound_h1"] >= 0.8
    for sub in result["long_wave_sub_bands"]:
        assert sub["mode0_water_h1"] == pytest.approx(0.8, rel=1e-6)
        assert all(
            v[0] == pytest.approx(0.8) for v in sub["daily_mode0_water_h1_h2"].values()
        )
    assert result["gain_applied"] is False and result["inputs_unchanged"] is True
    summary = Path(result["files"]["summary"])
    assert summary.is_file() and Path(result["files"]["plot"]).is_file()
    assert json.loads(summary.read_text())["reference"]["pressure_multiplier"] == 0.8
    c = band["median_apparent_phase_velocity_m_s"]["h1"]
    assert c == pytest.approx(
        band["median_model_phase_velocity_m_s"]["mode0_water"], rel=1e-6
    )


def test_run_check_flags_reference_above_structure_free_bound(tmp_path):
    pytest.importorskip("disba")
    run, reference = write_synthetic_run(tmp_path / "saved", g_true=0.5)
    layers = REFERENCE / "rr52_historical_rayleigh_layers.csv"
    result = run_check(run, reference, layers, tmp_path / "out")
    bound = result["bands"]["long_wave"]["structure_free_upper_bound"]
    # The vertical-incidence bound loosens towards 0.9 f_r; low bins still refuse.
    assert bound["fraction_reference_violates"] > 0.5
    assert bound["median_upper_bound_h2"] < 0.7


def test_run_check_refuses_convention_drift_and_output_beside_inputs(tmp_path):
    pytest.importorskip("disba")
    run, reference = write_synthetic_run(
        tmp_path / "saved", g_true=0.8, transfer_direction="displacement_to_pressure"
    )
    layers = REFERENCE / "rr52_historical_rayleigh_layers.csv"
    with pytest.raises(ValueError, match="convention"):
        run_check(run, reference, layers, tmp_path / "out")
    with pytest.raises(ValueError, match="beside"):
        run_check(run, reference, layers, run.parent / "check")
    with pytest.raises(ValueError, match="hash"):
        run_check(run, reference, layers, tmp_path / "out", reference_sha256="0" * 64)


def test_run_check_without_reference_uses_caller_water(tmp_path):
    pytest.importorskip("disba")
    run, _ = write_synthetic_run(tmp_path / "saved", g_true=0.8)
    layers = REFERENCE / "rr52_historical_rayleigh_layers.csv"
    water = (DEPTH, DENSITY, SPEED)
    result = run_check(run, None, layers, tmp_path / "out", water=water)
    band = result["bands"]["long_wave"]
    assert result["reference"] is None
    assert band["model_conditional_multiplier"]["mode0_water_h1"]["median"] == (
        pytest.approx(0.8, rel=1e-6)
    )
    assert band["fraction_phase_within_90deg_of_zero"] == 1.0
    assert "fraction_phase_within_90deg_of_expected" not in band
    assert "fraction_reference_consistent" not in band["structure_free_upper_bound"]
    assert band["median_apparent_phase_velocity_m_s"] == {}
    with pytest.raises(ValueError, match="exactly one"):
        run_check(run, None, layers, tmp_path / "out")


def test_phase_screen_separates_pole_polarity_and_unexplained():
    pytest.importorskip("scripts.summarize_yv_microseism")
    from scripts.summarize_yv_microseism import phase_screen, reflected_pole_phase

    f = np.linspace(0.06, 0.3, 200)
    coherence = np.full_like(f, 0.9)
    pole = [0.1]
    cases = {
        "nominal": (np.zeros_like(f), []),
        "reflected_pole": (reflected_pole_phase(f, pole), pole),
        "reversed_reflected_pole": (
            (180 + reflected_pole_phase(f, pole) + 180) % 360 - 180,
            pole,
        ),
        "unexplained_phase": (np.full_like(f, 100.0), pole),
    }
    for label, (phase, poles) in cases.items():
        assert phase_screen(f, phase, coherence, poles, 0.08)["label"] == label
    sparse = phase_screen(f, np.zeros_like(f), np.full_like(f, 0.1), [], 0.08)
    assert sparse["label"] == "insufficient_coherent_bins"
    # Without a recorded pole, pole hypotheses are not offered.
    result = phase_screen(f, np.zeros_like(f), coherence, [], 0.08)
    assert set(result["rms_residual_deg"]) == {"nominal", "reversed_nominal"}


def test_positive_real_poles_read_from_snapshot_in_rad_per_second(tmp_path):
    from obspy import UTCDateTime
    from obspy.core.inventory import Channel, Inventory, Network, Station
    from obspy.core.inventory.response import Response

    pytest.importorskip("scripts.compare_microseism_methods")
    from scripts.compare_microseism_methods import positive_real_poles

    channels = []
    for code, kind in (("LDH", "LAPLACE (RADIANS/SECOND)"), ("BDH", "LAPLACE (HERTZ)")):
        response = Response.from_paz(
            zeros=[], poles=[0.1 + 0j, -0.2 + 0j, 0.3 + 0.4j], stage_gain=1
        )
        response.response_stages[0].pz_transfer_function_type = kind
        channels.append(
            Channel(
                code,
                "",
                0,
                0,
                -1000,
                0,
                start_date=UTCDateTime(2020, 1, 1),
                response=response,
            )
        )
    path = tmp_path / "inventory.xml"
    Inventory(
        [Network("XX", stations=[Station("S", 0, 0, 0, channels=channels)])]
    ).write(str(path), format="STATIONXML")
    assert positive_real_poles(path, "XX.S..LDH") == [pytest.approx(0.1)]
    assert positive_real_poles(path, "XX.S..BDH") == [pytest.approx(0.2 * np.pi)]


def test_vertical_pole_enters_with_opposite_sign():
    pytest.importorskip("scripts.summarize_yv_microseism")
    from scripts.summarize_yv_microseism import phase_screen, reflected_pole_phase

    f = np.linspace(0.06, 0.3, 200)
    phase = -reflected_pole_phase(f, [0.1])  # vertical pole only
    result = phase_screen(f, phase, np.full_like(f, 0.9), [], 0.08, [0.1])
    assert result["label"] == "reflected_pole"
