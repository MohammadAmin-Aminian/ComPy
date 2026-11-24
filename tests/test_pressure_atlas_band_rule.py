"""Tests of compy_pressure.atlas.band_rule: the two limits per window of rule W."""

import numpy as np
import pandas as pd
import pytest

from compy_pressure.atlas import band_rule, config, native

GRID = np.asarray(native.GRID)


def _highpass_db(tau_s):
    """10 log10 |H|^2 of a first-order high-pass of time constant tau_s on the grid."""
    x = (2 * np.pi / GRID) * tau_s
    return 10 * np.log10(x**2 / (1 + x**2))


def test_welch_limit_is_one_welch_bin_in_a_quarter_octave():
    f = band_rule.welch_limit_hz(2700.0)
    width = f * (2**0.125 - 2**-0.125)
    assert width == pytest.approx(1 / 2700.0, rel=1e-12)
    assert 1 / f == pytest.approx(468.5, abs=0.5)


def test_response_band_of_a_first_order_high_pass_is_its_corner():
    tau = 80.0
    f_lo, f_hi = band_rule.response_band_hz(_highpass_db(tau))
    # -3 dB of a first-order high-pass: f = 1 / (2 pi tau); the grid is 1/8 octave, so allow the interpolation error
    assert f_lo == pytest.approx(1 / (2 * np.pi * tau), rel=0.03)
    assert np.isnan(f_hi)  # flat to the shortest period of the grid


def test_response_band_of_a_flat_response_has_no_limit():
    f_lo, f_hi = band_rule.response_band_hz(np.zeros(len(GRID)))
    assert np.isnan(f_lo) and np.isnan(f_hi)


def test_response_band_short_period_side():
    # a first-order low-pass with corner at 2 Hz: -3 dB at 2 Hz
    x = (1 / GRID) / 2.0
    db = -10 * np.log10(1 + x**2)
    f_lo, f_hi = band_rule.response_band_hz(db)
    assert f_hi == pytest.approx(2.0, rel=0.03)
    assert np.isnan(f_lo)


def test_empirical_cap_is_where_a_type_leaves_the_reference():
    P = GRID
    n = 12
    ref = np.tile(
        -30 + 20 * np.log10(np.maximum(P, 0.3) / 0.3) * 0, (n, 1)
    )  # flat reference
    apg = np.tile(
        np.where(P < 1.0, -30 + 25 * (1 - P), -30.0), (n, 1)
    )  # +25 dB at 0 s, 0 dB at 1 s
    med = np.vstack([ref, apg])
    typ = ["Hydrophone"] * n + ["APG"] * n
    caps = band_rule.empirical_caps(med, P, typ, np.full(2 * n, 3000.0))
    assert caps["Hydrophone"] == np.inf
    # the excess exceeds 6 dB where 25 (1 - P) > 6, i.e. P < 0.76 s: the cap is the first node at or above 0.76 s
    node = P[P >= 0.76].min()
    assert caps["APG"] == pytest.approx(1 / node, rel=1e-9)


def test_empirical_cap_ignores_shallow_gauges_and_small_groups():
    P = GRID
    med = np.tile(-30.0, (9, len(P)))
    typ = ["Hydrophone"] * 8 + ["DPG"]
    caps = band_rule.empirical_caps(med, P, typ, np.full(9, 3000.0))
    assert np.isnan(caps["DPG"])  # one gauge: fewer than 8


def test_channel_limits_dpg_and_apg_and_hydrophone():
    flat = np.zeros(len(GRID))
    welch = band_rule.welch_limit_hz()
    # DPG with a 80-s high-pass: response corner at 0.002 Hz is under 1/300 s: the 300 s limit sets the low limit
    d = band_rule.channel_limits("DPG", 62.5, _highpass_db(80.0))
    assert d.lo_by == ["DPG time constant"]
    assert d.f_lo_hz[0] == pytest.approx(1 / 300.0)
    assert d.f_hi_hz[0] == pytest.approx(
        13.0
    )  # the self-noise cap, under 0.4 * 62.5 / 2**0.125 = 22.9
    assert d.hi_by == ["self-noise"]
    # a DPG with a short time constant: the response sets the limit
    d2 = band_rule.channel_limits("DPG", 62.5, _highpass_db(18.0))
    assert d2.lo_by == ["response"] and d2.f_lo_hz[0] > 1 / 300.0
    # APG: flat to DC, so only the estimator sets the low limit; 0.745 Hz cap
    a = band_rule.channel_limits("APG", 125.0, flat)
    assert a.lo_by == ["estimator"] and a.f_lo_hz[0] == pytest.approx(welch)
    assert a.f_hi_hz[0] == pytest.approx(0.745)
    # hydrophone: no cap, the sampling rate sets the high limit; the validated cut-off replaces the response
    h = band_rule.channel_limits(
        "Hydrophone", 100.0, _highpass_db(0.1), validated_pz_hz=0.05
    )
    assert h.lo_by == ["validated P-Z phase"] and h.f_lo_hz[0] == pytest.approx(0.05)
    assert h.hi_by == ["sampling"] and h.f_hi_hz[0] == pytest.approx(
        0.4 * 100 / 2**0.125
    )
    h2 = band_rule.channel_limits("Hydrophone", 100.0, _highpass_db(0.1))
    assert h2.lo_by == ["response"]


def test_channel_limits_take_the_response_of_each_epoch():
    db = np.vstack([_highpass_db(18.0), _highpass_db(80.0)])
    lim = band_rule.channel_limits("DPG", 50.0, db)
    assert len(lim.f_lo_hz) == 2 and lim.f_lo_hz[0] > lim.f_lo_hz[1]
    f_lo, f_hi = band_rule.window_limits(lim, [0, 1, 1], [True, True, False])
    assert f_lo[0] == lim.f_lo_hz[0] and f_lo[1] == lim.f_lo_hz[1]
    assert (
        f_lo[2] >= 1 / 300.0
    )  # a window without long context has no level beyond 300 s


def test_band_mask_keeps_bins_with_centre_in_the_band():
    f_lo, f_hi = np.array([0.01, 1.0]), np.array([1.0, 2.0])
    m = band_rule.band_mask(f_lo, f_hi)
    f = 1 / GRID
    assert m.shape == (2, len(GRID))
    assert np.all((f[m[0]] >= 0.01) & (f[m[0]] <= 1.0))
    assert np.all((f[m[1]] >= 1.0) & (f[m[1]] <= 2.0))
    assert m[0].any() and m[1].any()


def test_add_limits_reads_the_response_of_each_native_file(tmp_path):
    for name, db in (("dpg", _highpass_db(80.0)), ("hyd", _highpass_db(0.1))):
        np.savez(tmp_path / f"{name}.npz", periods=GRID, response_db=db[None, :])
    t = pd.DataFrame(
        dict(
            seed_id=["A", "B", "C"],
            type=["DPG", "Hydrophone", "APG"],
            fs=[62.5, 100.0, 125.0],
            f_lo_hz=[0.003, 0.05, 0.0],
            f_lo_source=[
                "response 20 dB limit",
                "validated hydrophone cut-off (supplement)",
                "",
            ],
            f_hi_hz=[25.0, 40.0, 60.0],
            product=[
                str(tmp_path / "dpg.npz"),
                str(tmp_path / "hyd.npz"),
                str(tmp_path / "missing.npz"),
            ],
        )
    )
    out = band_rule.add_limits(t)
    assert (
        out.f_lo_hz[0] == pytest.approx(1 / 300.0)
        and out.band_rule_lo_by[0] == "DPG time constant"
    )
    assert (
        out.f_lo_hz[1] == pytest.approx(0.05)
        and out.band_rule_lo_by[1] == "validated P-Z phase"
    )
    assert np.isnan(out.f_lo_hz[2]) and out.band_rule_hi_by[2] == "no native file"
    assert list(out.f_lo_hz_before) == [0.003, 0.05, 0.0]
    with pytest.raises(ValueError):
        band_rule.add_limits(t.drop(columns="product"))


def test_config_selects_the_rule_and_defaults_keep_the_atlas_band():
    cfg = config.load()
    assert cfg["population"]["band_rule"] == "atlas"
    assert config.load(population={"band_rule": "W"})["population"]["band_rule"] == "W"


def test_band_taper_has_hard_edges_for_rule_w_and_the_octave_taper_otherwise():
    from compy_pressure.atlas import models

    f = 1 / GRID
    hard = models.band_taper(0.01, 1.0, GRID, edge_octaves=0.0)
    assert set(np.unique(hard)) == {0.0, 1.0}
    assert np.all(hard[(f >= 0.01) & (f <= 1.0)] == 1.0) and np.all(
        hard[(f < 0.01) | (f > 1.0)] == 0.0
    )
    soft = models.band_taper(
        0.01, 1.0, GRID
    )  # default: unchanged behaviour of the November 2026 atlas
    inside_edge = (f >= 0.01) & (f < 0.02)
    assert np.all(soft[inside_edge] < 1.0) and np.all(
        soft[(f >= 0.02) & (f <= 0.5)] == 1.0
    )


# ---------------------------------------------------------------------------------------------- rule V
def _two_pole_db(t_short, t_long):
    return -(band_rule.correction_db(GRID, t_short) + band_rule.correction_db(GRID, t_long))


def test_dpg_poles_find_one_and_two_poles():
    one = band_rule.dpg_poles_s(_highpass_db(80.0))
    assert one[0] == pytest.approx(2 * np.pi * 80.0, rel=0.1) and np.isnan(one[1])
    ts, tl = band_rule.dpg_poles_s(_two_pole_db(30.0, 600.0))
    assert ts == pytest.approx(30.0, rel=0.15) and tl == pytest.approx(600.0, rel=0.25)


def test_tolerance_period_grows_with_the_corner_and_is_infinite_for_small_errors():
    t3 = band_rule.tolerance_period_s(500.0, 2.5, 3.0)
    assert 450.0 < t3 < 600.0  # about the corner: the error is 2.5 dB there for k = 3
    assert band_rule.tolerance_period_s(30.0, 2.5, 3.0) < 35.0
    assert np.isinf(band_rule.tolerance_period_s(30.0, 2.5, 1.2, two_sided=True))  # 20 log10 1.2 = 1.6 dB < 2.5 dB


def test_rule_v_dpg_limits_verified_and_unverified():
    cfg = {"rule": "V"}
    db = _two_pole_db(29.5, 595.0)[None, :]
    plain = band_rule.channel_limits("DPG", 40.0, db, cfg=cfg)
    assert 1 / plain.f_lo_hz[0] < 40.0 and plain.lo_by == ["response"]  # a 30 s corner nobody checked: its -3 dB point, as in rule W
    ok = band_rule.channel_limits("DPG", 40.0, db, cfg=cfg, dpg_validated=True)
    assert 1 / ok.f_lo_hz[0] == pytest.approx(468.5, abs=1.0) and ok.lo_by == ["estimator"]  # verified: up to the estimator limit
    slow = band_rule.channel_limits("DPG", 40.0, _highpass_db(80.0), cfg=cfg)
    assert 1 / slow.f_lo_hz[0] == pytest.approx(300.0, abs=1.0)  # 80 s time constant: unverified, 300 s as in rule W


def test_rule_v_floor_limit_apg_depth_and_hydrophone_correction_cap():
    cfg = {"rule": "V"}
    flat = np.zeros(len(GRID))
    d = band_rule.channel_limits("DPG", 62.5, flat, cfg=cfg, f_hi_data_hz=3.0)
    assert d.f_hi_hz[0] == pytest.approx(3.0) and d.hi_by == ["noise floor"]
    assert band_rule.channel_limits("DPG", 62.5, flat, cfg=cfg).hi_by == ["sampling"]  # no 13 Hz type cap in rule V
    shallow = band_rule.channel_limits("APG", 125.0, flat, cfg=cfg, depth_m=175.0)
    deep = band_rule.channel_limits("APG", 125.0, flat, cfg=cfg, depth_m=3000.0)
    assert 1 / shallow.f_lo_hz[0] == pytest.approx(468.5, abs=1.0) and 1 / deep.f_lo_hz[0] == pytest.approx(468.5, abs=1.0)  # the 900 s option is off
    opt = band_rule.channel_limits("APG", 125.0, flat, cfg={"rule": "V", "apg_deep_depth_m": 1500.0}, depth_m=3000.0)
    assert 1 / opt.f_lo_hz[0] == pytest.approx(900.0, abs=1.0)
    # a hydrophone whose validated cut-off (200 s) lies beyond the period where the correction reaches 10 dB is cut there
    hp = _highpass_db(2.0)  # corner at 12.6 s: the correction reaches 10 dB at 3 x 12.6 = 38 s
    h = band_rule.channel_limits("Hydrophone", 100.0, hp, cfg=cfg, validated_pz_hz=1 / 200.0)
    assert h.lo_by == ["response correction cap"] and 30.0 < 1 / h.f_lo_hz[0] < 45.0
