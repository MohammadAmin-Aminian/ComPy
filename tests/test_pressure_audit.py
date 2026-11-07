"""Pressure-gauge audit (November 2026): water-layer model, table scan, figures."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from compy_pressure import audit_plots, table_audit, water_layer  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "pressure/docs/audit_2026-11/data"
DISPERSION = (
    ROOT / "compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv"
)
F = np.arange(0.03, 0.0701, 0.0025)
C = np.full(F.shape, 3500.0)


def test_ratio_is_the_zha_webb_expression():
    x = 2 * np.pi * F * 4000 * np.sqrt(1500.0**-2 - 3500.0**-2)
    np.testing.assert_allclose(
        water_layer.ratio(F, C, 4000), 1028 * 4000 * np.tan(x) / x, rtol=1e-12
    )
    slow = water_layer.ratio(F, 1000.0, 4000)  # c < alpha: tanh branch
    y = 2 * np.pi * F * 4000 * np.sqrt(1000.0**-2 - 1500.0**-2)
    np.testing.assert_allclose(slow, 1028 * 4000 * np.tanh(y) / y, rtol=1e-12)


@pytest.mark.parametrize("c", [3500.0, 1000.0])
def test_log_sensitivities_match_finite_differences(c):
    h, a, eps = 4500.0, 1500.0, 1e-6
    s = water_layer.log_sensitivity(F, c, h, a)

    def ln(**kw):
        return np.log(
            water_layer.ratio(F, kw.get("c", c), kw.get("h", h), 1028, kw.get("a", a))
        )

    np.testing.assert_allclose(
        s["sound_speed"],
        (ln(a=a * (1 + eps)) - ln(a=a * (1 - eps))) / (2 * eps),
        rtol=1e-5,
    )
    np.testing.assert_allclose(
        s["phase_velocity"],
        (ln(c=c * (1 + eps)) - ln(c=c * (1 - eps))) / (2 * eps),
        rtol=1e-5,
    )
    np.testing.assert_allclose(
        s["depth"], (ln(h=h * (1 + eps)) - ln(h=h * (1 - eps))) / (2 * eps), rtol=1e-5
    )


def test_guard_matches_rayleigh_calibration_rule():
    x = 2 * np.pi * F * 5600 * np.sqrt(1500.0**-2 - 3500.0**-2)
    assert np.array_equal(water_layer.usable(F, C, 5600), x < 0.9 * np.pi / 2)
    assert not water_layer.usable(F, C, 5600).all()


def test_gain_bias_is_one_for_the_true_model_and_follows_density():
    assert water_layer.gain_bias(F, C, 4000) == pytest.approx(1.0, abs=1e-12)
    assert water_layer.gain_bias(F, C, 4000, true_density=1038.28) == pytest.approx(
        1028 / 1038.28, rel=1e-12
    )
    up = water_layer.gain_bias(F, C, 5000, true_sound_speed=1515.0)
    log = water_layer.gain_bias(F, C, 5000, true_sound_speed=1515.0, weights="log")
    assert up > log > 1.0  # least squares weights the near-resonance bins most
    assert water_layer.column_mean_density(0) == pytest.approx(1027.8)
    assert 1037 < water_layer.column_mean_density(4500) < 1039


def test_fit_gain_matches_rayleigh_calibration():
    pytest.importorskip("obspy")
    from compy_pressure.rayleigh_calibration import fit_multiplier, water_ratio

    np.testing.assert_allclose(
        water_ratio(F, C, 4000, 1028, 1500), water_layer.ratio(F, C, 4000), rtol=1e-12
    )
    observed = water_layer.ratio(F, C, 4000) * (1 + 0.05 * np.sin(F * 300)) / 0.8
    expected = water_layer.ratio(F, C, 4000)
    result = fit_multiplier(observed, expected)[0]
    assert result["pressure_multiplier"] == pytest.approx(
        water_layer.fit_gain(observed, expected), rel=1e-12
    )


# --------------------------------------------------------------------------
def _table(tmp_path, rows):
    cols = sorted({k for r in rows for k in r} | {"Seed_id", "channel_start_time", "type", "basin",
                                                  "depth_bin", "final_class"})  # fmt: skip
    path = tmp_path / "t.csv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    return table_audit.load_table(path)


def test_scan_finds_outliers_disagreement_polarity_and_confounding(tmp_path):
    rows = []
    for i in range(8):  # deep network, low gains, no within-network depth trend
        rows.append(dict(Seed_id=f"AA.S{i}..BDH", type="DPG", basin="Indian", Profondeur=4000 + 50 * (i % 3),
                         microseism_g=0.7 * (1 + 0.02 * (i % 2)), microseism_coherence_squared=0.9,
                         rayleigh_g=0.7, rayleigh_sign=-1, local_p_g=0.7, local_p_polarity=-1))  # fmt: skip
    for i in range(8):  # shallow network, gains near one
        rows.append(dict(Seed_id=f"BB.S{i}..BDH", type="DPG", basin="Pacific", Profondeur=1000 + 50 * (i % 3),
                         microseism_g=1.0 * (1 + 0.02 * (i % 2)), microseism_coherence_squared=0.9,
                         rayleigh_g=1.0, rayleigh_sign=1 if i else -1, local_p_g=1.0,
                         local_p_polarity=1))  # fmt: skip
    rows[3]["microseism_g"] = 7.0  # outlier and disagreement with Rayleigh
    rows[9]["local_p_polarity"] = -1  # sign disagrees with Rayleigh +1 on BB.S1
    for i in range(6):
        rows.append(dict(Seed_id=f"CC.S{i}..BDH", type="DPG", basin="Indian", Profondeur=2000 + 600 * i,
                         microseism_g=0.7 + 0.02 * i, microseism_coherence_squared=0.9))  # fmt: skip
    t = _table(tmp_path, rows)
    found, info = table_audit.scan_table(t)
    kinds = {(a.kind, a.station) for a in found}
    assert ("gain_outlier", "AA.S3..BDH") in kinds
    assert ("method_disagreement", "AA.S3..BDH") in kinds
    assert ("polarity", "BB.S1..BDH") in kinds
    assert ("polarity", "BB") in kinds  # Rayleigh sign mixed within BB
    assert info["depth"]["microseism"]["across_p"] < 0.01
    assert abs(info["depth"]["microseism"]["within_slope_db_per_km"]) < 1.0


def test_invalid_offsets_are_flagged_and_not_used_as_gains(tmp_path):
    t = _table(
        tmp_path, [dict(Seed_id="AA.S1..BDH", type="DPG", final_class="nonworking")]
    )
    t["offset_db"][:] = -3087.5
    t["peer_multiplier"][:] = 2.4e154
    assert np.isnan(table_audit.gains(t)["peer_offset"]).all()
    assert [a.kind for a in table_audit.offset_invalid(t)] == ["offset_invalid"]


@pytest.mark.skipif(
    not (DATA / "pressure_qc_calibration_table.csv").exists(),
    reason="audit data absent",
)
def test_archived_master_table_scan():
    t = table_audit.load_table(
        DATA / "pressure_qc_calibration_table.csv",
        DATA / "metadata_offsets_channels.csv",
    )
    assert len(t["Seed_id"]) == 846
    found, info = table_audit.scan_table(t)
    assert info["n_with_gain"]["rayleigh"] == 54
    assert sum(a.kind == "offset_invalid" for a in found) == 5
    assert any(a.kind == "network_gain_structure" for a in found)
    for fig in (audit_plots.plot_gain_consistency(t), audit_plots.plot_gain_depth(t),
                audit_plots.plot_infragravity_residuals(t), audit_plots.plot_network_scores(t)):  # fmt: skip
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


def test_water_layer_figure_renders():
    d = np.loadtxt(DISPERSION, delimiter=",", skiprows=1)
    fig = audit_plots.plot_water_layer_sensitivity(d[:, 0], d[:, 1], depths=np.linspace(1000, 5500, 10),
                                                   freqs=np.linspace(0.02, 0.08, 13))  # fmt: skip
    assert isinstance(fig, plt.Figure)
    plt.close(fig)
