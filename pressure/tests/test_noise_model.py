"""Test model statistics and rejection paths without approving real instruments."""

from types import SimpleNamespace

import numpy as np
import pytest
from obspy import UTCDateTime
from test_quality_review import adopted_review, create_review_case, run_review

from compy_pressure import noise_model as noise
from compy_pressure.plot_network_ppsd import accumulate_histogram


@pytest.fixture
def review_case(tmp_path):
    return create_review_case(tmp_path)


def state(level, n, edges=(0.5, 1.5)):
    result = {}
    accumulate_histogram(
        result, np.array([1.0]), np.full((n, 1), level), np.array(edges)
    )
    return result


def fake_ppsd():
    periods = np.array([0.5, 1.0, 2.0, 10.0])
    return SimpleNamespace(
        id="YV.RR52.00.BDH",
        special_handling="hydrophone",
        period_bin_centers=periods,
        period_bin_left_edges=periods / np.sqrt(2),
        period_bin_right_edges=periods * np.sqrt(2),
        period_xedges=np.array([0.35, 0.7, 1.4, 4.5, 14.0]),
        psd_values=np.array([[-20.0, -10.0, 0.0, 10.0], [-30.0, -20.0, -10.0, 0.0]]),
        times_processed=[UTCDateTime("2020-01-01"), UTCDateTime("2020-01-01") + 20],
        ppsd_length=10.0,
    )


def test_station_balancing_does_not_weight_long_records_more():
    _, _, db, probability, coverage = noise.station_distribution(
        [state(0, 100), state(10, 1)]
    )
    assert coverage.tolist() == [2]
    assert probability[0, 0] == 0.5
    assert probability[-1, 0] == 0.5
    curves = noise.quantile_curves(probability, db, coverage, (5, 95))
    np.testing.assert_array_equal(curves[:, 0], [0.5, 10.5])


def test_station_distribution_keeps_frequency_coverage_explicit():
    _, _, _, probability, coverage = noise.station_distribution(
        [state(0, 10), state(10, 1, (0.75, 1.75))]
    )
    assert coverage.tolist() == [1, 2, 1]
    np.testing.assert_allclose(probability.sum(axis=0), 1)


@pytest.mark.parametrize("percentiles", [(0, 95), (5, 100), (95, 5), (5, np.nan)])
def test_invalid_percentiles_refused(percentiles):
    with pytest.raises(ValueError, match="percentiles"):
        noise.quantile_curves(
            np.array([[1.0]]), np.array([0, 1]), np.array([1]), percentiles
        )


def test_actual_review_plan_does_not_manufacture_accepted_stations(review_case):
    root, _, _ = review_case
    run_review(review_case)
    path = root / "results/REVIEWS/review/review.json"
    plan = noise.build_noise_model(
        path, root / "results/REVIEWS/plan", root=root, plan=True
    )
    assert plan["channels_approved"] == 0
    assert not plan["models"]
    with pytest.raises(ValueError, match="No approved channels"):
        noise.build_noise_model(path, root / "results/REVIEWS/model", root=root)
    assert not (root / "results/REVIEWS/model").exists()


def test_approved_synthetic_pipeline_checks_sources_and_produces_model(
    review_case, monkeypatch
):
    root, _, _ = review_case
    adopted_review(review_case)
    monkeypatch.setattr(noise.PPSD, "load_npz", lambda _: fake_ppsd())
    path = root / "results/REVIEWS/adopted/review.json"
    before = (root / "master.csv").read_bytes()
    report = noise.build_noise_model(path, root / "results/REVIEWS/model", root=root)
    assert report["channels_approved"] == 1
    assert report["models"][0]["retained_windows"] == 2
    assert not report["gain_applied"]
    assert (root / "master.csv").read_bytes() == before
    assert (root / "results/REVIEWS/model/YV_BDH_noise_curves.png").is_file()
    (root / "ppsd.npz").write_bytes(b"changed input")
    with pytest.raises(ValueError, match="hash changed"):
        noise.build_noise_model(path, root / "results/REVIEWS/stale", root=root)


def test_full_time_windows_and_smoothing_support_are_required(review_case):
    report = adopted_review(review_case)
    row = report["rows"][0]
    row["review_decision"]["approved_frequency_band_hz"] = [0.2, 2]
    row["review_decision"]["approved_time_range_utc"] = [
        "2020-01-01T00:00:05Z",
        "2020-01-01T00:00:30Z",
    ]
    f, values, _, _ = noise.selected_samples(fake_ppsd(), row)
    np.testing.assert_array_equal(f, [0.5, 1])
    assert values.shape == (1, 2)
