import numpy as np
import pandas as pd
import pytest

from compy_pressure import group_noise_model as gnm
from compy_pressure.plot_network_ppsd import accumulate_histogram


def test_local_season_flips_south_of_the_equator():
    assert [gnm.local_season(m, 45.0) for m in (1, 4, 7, 10)] == list(gnm.SEASONS)
    assert [gnm.local_season(m, -20.0) for m in (7, 10, 1, 4)] == list(gnm.SEASONS)
    assert gnm.local_season(12, 10.0) == "winter"
    assert gnm.local_season(12, -10.0) == "summer"


def test_one_per_station_keeps_widest_band_per_sensor_type():
    table = pd.DataFrame(
        {
            "Seed_id": ["XX.A..BDH", "XX.A..LDH", "XX.A..HXH", "XX.B..BDH"],
            "type": ["DPG", "DPG", "Hydrophone", "DPG"],
            "band_min_hz": [0.005, 0.005, 1.0, 0.01],
            "band_max_hz": [8.0, 0.4, 50.0, 8.0],
        }
    )
    kept = gnm.one_per_station(table)
    assert sorted(kept.Seed_id) == ["XX.A..BDH", "XX.A..HXH", "XX.B..BDH"]
    assert set(kept.station) == {"XX.A", "XX.B"}


class FakePPSD:
    """Overlapping one-octave smoothing bins on a 1/2-octave grid, as in ObsPy."""

    def __init__(self):
        periods = 2.0 ** np.arange(0, 6, 0.5)  # 1 s ... 45 s
        self.period_bin_centers = periods
        self.period_bin_left_edges = periods / np.sqrt(2)
        self.period_bin_right_edges = periods * np.sqrt(2)
        self.period_xedges = np.r_[periods / 2**0.25, periods[-1] * 2**0.25]


def test_band_selection_keeps_whole_bins_with_display_edges():
    f, edges, columns = gnm.band_selection(FakePPSD(), 0.05, 0.6)
    assert np.all(np.diff(f) > 0)
    assert len(edges) == len(f) + 1
    assert np.all((f > edges[:-1]) & (f < edges[1:]))
    # Each kept bin's full smoothing interval lies inside the band.
    periods = FakePPSD().period_bin_centers[columns]
    assert np.all(1 / (periods * np.sqrt(2)) >= 0.05)
    assert np.all(1 / (periods / np.sqrt(2)) <= 0.6)
    with pytest.raises(ValueError):
        gnm.band_selection(FakePPSD(), 10.0, 20.0)


def state(level, seed):
    f = np.array([0.1, 0.2, 0.4])
    edges = np.array([0.07, 0.14, 0.28, 0.56])
    values = level + np.random.default_rng(seed).normal(0, 1, (200, 3))
    s = {}
    accumulate_histogram(s, f, values, edges)
    return s


def test_build_models_weights_stations_equally_and_groups_by_column():
    results = [
        {"Seed_id": "XX.A..BDH", "states": {"all": state(0, 1), "winter": state(0, 2)}},
        {
            "Seed_id": "XX.B..BDH",
            "states": {"all": state(20, 3), "summer": state(20, 4)},
        },
    ]
    rows = pd.DataFrame(
        {"Seed_id": ["XX.A..BDH", "XX.B..BDH"], "basin": ["Atlantic", "Pacific"]}
    )
    everything = gnm.build_models(results, rows, "all")["all"]
    curves, coverage, stations = (
        everything.curves,
        everything.coverage,
        everything.stations,
    )
    assert stations == 2 and list(coverage) == [2, 2, 2]
    # Two equally weighted stations 20 dB apart: the median lies between them.
    assert np.all((curves[50] > 0) & (curves[50] < 21))
    assert curves[5][0] < 2 and curves[95][0] > 18
    basins = gnm.build_models(results, rows, "basin")
    assert list(basins) == ["Atlantic", "Pacific"]
    assert abs(basins["Pacific"].curves[50][0] - 20) <= 1.5
    seasons = gnm.build_models(results, rows, "season")
    assert list(seasons) == ["winter", "summer"]


def test_heatmap_and_line_figures_are_written(tmp_path):
    f = np.array([0.1, 0.2, 0.4])
    results = [
        {
            "Seed_id": f"XX.S{i}..BDH",
            "states": {"all": state(10 * i, i)},
            "medians": {"all": (f, np.full(3, 10.0 * i))},
        }
        for i in range(4)
    ]
    rows = pd.DataFrame({"Seed_id": [r["Seed_id"] for r in results]})
    models = gnm.build_models(results, rows, "all")
    assert len(models["all"].lines) == 4
    assert models["all"].probability.shape[1] == len(models["all"].frequency)
    gnm.plot_heatmaps(models, tmp_path / "heat.png", "test")
    gnm.plot_lines(models, tmp_path / "lines.png", "test")
    assert (tmp_path / "heat.png").stat().st_size > 0
    assert (tmp_path / "lines.png").stat().st_size > 0
