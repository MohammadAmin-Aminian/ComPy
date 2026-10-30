"""Read-only network PPSD plots: no raw processing, mixed units or overwrites."""

import csv
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pytest

from compy_pressure import plot_network_ppsd as plot


def fake_ppsd(seed="YV.RR52.00.BDH", pressure=True):
    return SimpleNamespace(
        id=seed,
        special_handling="hydrophone" if pressure else None,
        period_bin_centers=np.array([0.5, 1, 2, 10]),
        period_xedges=np.array([0.35, 0.7, 1.4, 4.5, 14.0]),
        psd_values=np.array([[-20, -10, 0, 10], [-30, -20, -10, 0]]),
        times_processed=["2020-01-01", "2020-01-02"],
    )


def row(path, seed="YV.RR52.00.BDH", pressure=True):
    return {
        "Seed_id": seed,
        "channel_start_time": "2020-01-01",
        "channel_end_time": "2021-01-01",
        "inventory_file": "inventory.xml",
        "ppsd_npz_file": str(path),
        "ppsd_status": "complete",
        "ppsd_units": "Pa^2/Hz" if pressure else "m^2/s^4/Hz",
    }


def master(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_frequency_order_and_exact_window_percentiles():
    f, values, pressure = plot.spectrum(fake_ppsd(), row("unused"))
    np.testing.assert_array_equal(f, [0.1, 0.5, 1, 2])
    np.testing.assert_array_equal(values[:, 0], [10, 0])
    np.testing.assert_array_equal(np.percentile(values, 50, axis=0), [5, -5, -15, -25])
    assert pressure


@pytest.mark.parametrize("fault", ["id", "units", "nan", "duplicate_times"])
def test_invalid_inputs_refused(fault):
    p, r = fake_ppsd(), row("unused")
    if fault == "id":
        p.id = "YV.OTHER.00.BDH"
    elif fault == "units":
        r["ppsd_units"] = "m^2/s^4/Hz"
    elif fault == "nan":
        p.psd_values = p.psd_values.astype(float)
        p.psd_values[0, 0] = np.nan
    else:
        p.times_processed = ["same", "same"]
    with pytest.raises(ValueError):
        plot.spectrum(p, r)


def test_discovery_deduplicates_paths_and_refuses_conflicting_epochs(tmp_path):
    r = row(tmp_path / "saved.npz")
    path = master(tmp_path / "master.csv", [r, r])
    assert len(plot.saved_rows(path, "YV", ["BDH"])) == 1
    master(path, [r, r | {"channel_start_time": "different"}])
    with pytest.raises(ValueError, match="conflicting"):
        plot.saved_rows(path, "YV", ["BDH"])


@pytest.mark.parametrize("pressure", [True, False])
def test_reference_units_kept_separate(pressure):
    fig, ax = plt.subplots()
    plot.add_references(ax, pressure, 0.001, 30)
    labels = [line.get_label() for line in ax.lines]
    if pressure:
        assert any("Brown" in label for label in labels)
        assert not any("Peterson" in label for label in labels)
        brown = [line for line in ax.lines if "Brown" in line.get_label()]
        assert all(min(line.get_xdata()) >= 0.2 for line in brown)
        assert all(min(line.get_xdata()) < 0.3 for line in brown)
        supplied = [line for line in ax.lines if "Webb 1998" in line.get_label()]
        assert len(supplied) == 2
        for line in supplied:
            np.testing.assert_allclose(line.get_xdata(), [0.005, 0.05, 0.1, 0.2])
        np.testing.assert_allclose(supplied[0].get_ydata(), [20, -25, -20, 10])
        np.testing.assert_allclose(supplied[1].get_ydata(), [50, 5, 20, 50])
    else:
        assert labels == ["Peterson NLNM", "Peterson NHNM"]
    plt.close(fig)


@pytest.mark.parametrize("all_curves", [False, True])
def test_plot_fresh_outputs_missing_files_and_protected_inputs(
    tmp_path, monkeypatch, all_curves
):
    saved = tmp_path / "saved.npz"
    saved.write_bytes(b"unchanged test input")
    rows = [row(saved), row(tmp_path / "missing.npz", "YV.RR40.00.BDH")]
    path = master(tmp_path / "master.csv", rows)
    before = path.read_bytes()
    monkeypatch.setattr(plot.PPSD, "load_npz", lambda _: fake_ppsd())
    out = tmp_path / "new_plots"
    report = plot.plot_network(
        path, channels=["BDH"], output=out, all_curves=all_curves
    )
    assert path.read_bytes() == before
    assert saved.read_bytes() == b"unchanged test input"
    assert report["protected_inputs_unchanged"]
    assert [r["plot_status"] for r in report["inputs"]] == ["missing_npz", "plotted"]
    assert len(report["figures"]) == 2
    assert (out / "YV_BDH_station_ppsd.png").is_file()
    assert (out / "YV_BDH_station_ppsd_full_range.png").is_file()
    assert (out / "plot_manifest.json").is_file()
    with pytest.raises(FileExistsError):
        plot.plot_network(path, channels=["BDH"], output=out)


def test_invalid_pressure_limits_refused_before_io():
    with pytest.raises(ValueError, match="Pressure limits"):
        plot.plot_network("unused", pressure_limits=(1, -1))


def test_histogram_pools_all_windows_and_preserves_extremes():
    state = {}
    f, edges = np.array([1.0, 2.0]), np.array([0.5, 1.5, 2.5])
    plot.accumulate_histogram(state, f, np.array([[-3077.0, 0.0], [5.0, 1.0]]), edges)
    plot.accumulate_histogram(state, f, np.array([[10.0, 2.0]]), edges)
    assert state["windows"] == 3
    np.testing.assert_array_equal(state["counts"].sum(axis=0), [3, 3])
    assert state["lo"] == -3077
    assert state["counts"][0, 0] == 1
    with pytest.raises(ValueError, match="Invalid native"):
        plot.accumulate_histogram(state, f * 2, np.array([[1.0, 2.0]]), edges)


def test_histogram_mixed_grids_preserves_frequency_dependent_coverage():
    state = {}
    plot.accumulate_histogram(
        state, np.array([1.0, 2.0]), np.array([[0.0, 5.0]]), np.array([0.5, 1.5, 2.5])
    )
    plot.accumulate_histogram(
        state, np.array([2.0, 3.0]), np.array([[10.0, 15.0]]), np.array([1.5, 2.5, 3.5])
    )
    np.testing.assert_array_equal(state["counts"].sum(axis=0), [1, 2, 1])
    assert state["windows"] == 2


def test_probability_has_no_summary_computation(tmp_path, monkeypatch):
    saved = tmp_path / "saved.npz"
    saved.write_bytes(b"unchanged")
    path = master(tmp_path / "master.csv", [row(saved)])
    monkeypatch.setattr(plot.PPSD, "load_npz", lambda _: fake_ppsd())

    def forbidden(*args, **kwargs):
        raise AssertionError("No percentile summaries in probability mode")

    monkeypatch.setattr(plot.np, "percentile", forbidden)
    report = plot.plot_network(
        path, channels=["BDH"], output=tmp_path / "heatmap", probability=True
    )
    assert report["mode"] == "probability_heatmap"
    assert report["percentiles"] == []
    assert report["single_color"]
    assert report["probability_display"]["colormap"] == "pqlx"
    assert (tmp_path / "heatmap/YV_BDH_probability_ppsd.png").is_file()


def test_display_threshold_does_not_change_histogram(tmp_path, monkeypatch):
    saved = tmp_path / "saved.npz"
    saved.write_bytes(b"unchanged")
    path = master(tmp_path / "master.csv", [row(saved)])
    monkeypatch.setattr(plot.PPSD, "load_npz", lambda _: fake_ppsd())
    for name, threshold in [("all", 0.0), ("display_cutoff", 0.1)]:
        report = plot.plot_network(
            path,
            channels=["BDH"],
            output=tmp_path / name,
            probability=True,
            min_percentage=threshold,
            max_percentage=15.0,
        )
        assert report["probability_display"]["minimum_shown_percent"] == threshold
    with (
        np.load(tmp_path / "all/YV_BDH_probability.npz") as a,
        np.load(tmp_path / "display_cutoff/YV_BDH_probability.npz") as b,
    ):
        for key in a.files:
            np.testing.assert_array_equal(a[key], b[key])


@pytest.mark.parametrize("limits", [(-1, 30), (30, 30), (0, 101), (0, np.nan)])
def test_invalid_probability_limits(limits):
    with pytest.raises(ValueError, match="occurrence"):
        plot.plot_network("unused", min_percentage=limits[0], max_percentage=limits[1])


def test_curves_only_draws_every_window_without_summary_or_station_colors(
    tmp_path, monkeypatch
):
    saved = tmp_path / "saved.npz"
    saved.write_bytes(b"unchanged")
    path = master(tmp_path / "master.csv", [row(saved)])
    monkeypatch.setattr(plot.PPSD, "load_npz", lambda _: fake_ppsd())

    def forbidden(*args, **kwargs):
        raise AssertionError("No percentiles may be computed in curves-only mode")

    monkeypatch.setattr(plot.np, "percentile", forbidden)
    original = plot.LineCollection
    collections = []

    def capture(segments, **kwargs):
        collections.append((segments.copy(), kwargs))
        return original(segments, **kwargs)

    monkeypatch.setattr(plot, "LineCollection", capture)
    report = plot.plot_network(
        path,
        channels=["BDH"],
        output=tmp_path / "curves_only",
        curves_only=True,
        single_color=True,
    )
    assert report["mode"] == "all_windows_only"
    assert report["percentiles"] == []
    assert report["single_color"]
    segments, options = collections[0]
    assert segments.shape == (2, 4, 2)
    np.testing.assert_array_equal(segments[:, :, 1], fake_ppsd().psd_values[:, ::-1])
    assert options["colors"] == ["#4477AA"]
