"""Network median PSD figure: limits, off-scale channels and exclusive output."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from compy_pressure.pressure_test import load_references  # noqa: E402
from scripts.plot_network_median_ppsd import (  # noqa: E402
    draw_pressure_limits,
    off_scale,
    plot,
)


def curve(station, level, pressure=True):
    f = np.logspace(-3, 1.4, 40)
    return {
        "station": station,
        "seed_id": f"YV.{station}.00.BDH",
        "channel_start_time": "2012-10-01T00:00:00.000000Z",
        "windows": 10,
        "pressure": pressure,
        "frequency_hz": f,
        "median_db": np.full(f.shape, float(level)),
    }


def test_off_scale_only_when_no_part_is_visible():
    assert off_scale([-3076.5] * 5, (-80, 100))
    assert not off_scale([-3076.5, 0.0], (-80, 100))
    assert off_scale([150.0, 160.0], (-80, 100))


def test_pressure_limits_use_reference_data_without_extrapolation():
    metadata, brown = load_references()
    figure, axis = plt.subplots()
    handles = draw_pressure_limits(axis, 0.02, 0.5)
    labels = [h.get_label() for h in handles]
    assert sum("Brown" in label for label in labels) == 2
    assert sum("Webb 1998" in label for label in labels) == 2
    brown_x = handles[0].get_xdata()
    # Brown starts at the 0.2 Hz split; Webb covers lower frequencies.
    assert brown_x.min() >= 0.2 and brown_x.max() <= 0.5
    inside = [p for p in metadata["webb"]["points"] if 0.02 <= p["frequency_hz"] <= 0.5]
    assert len(handles[2].get_xdata()) == len(inside)
    plt.close(figure)


def test_plot_flags_off_scale_channels_and_refuses_existing_output(tmp_path):
    curves = [curve("RR01", 10), curve("RR12", -3076.5), curve("RR52", -20)]
    ylim, flagged, stem = plot(curves, "YV", "BDH", None, tmp_path / "out")
    assert ylim == (-80.0, 100.0)
    assert [c["station"] for c in flagged] == ["RR12"]
    assert (tmp_path / "out" / f"{stem}.png").is_file()
    with pytest.raises(FileExistsError):
        plot(curves, "YV", "BDH", None, tmp_path / "out")
