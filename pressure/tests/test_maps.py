import pytest

pytest.importorskip("cartopy")

import matplotlib.pyplot as plt  # noqa: E402

from compy_pressure.maps import PLATE, add_land  # noqa: E402


def test_add_land_sets_the_requested_extent():
    # Outlines are only fetched when drawn, so this needs no network access.
    figure, axis = plt.subplots(subplot_kw={"projection": PLATE})
    add_land(axis, (40, 80, -40, 0), labels=False)
    assert axis.get_extent(crs=PLATE) == pytest.approx((40, 80, -40, 0))
    plt.close(figure)
