"""Constant-data PSD windows are flagged and left out of plot scaling."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import Trace, UTCDateTime
from obspy.core.inventory import Channel, Inventory, Network, Site, Station
from obspy.core.inventory.response import Response
from obspy.signal import PPSD

from compy_pressure import plot_ppsd
from compy_pressure.screening import dropout_windows, screened_values


def pressure_ppsd():
    response = Response.from_paz(
        zeros=[], poles=[], stage_gain=1e4, input_units="PA", output_units="COUNTS"
    )
    response.instrument_sensitivity.frequency = 1
    start = UTCDateTime("2020-01-01")
    channel = Channel(
        code="BDH",
        location_code="",
        latitude=0,
        longitude=0,
        elevation=-3000,
        depth=0,
        sample_rate=8,
        start_date=start,
        response=response,
    )
    station = Station(
        code="T",
        latitude=0,
        longitude=0,
        elevation=-3000,
        site=Site(name="T"),
        channels=[channel],
    )
    inventory = Inventory([Network(code="XX", stations=[station])])
    data = np.random.default_rng(1).normal(0, 1000, 8 * 4096).astype(np.int32)
    data[8 * 2048 :] = 0  # second half dead: constant data
    trace = Trace(data)
    trace.stats.update(dict(network="XX", station="T", channel="BDH", sampling_rate=8))
    trace.stats.starttime = start
    ppsd = PPSD(
        trace.stats,
        inventory,
        ppsd_length=512,
        overlap=0.5,
        special_handling="hydrophone",
        db_bins=(-250, 100, 1.0),
        period_limits=(0.25, 50),
    )
    ppsd.add(trace)
    return ppsd


def test_constant_windows_are_flagged_and_screened():
    ppsd = pressure_ppsd()
    dropouts = dropout_windows(ppsd)
    values = np.asarray(ppsd.psd_values)
    assert dropouts.any() and not dropouts.all()
    assert (values[dropouts] < -250).all()
    assert (values[~dropouts] > -250).all()
    assert len(screened_values(ppsd)) == (~dropouts).sum()


def test_plots_scale_to_real_windows(tmp_path):
    ppsd = pressure_ppsd()
    values = np.asarray(ppsd.psd_values, dtype=float)
    values[dropout_windows(ppsd)] = np.nan
    for factory in (plot_ppsd.plot_histogram, plot_ppsd.plot_all_lines):
        figure = factory(ppsd, values)
        low, high = figure.axes[0].get_ylim()
        assert low > -300 and high <= 120
        plt.close(figure)
    figure = plot_ppsd.plot_all_lines(ppsd, values)
    assert "constant-data windows left out" in figure.axes[0].get_title()
    plt.close(figure)
    plt.close(plot_ppsd.plot_spectrogram(ppsd, values))
