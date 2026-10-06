import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import Stream, Trace, UTCDateTime
from obspy.core.inventory import Channel, Inventory, Network, Station
from obspy.core.inventory.response import (
    InstrumentSensitivity,
    PolesZerosResponseStage,
    Response,
)

from Pressure_calibration import _date_ticks, plot_spectrogram


def test_short_spectrogram_does_not_divide_by_zero():
    rng = np.random.default_rng(4)
    stream = Stream(
        [
            Trace(rng.normal(size=20000), header={"channel": c, "sampling_rate": 2})
            for c in ["BHZ", "BDH"]
        ]
    )
    plot_spectrogram(stream)
    plt.close("all")


def test_single_bin_date_ticks():
    positions, dates = _date_ticks(UTCDateTime(0), UTCDateTime(60), 1)
    assert positions.tolist() == [0]
    assert len(dates) == 1


def test_decimated_response_includes_fir_filter():
    path = Path(__file__).parents[1] / "_Example/0_Download_data.py"
    spec = importlib.util.spec_from_file_location("download_example", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    response = Response(
        instrument_sensitivity=InstrumentSensitivity(1, 1, "M/S", "COUNTS"),
        response_stages=[
            PolesZerosResponseStage(
                stage_sequence_number=1,
                stage_gain=1,
                stage_gain_frequency=1,
                input_units="M/S",
                output_units="COUNTS",
                pz_transfer_function_type="LAPLACE (RADIANS/SECOND)",
                normalization_frequency=1,
                zeros=[0j],
                poles=[-1 + 0j],
            )
        ],
    )
    channel = Channel("BHZ", "", 0, 0, 0, 0, sample_rate=20, response=response)
    inventory = Inventory(
        [Network("XX", stations=[Station("TEST", 0, 0, 0, channels=[channel])])], "test"
    )
    original = Stream(
        [
            Trace(
                np.sin(np.arange(1000)),
                header={
                    "network": "XX",
                    "station": "TEST",
                    "channel": "BHZ",
                    "sampling_rate": 20,
                },
            )
        ]
    )
    decimated, updated = module.decimate_with_responses(original, inventory, [2])
    assert decimated[0].stats.sampling_rate == 10
    assert len(decimated[0].stats.response.response_stages) == 2
    assert decimated[0].stats.response == updated.get_response(
        decimated[0].id, decimated[0].stats.starttime
    )
    assert len(inventory[0][0][0].response.response_stages) == 1
    assert original[0].stats.sampling_rate == 20


def test_display_smoothing_handles_short_arrays():
    from compy_numerics import plot_smooth

    for count in [1, 2, 3, 4, 7, 20]:
        data = np.arange(count, dtype=float)
        np.testing.assert_allclose(plot_smooth(data, 31, 3), data, atol=1e-10)
