"""Download a chosen station interval, decimate and remove channel responses.

Example: python _Example/0_Download_data.py YV RR52 2012-12-01 2012-12-02 output.mseed
Review the decimation factors for your input sampling rate before running.
"""

import argparse
from copy import deepcopy
from pathlib import Path

from obspy import UTCDateTime
from obspy.core.inventory import Inventory
from obspy.clients.fdsn import Client
from tiskitpy import Decimator
from tiskitpy.decimate import FIRFilter


def decimate_with_responses(stream, inventory, factors):
    """Keep channel IDs and added FIR response stages consistent with decimation."""
    decimator = Decimator(factors)
    decimated = decimator.decimate(stream, keep_dtype=False)
    updated = Inventory([], source=inventory.source)
    for original, output in zip(stream, decimated):
        metadata = inventory.select(
            network=original.stats.network,
            station=original.stats.station,
            location=original.stats.location,
            channel=original.stats.channel,
            time=original.stats.starttime,
        ).copy()
        channels = [
            channel
            for network in metadata
            for station in network
            for channel in station
        ]
        if len(channels) != 1:
            raise ValueError(f"Expected one response epoch for {original.id}")
        channel = channels[0]
        if channel.end_date is not None and channel.end_date < original.stats.endtime:
            raise ValueError(
                f"Response epoch ends inside {original.id}; split at the epoch boundary"
            )
        if channel.sample_rate != original.stats.sampling_rate:
            raise ValueError("Waveform and inventory sampling rates differ")
        response = channel.response
        if response is None or not response.response_stages:
            raise ValueError(f"Missing response stages for {original.id}")
        rate = original.stats.sampling_rate
        sequence = response.response_stages[-1].stage_sequence_number
        for factor in decimator.decimates:
            sequence += 1
            response.response_stages.append(
                FIRFilter.from_SAC(factor).to_obspy(
                    rate, sequence, response.response_stages[-1].output_units
                )
            )
            rate /= factor
        # Evalresp does not recognize pressure units when recalculating gain.
        # The numerical sensitivity is independent of this temporary unit label.
        units = response.instrument_sensitivity.input_units
        if units.upper() == "PA":
            response.instrument_sensitivity.input_units = "M/S"
        try:
            response.recalculate_overall_sensitivity()
        finally:
            response.instrument_sensitivity.input_units = units
        channel.code = output.stats.channel
        channel.location_code = output.stats.location
        channel.sample_rate = rate
        output.stats.response = deepcopy(response)
        updated += metadata
    return decimated, updated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("network")
    parser.add_argument("station")
    parser.add_argument("start")
    parser.add_argument("end")
    parser.add_argument("output", type=Path)
    parser.add_argument("--server", default="RESIF")
    parser.add_argument("--factors", type=int, nargs="+", default=[5, 5, 2])
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    start, end = UTCDateTime(args.start), UTCDateTime(args.end)
    if end <= start:
        raise ValueError("end must follow start")
    client = Client(args.server)
    inventory = client.get_stations(
        network=args.network,
        station=args.station,
        location="*",
        channel="*",
        starttime=start,
        endtime=end,
        level="response",
    )
    stream = client.get_waveforms(args.network, args.station, "*", "*", start, end)
    stream.merge(method=0)
    stream, inventory = decimate_with_responses(stream, inventory, args.factors)
    for trace in stream:
        trace.remove_response(
            output="DEF" if trace.stats.channel.endswith("H") else "DISP"
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    stream.write(str(args.output), format="MSEED")


if __name__ == "__main__":
    main()
