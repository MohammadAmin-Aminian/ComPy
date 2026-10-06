"""Download a chosen station interval, decimate and remove channel responses.

Example: python _Example/0_Download_data.py YV RR52 2012-12-01 2012-12-02 output.mseed
Review the decimation factors for your input sampling rate before running.
"""

import argparse
from pathlib import Path

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from tiskitpy import Decimator


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
    stream = client.get_waveforms(
        args.network, args.station, "*", "*", start, end, attach_response=True
    )
    stream.merge(method=0)
    stream = Decimator(args.factors).decimate(stream)
    for trace in stream:
        trace.remove_response(
            output="DEF" if trace.stats.channel.endswith("H") else "DISP"
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    stream.write(str(args.output), format="MSEED")


if __name__ == "__main__":
    main()
