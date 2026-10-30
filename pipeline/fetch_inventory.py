"""Stage 0: StationXML with the seismometer and pressure responses.

Requests level=response metadata for config.NETWORK, config.STATIONS and the
channels Z, H1, H2 and P over config.YEARS from an FDSN service (default the
IRIS federator) and writes config.INVENTORY. Skip this stage if the operator
supplied a complete StationXML.

Usage: python fetch_inventory.py [--client IRIS]
"""

import argparse

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.routing.routing_client import RoutingClient

from config import CHANNELS, INVENTORY, NETWORK, STATIONS, YEARS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--client", default="iris-federator",
                        help="FDSN client name, or iris-federator / eida-routing")  # fmt: skip
    args = parser.parse_args()
    client = RoutingClient(args.client) if "routing" in args.client or "federator" in args.client \
        else Client(args.client)  # fmt: skip
    inventory = client.get_stations(
        network=NETWORK,
        station=",".join(STATIONS),
        channel=",".join(CHANNELS),
        starttime=UTCDateTime(YEARS[0], 1, 1),
        endtime=UTCDateTime(YEARS[-1] + 1, 1, 1),
        level="response",
    )
    INVENTORY.parent.mkdir(parents=True, exist_ok=True)
    inventory.write(str(INVENTORY), format="STATIONXML")
    n = sum(len(s.channels) for net in inventory for s in net)
    print(f"{n} channels -> {INVENTORY}")


if __name__ == "__main__":
    main()
