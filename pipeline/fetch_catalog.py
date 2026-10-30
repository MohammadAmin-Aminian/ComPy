"""Stage 0: download the USGS earthquake catalogue used by stages B and D.

Writes config.CATALOG (QuakeML, M >= 5.85 over config.CATALOG_SPAN): stage B
removes earthquake spans of 1.5 days per magnitude unit above 5.85 (Deen et
al. 2017) and stage D uses the M >= 7 events for the DPG calibration.
"""

from obspy import UTCDateTime
from obspy.clients.fdsn import Client

from config import CATALOG, CATALOG_SPAN

MIN_MAGNITUDE = 5.85


def main():
    catalog = Client("USGS").get_events(
        starttime=UTCDateTime(CATALOG_SPAN[0]),
        endtime=UTCDateTime(CATALOG_SPAN[1]),
        minmagnitude=MIN_MAGNITUDE,
    )
    CATALOG.parent.mkdir(parents=True, exist_ok=True)
    catalog.write(str(CATALOG), format="QUAKEML")
    big = sum(1 for e in catalog if e.preferred_magnitude().mag >= 7)
    print(f"{len(catalog)} events (M >= 7: {big}) -> {CATALOG}")


if __name__ == "__main__":
    main()
