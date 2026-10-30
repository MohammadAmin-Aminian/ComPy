"""Stage D: ComPy's own teleseismic DPG gain for one station (gain 'compy').

Runs Pressure_calibration.calculate_spectral_ratio (ComPy 2.2, defaults:
M >= 7 events, 0.02–0.06 Hz, magnitude-squared coherence median and mean
> 0.97, 512-sample Welch segments) on the decimated raw counts (1.25 Hz)
with the YV StationXML responses. Event start times come from the saved
USGS catalogue; events whose 2-hour window contains a data gap are skipped
because gaps are zero-filled for the calibration call.

Output: work/<STA>/calibration/compy_gain.json and a diagnostic figure.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from obspy import Stream, read, read_inventory

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "ComPy"))

import Pressure_calibration as dpg  # noqa: E402
import tiskitpy as tiskit  # noqa: E402

WORK = ROOT / "work"
CATALOG = WORK / "catalog_usgs_M5.85_20121001_20131201.qml"
INVENTORY = "/Volumes/DATA_bruitfm/INVENTORIES/YV.xml"
MIN_MAG = 7.0
EVENT_WINDOW_S = 7200


def main():
    import logging

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("station")
    args = parser.parse_args()
    logging.disable(logging.INFO)
    station = args.station
    out = WORK / station / "calibration"
    out.mkdir(parents=True, exist_ok=True)
    st = Stream()
    for f in sorted((WORK / station / "counts").glob(f"{station}.*.mseed")):
        day = read(str(f))
        st += day.select(channel="BHZ") + day.select(channel="BDH")
    st.merge(method=0)
    start = max(tr.stats.starttime for tr in st)
    end = min(tr.stats.endtime for tr in st)
    spans = tiskit.TimeSpans.from_eqs(
        (start, end), minmag=MIN_MAG, days_per_magnitude=0.5,
        eq_file=str(CATALOG), save_eq_file=False,
    )  # fmt: skip
    usable = []
    for t0 in spans.start_times:
        t1 = t0 + EVENT_WINDOW_S
        if t0 < start or t1 > end:
            continue
        per_channel = []
        for channel in ("BHZ", "BDH"):
            parts = [(tr.stats.starttime, tr.stats.endtime)
                     for tr in st.select(channel=channel).split()]  # fmt: skip
            per_channel.append(any(a <= t0 and b >= t1 for a, b in parts))
        if all(per_channel):
            usable.append(t0)
    st.merge(method=1, fill_value=0)
    st.trim(start, end)
    inv = read_inventory(INVENTORY).select(station=station)
    event_spans = tiskit.TimeSpans(
        start_times=usable, end_times=[t + EVENT_WINDOW_S for t in usable]
    )
    result = {"station": station, "events_considered": len(usable),
              "min_magnitude": MIN_MAG}  # fmt: skip
    try:
        details = dpg.calculate_spectral_ratio(
            st, mag=MIN_MAG, inventory=inv, event_spans=event_spans,
            plot_condition=True, return_details=True,
        )  # fmt: skip
        import matplotlib.pyplot as plt

        plt.gcf().suptitle(
            f"YV.{station} ComPy DPG calibration, g = {details['gain']:.3f}"
        )
        plt.savefig(out / f"{station}_compy_calibration.png", dpi=110)
        plt.close("all")
        result.update({
            "gain": float(details["gain"]),
            "event_gains": details["event_gains"].tolist(),
            "gain_std": details["gain_std"],
            "events_used": [str(t) for t in details["event_starts"]],
            "band_hz": details["band"],
            "depth_m": details["depth"],
        })  # fmt: skip
    except ValueError as exc:
        result["error"] = str(exc)
    (out / "compy_gain.json").write_text(json.dumps(result, indent=1))
    print(
        json.dumps({k: v for k, v in result.items() if k != "event_gains"}), flush=True
    )


if __name__ == "__main__":
    main()
