"""Stage D: ComPy's teleseismic DPG gain for one station (gain 'compy').

Runs the Rayleigh-wave calibration of compy_pressure (from BRUIT-FM; the same code as ``compy-pressure rayleigh``):

- data: the station's native-rate SDS waveforms (config.SDS) and the network
  StationXML (config.INVENTORY);
- events: the saved USGS catalogue (converted to the ComCat GeoJSON fields it
  reads), Mw >= 6 teleseisms at most 150 km deep, with a magnitude schedule;
- method: a group-velocity window per event; pressure and vertical
  acceleration spectra with the full responses removed; bins with squared
  coherence >= 0.95 in CALIBRATION_BAND; a gain fitted per event against the
  water-layer prediction for the station depth;
- tilt cleaning: CALIBRATION_TILT ("tiskitpy" needs horizontals with valid
  orientations; "none" uses pressure and vertical only);
- gain: the median of at least three accepted events.

Nothing is applied to the data.

Outputs, read by stages C, F and G:

- work/<STA>/calibration/compy_gain.json: gain, gain_std, event_gains,
  events_used, band_hz, depth_m, or error;
- <STA>_compy_calibration.npz: frequencies, coherences, measured_ratio,
  theoretical_ratio, band, gain, event_gains, depth;
- <STA>_compy_calibration.png.

The full run (per-event spectra and figures, comparison.json) is in
work/<STA>/calibration/compy_pressure/.
"""

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

import config
import numpy as np
from config import CATALOG, COMPY, INVENTORY, NETWORK, SDS, WORK, P, Z
from obspy import read_events

sys.path.insert(0, str(COMPY))

from compy_pressure import network as pressure_network  # noqa: E402
from compy_pressure.rayleigh_network import (  # noqa: E402
    calibration_master,
    run_network,
)

BAND = tuple(getattr(config, "CALIBRATION_BAND", (0.02, 0.06)))
TILT = getattr(config, "CALIBRATION_TILT", "tiskitpy")
DISPERSION = Path(
    getattr(
        config,
        "CALIBRATION_DISPERSION",
        COMPY / "compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv",
    )
)
DENSITY = 1028.0  # compy_pressure's default water density, kg/m^3


def geojson_catalog(quakeml, destination):
    """Write the ComCat GeoJSON fields compy_pressure reads from a QuakeML file."""
    features = []
    for event in read_events(str(quakeml)):
        origin = event.preferred_origin() or event.origins[0]
        magnitude = event.preferred_magnitude() or event.magnitudes[0]
        resource = str(event.resource_id)
        # USGS ids end in "eventid=us1000abcd&format=quakeml".
        match = re.search(r"eventid=([A-Za-z0-9_-]+)", resource)
        name = match[1] if match else re.sub(r"[^A-Za-z0-9_-]", "_", resource)[-40:]
        features.append(
            {
                "id": name,
                "properties": {
                    "time": origin.time.timestamp * 1000,
                    "mag": magnitude.mag,
                    "magType": magnitude.magnitude_type or "",
                    "url": resource,
                },
                "geometry": {
                    "coordinates": [
                        origin.longitude,
                        origin.latitude,
                        (origin.depth or 0) / 1000,
                    ]
                },
            }
        )
    destination.write_text(
        json.dumps({"type": "FeatureCollection", "features": features})
    )
    return destination


def figure_arrays(summary, depth):
    """The arrays of stage F's calibration figure, on one frequency grid."""
    events = [e for e in summary["events"] if e["status"] == "candidate"]
    spectra = [np.load(e["files"]["spectra"]) for e in events]
    f = spectra[0]["frequency_hz"]

    def on_grid(s, key):
        if len(s["frequency_hz"]) == len(f) and np.allclose(s["frequency_hz"], f):
            return s[key]
        return np.interp(f, s["frequency_hz"], s[key])

    settings = summary["settings"]
    band = np.array(settings.get("band_hz") or settings["band"], float)
    inside = (f >= band[0] - 1e-9) & (f <= band[1] + 1e-9)
    scale = DENSITY * depth
    return {
        "frequencies": f,
        "coherences": np.array([on_grid(s, "coherence_squared") for s in spectra]),
        "measured_ratio": np.median(
            [on_grid(s, "psd_ratio_pa") for s in spectra], axis=0
        )
        / scale,
        "theoretical_ratio": on_grid(spectra[0], "expected_pa")[inside] / scale,
        "band": band,
        "gain": summary["median_pressure_multiplier"],
        "event_gains": np.array([e["pressure_multiplier"] for e in events]),
        "depth": depth,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("station")
    args = parser.parse_args()
    station = args.station
    out = WORK / station / "calibration"
    run = out / "compy_pressure"
    run.mkdir(parents=True, exist_ok=True)
    result = {
        "station": station,
        "method": "compy_pressure Rayleigh-wave calibration (BRUIT-FM)",
        "band_hz": list(BAND),
        "tilt_method": TILT,
        "dispersion_file": str(DISPERSION),
    }
    try:
        # The pipeline makes no PPSDs: the StationXML is each channel's response
        # snapshot. Without tilt cleaning only the pressure and vertical are
        # kept, which also avoids pairing errors from wrong horizontal dips
        # (XE: BH1 and BH2 at -90 degrees look like verticals).
        master_dir, master = calibration_master(
            INVENTORY,
            run / "master",
            NETWORK,
            station,
            channels=(P, Z) if TILT == "none" else None,
        )
        catalog = geojson_catalog(CATALOG, run / "catalog.geojson")
        files = pressure_network.scan_files(SDS, [NETWORK], [station])
        files = pressure_network.apply_location_aliases(files, set(master.Seed_id))
        report = run_network(
            master_dir,
            run / "results",
            files,
            update_master=False,
            data_dir=SDS,
            options=dict(
                catalog=str(catalog),
                dispersion_file=str(DISPERSION),
                band=list(BAND),
                tilt_method=TILT,
            ),
        )
        # One gain per station: the pressure epoch with most accepted events.
        entries = [
            e
            for e in report["entries"].values()
            if e["seed_id"].startswith(f"{NETWORK}.{station}.")
        ]
        if not entries:
            raise ValueError("No pressure channel with SDS waveforms for this station")
        best = max(entries, key=lambda e: e["update"]["rayleigh_candidate_events"])
        summary = json.loads(Path(best["update"]["rayleigh_result_file"]).read_text())
        row = master[master.Seed_id.eq(best["seed_id"])].iloc[0]
        depth = float(row.Profondeur)
        accepted = [e for e in summary["events"] if e["status"] == "candidate"]
        result.update(
            seed_id=best["seed_id"],
            status=summary["status"],
            events_considered=len(summary["events"]),
            candidate_events=len(accepted),
            depth_m=depth,
            result_file=best["update"]["rayleigh_result_file"],
        )
        if summary.get("median_pressure_multiplier") is None:
            raise ValueError(summary.get("reason") or summary["status"])
        result.update(
            gain=float(summary["median_pressure_multiplier"]),
            gain_std=summary.get("event_multiplier_std"),
            event_gains=[float(e["pressure_multiplier"]) for e in accepted],
            events_used=[e["origin_utc"] for e in accepted],
        )
        np.savez(
            out / f"{station}_compy_calibration.npz", **figure_arrays(summary, depth)
        )
        figure = summary["files"].get("paper_comparison")
        if figure:
            shutil.copy2(figure, out / f"{station}_compy_calibration.png")
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        result["error"] = str(exc)
    (out / "compy_gain.json").write_text(json.dumps(result, indent=1))
    print(
        json.dumps({k: v for k, v in result.items() if k != "event_gains"}), flush=True
    )


if __name__ == "__main__":
    main()
