"""Pressure-gauge gain from teleseismic Rayleigh waves, for one station.

The Rayleigh-wave calibration of compy_pressure, which replaces
Pressure_calibration.calculate_spectral_ratio (removed after ComPy 2.4.0):

    python _Example/3_DPG_Calibration.py SDS STATIONXML YV RR52 events.geojson

- SDS: an SDS archive with the station's pressure and seismometer channels
  at their native sample rate;
- STATIONXML: a StationXML file, or a folder of them, with full responses;
- events.geojson: a USGS ComCat catalogue covering the deployment, e.g.
  https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson&minmagnitude=6
  &starttime=2012-10-01&endtime=2013-12-01

Each Mw >= 6 teleseism gets a group-velocity window. The spectra of
pressure and vertical acceleration (responses removed) are compared with the
water-layer prediction in the bins whose squared coherence is at least 0.95,
after TiSKitPy tilt cleaning (``--tilt none``: pressure and vertical only).
The station gain is the median of at least three accepted events. It
multiplies the nominal pressure (calibrated = gain x nominal) and is never
applied to the data. Spectra, figures and comparison.json are written below
--output.
"""

import argparse
from pathlib import Path

import compy_pressure
from compy_pressure.network import apply_location_aliases, scan_files
from compy_pressure.rayleigh_network import calibration_master, run_network

DISPERSION = (
    Path(compy_pressure.__file__).parent
    / "reference_data/rr52_historical_rayleigh_dispersion.csv"
)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("data_dir", type=Path, help="SDS archive")
    parser.add_argument("inventory", type=Path, help="StationXML file or folder")
    parser.add_argument("network")
    parser.add_argument("station")
    parser.add_argument("catalog", type=Path, help="USGS ComCat GeoJSON")
    parser.add_argument("--band", type=float, nargs=2, default=[0.02, 0.06])
    parser.add_argument("--tilt", choices=("tiskitpy", "none"), default="tiskitpy")
    parser.add_argument("--pressure", default="BDH", help="channel code (--tilt none)")
    parser.add_argument("--vertical", default="BHZ", help="channel code (--tilt none)")
    parser.add_argument(
        "--dispersion-file",
        type=Path,
        default=DISPERSION,
        help="Rayleigh phase velocities (default: RR52, BRUIT-FM's default)",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or Path(f"calibration_{args.network}.{args.station}")
    channels = (args.pressure, args.vertical) if args.tilt == "none" else None
    master_dir, master = calibration_master(
        args.inventory, output / "master", args.network, args.station, channels=channels
    )
    files = scan_files(args.data_dir, [args.network], [args.station])
    files = apply_location_aliases(files, set(master.Seed_id))
    report = run_network(
        master_dir,
        output / "results",
        files,
        update_master=False,
        data_dir=args.data_dir,
        options=dict(
            catalog=str(args.catalog),
            dispersion_file=str(args.dispersion_file),
            band=list(args.band),
            tilt_method=args.tilt,
        ),
    )
    for entry in report["entries"].values():
        update = entry["update"]
        gain = update["rayleigh_pressure_multiplier"]
        print(
            f"{entry['seed_id']}: {entry['status']}, "
            f"{update['rayleigh_candidate_events']} accepted events, "
            f"gain {gain if gain != '' else 'not estimated'}\n"
            f"  details: {update['rayleigh_result_file']}"
        )


if __name__ == "__main__":
    main()
