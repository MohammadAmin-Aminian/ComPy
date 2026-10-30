"""Pressure-gauge QC with compy_pressure: master table, PPSDs, plots and QC tests.

Without arguments, this writes a small synthetic SDS archive and StationXML
(one pressure gauge and one vertical, eight hours of white noise) to --output
and runs every step offline in a few seconds. With --data-dir, --inventory-dir and --network
it runs the same steps on real data:

    python _Example/5_Pressure_QC.py --output qc_run \
        --data-dir /path/to/SDS --inventory-dir /path/to/stationxml --network YV

The command-line equivalent is ``compy-pressure inventory|ppsd|plots|pressure|
nearby|infragravity ...``; see docs/PRESSURE.md. Each step resumes from its
checkpoints when it is rerun. Nothing changes a response or the data.
"""

import argparse
from pathlib import Path

import numpy as np
from obspy import Trace, UTCDateTime
from obspy.core.inventory import Channel, Inventory, Network, Site, Station
from obspy.core.inventory.response import Response

from compy_pressure import network
from compy_pressure.inventory import load_master_dataset, synchronize_inventory
from compy_pressure.ppsd import is_pressure_unit

# QC test 1 records the fraction of median PSD bins inside the Webb/Brown
# seafloor limits; QC test 3 the median residual from the Godin et al. (2013)
# infragravity model at the station depth. Scores are given by a later policy.
COLUMNS = [
    "Seed_id",
    "ppsd_status",
    "pressure_expected_tested_bins",
    "pressure_expected_median_in_range_fraction",
    "infragravity_status",
    "infragravity_compared_bins",
    "infragravity_median_residual_db",
]


def write_synthetic(root, seconds=8 * 3600, rate=1.0):
    """One gauge (Pa) and one vertical (m/s) at 1000 m with flat responses."""
    start = UTCDateTime("2020-01-01")
    channels = []
    for code, units in (("BDH", "PA"), ("BHZ", "M/S")):
        response = Response.from_paz(
            zeros=[], poles=[], stage_gain=4, input_units=units, output_units="COUNTS"
        )
        response.instrument_sensitivity.frequency = 0.1  # below Nyquist
        channels.append(
            Channel(
                code=code,
                location_code="00",
                latitude=1,
                longitude=2,
                elevation=-1000,
                depth=0,
                sample_rate=rate,
                start_date=start,
                end_date=start + seconds,
                response=response,
            )
        )
        trace = Trace(
            np.random.default_rng(9).integers(-1000, 1000, int(seconds * rate))
        )
        trace.data = trace.data.astype(np.int32)
        trace.stats.update(
            dict(
                network="XX",
                station="SYN",
                location="00",
                channel=code,
                sampling_rate=rate,
                starttime=start,
            )
        )
        folder = root / "SDS" / "2020" / "XX" / "SYN" / f"{code}.D"
        folder.mkdir(parents=True, exist_ok=True)
        trace.write(str(folder / f"XX.SYN.00.{code}.D.2020.001"), format="MSEED")
    station = Station(
        code="SYN",
        latitude=1,
        longitude=2,
        elevation=-1000,
        site=Site(name="Synthetic"),
        channels=channels,
    )
    (root / "stationxml").mkdir(parents=True, exist_ok=True)
    Inventory([Network(code="XX", stations=[station])], source="Example").write(
        str(root / "stationxml" / "XX.xml"), format="STATIONXML"
    )
    return root / "SDS", root / "stationxml", ["XX"]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=Path("pressure_qc_example"))
    parser.add_argument("--data-dir", type=Path, help="SDS archive")
    parser.add_argument("--inventory-dir", type=Path, help="StationXML folder")
    parser.add_argument("--network", action="append")
    parser.add_argument("--workers", type=network.worker_count, default=2)
    args = parser.parse_args()
    output = args.output.resolve()
    if args.data_dir:
        data, inventories, networks = args.data_dir, args.inventory_dir, args.network
        # ComPy's defaults: 3 h windows with 50 % overlap, down to 1 mHz.
        ppsd = dict(ppsd_length=10800, minimum_frequency=0.001, chunk_seconds=86400)
        ppsd["overlap"] = 0.5
        stages = ["ppsd", "plots", "pressure", "nearby", "infragravity"]
    else:
        data, inventories, networks = write_synthetic(output / "inputs")
        # Two-hour windows fit the eight-hour record and reach 2 mHz.
        ppsd = dict(ppsd_length=7200, minimum_frequency=0.002, chunk_seconds=86400)
        ppsd["overlap"] = 0.5
        stages = ["ppsd", "plots", "pressure", "infragravity"]
    master, results = output / "master", output / "results"

    print(synchronize_inventory(inventories, master))
    known = set(load_master_dataset(master).Seed_id)
    files = network.apply_location_aliases(network.scan_files(data, networks), known)
    for stage in stages:
        network.run_stage(
            stage,
            master,
            results,
            files,
            workers=args.workers,
            options=ppsd if stage == "ppsd" else {},
            pressure_only=True,
        )

    table = load_master_dataset(master)
    pressure = table[table.input_units.map(is_pressure_unit)]
    print(pressure[[c for c in COLUMNS if c in pressure]].T.to_string(header=False))
    print(f"Master table: {master / 'master_inventory.csv'}; results: {results}")


if __name__ == "__main__":
    main()
