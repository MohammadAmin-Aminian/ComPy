"""The pressure calibration in ComPy: calibration_master and pipeline stage D."""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from obspy import UTCDateTime, read_inventory
from obspy.core.event import Catalog, Event, Magnitude, Origin, ResourceIdentifier
from obspy.core.inventory import Channel, Inventory, Network, Site, Station
from obspy.core.inventory.response import Response

from compy_pressure.rayleigh_calibration import catalog_events, run_calibration
from compy_pressure.rayleigh_network import calibration_master

PIPELINE = Path(__file__).resolve().parents[1] / "pipeline"


def write_stationxml(path, dip_horizontals=0.0):
    """One OBS: pressure, vertical and two horizontals (optionally mislabelled)."""
    start = UTCDateTime("2020-01-01")

    def channel(code, units, dip, azimuth):
        response = Response.from_paz(
            zeros=[], poles=[], stage_gain=4, input_units=units, output_units="COUNTS"
        )
        response.instrument_sensitivity.frequency = 0.1
        return Channel(
            code=code,
            location_code="",
            latitude=1,
            longitude=2,
            elevation=-4000,
            depth=0,
            dip=dip,
            azimuth=azimuth,
            sample_rate=1.0,
            start_date=start,
            end_date=start + 86400,
            response=response,
        )

    channels = [
        channel("BDH", "PA", -90, 0),
        channel("BHZ", "M/S", -90, 0),
        channel("BH1", "M/S", dip_horizontals, 0),
        channel("BH2", "M/S", dip_horizontals, 90),
    ]
    station = Station(
        "OBS1", 1, 2, -4000, site=Site(name="Synthetic"), channels=channels
    )
    Inventory([Network("XX", stations=[station])], source="Test").write(
        str(path), format="STATIONXML"
    )
    return path


def test_calibration_master_uses_each_stationxml_as_response_snapshot(tmp_path):
    xml = write_stationxml(tmp_path / "XX.xml")
    master_dir, master = calibration_master(xml, tmp_path / "cal", "XX", "OBS1")
    rows = master[master.Seed_id.str.startswith("XX.OBS1.")]
    assert sorted(rows.Seed_id.str[-3:]) == ["BDH", "BH1", "BH2", "BHZ"]
    for snapshot in rows.rayleigh_response_snapshot_file:
        assert Path(snapshot).read_bytes() == xml.read_bytes()
    # A rerun reuses the saved table; a changed StationXML is refused.
    again = calibration_master(xml, master_dir, "XX", "OBS1")[1]
    assert again.Seed_id.tolist() == master.Seed_id.tolist()
    write_stationxml(xml, dip_horizontals=-90)
    with pytest.raises(ValueError, match="changed"):
        calibration_master(xml, master_dir, "XX", "OBS1")
    with pytest.raises(ValueError, match="not in"):
        calibration_master(xml, tmp_path / "other", "XX", "NONE")


def test_calibration_master_keeps_only_the_named_channels(tmp_path):
    """XE-like metadata: horizontals at dip -90 would pair as verticals."""
    xml = write_stationxml(tmp_path / "XX.xml", dip_horizontals=-90)
    _, master = calibration_master(
        xml, tmp_path / "cal", "XX", "OBS1", channels=("BDH", "BHZ")
    )
    assert sorted(master.Seed_id) == ["XX.OBS1..BDH", "XX.OBS1..BHZ"]


@pytest.mark.parametrize(
    "station,channels", [("OTHER", ("BDH", "BHZ")), ("OBS1", None)]
)
def test_calibration_master_does_not_reuse_another_selection(
    tmp_path, station, channels
):
    xml = write_stationxml(tmp_path / "XX.xml")
    folder, _ = calibration_master(
        xml, tmp_path / "cal", "XX", "OBS1", channels=("BDH", "BHZ")
    )
    with pytest.raises(ValueError, match="inputs changed"):
        calibration_master(xml, folder, "XX", station, channels=channels)


def test_calibration_master_detects_a_removed_inventory_source(tmp_path):
    sources = tmp_path / "xml"
    sources.mkdir()
    write_stationxml(sources / "XX.xml")
    extra = write_stationxml(sources / "extra.xml")
    folder, _ = calibration_master(sources, tmp_path / "cal", "XX", "OBS1")
    extra.unlink()
    with pytest.raises(ValueError, match="inputs changed"):
        calibration_master(sources, folder, "XX", "OBS1")


def test_tilt_training_reads_the_explicit_response_snapshot(tmp_path):
    """Regression: the TiSKitPy path read <PPSD folder>/inventory.xml.

    Without PPSDs (calibration_master), tilt cleaning raised AttributeError.
    """
    from compy_pressure.rayleigh_calibration import physical_quiet_stream

    xml = write_stationxml(tmp_path / "XX.xml")
    _, master = calibration_master(xml, tmp_path / "cal", "XX", "OBS1")
    assert "ppsd_npz_file" not in master or not master.ppsd_npz_file.any()
    inventory = read_inventory(str(xml))
    rows, channels = [], []
    for code in ("BDH", "BHZ", "BH1", "BH2"):
        rows.append(master[master.Seed_id.eq(f"XX.OBS1..{code}")].iloc[0])
        channels.append(inventory.select(channel=code)[0][0][0])
    arrays = np.random.default_rng(1).normal(size=(4, 4096))
    stream = physical_quiet_stream(
        arrays, rows, channels, 1.0, UTCDateTime("2020-01-01T01:00:00")
    )
    # The tilt-training stream holds the seismometer channels only.
    assert [tr.stats.channel for tr in stream] == ["BHZ", "BH1", "BH2"]


def test_invalid_band_is_refused_before_any_data_is_read(tmp_path):
    for band in ((0.06, 0.02), (0.0, 0.06)):
        with pytest.raises(ValueError, match="Invalid band"):
            run_calibration(catalog=tmp_path / "none.geojson", band=band)


@pytest.fixture
def stage_d(monkeypatch, tmp_path):
    """pipeline/d_calibrate.py, imported with its config pointed at tmp_path."""
    monkeypatch.setenv("COMPY_RUN_DIR", str(tmp_path))
    monkeypatch.syspath_prepend(str(PIPELINE))
    monkeypatch.delitem(sys.modules, "config", raising=False)
    spec = importlib.util.spec_from_file_location(
        "d_calibrate", PIPELINE / "d_calibrate.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    yield module
    sys.modules.pop("config", None)


def test_stage_d_converts_the_usgs_quakeml_catalogue(stage_d, tmp_path):
    event = Event(
        resource_id=ResourceIdentifier(
            "quakeml:earthquake.usgs.gov/fdsnws/event/1/query?eventid=us1000abcd&format=quakeml"
        ),
        origins=[
            Origin(
                time=UTCDateTime("2018-08-19T00:19:40.67"),
                latitude=-18.1,
                longitude=-178.1,
                depth=600000,
            )
        ],
        magnitudes=[Magnitude(mag=8.2, magnitude_type="mww")],
    )
    quakeml = tmp_path / "catalog.qml"
    Catalog([event]).write(str(quakeml), format="QUAKEML")
    events = catalog_events(stage_d.geojson_catalog(quakeml, tmp_path / "c.geojson"))
    assert events[0]["id"] == "us1000abcd"
    assert UTCDateTime(events[0]["origin_utc"]) == UTCDateTime("2018-08-19T00:19:40.67")
    assert events[0]["depth_km"] == pytest.approx(600)
    assert (events[0]["magnitude"], events[0]["magnitude_type"]) == (8.2, "mww")


def test_stage_d_figure_arrays_match_stage_f(stage_d, tmp_path):
    f = np.linspace(0.01, 0.1, 10)
    events = []
    for i, gain in enumerate((0.69, 0.71, 0.70)):
        path = tmp_path / f"event{i}_spectra.npz"
        np.savez(
            path,
            frequency_hz=f,
            coherence_squared=np.full(len(f), 0.98),
            psd_ratio_pa=np.full(len(f), 2.0e7),
            expected_pa=np.full(len(f), 1.4e7),
        )
        events.append(
            {
                "status": "candidate",
                "pressure_multiplier": gain,
                "files": {"spectra": str(path)},
            }
        )
    events.append({"status": "unresolved"})
    summary = {
        "events": events,
        "settings": {"band_hz": [0.02, 0.06]},
        "median_pressure_multiplier": 0.70,
    }
    arrays = stage_d.figure_arrays(summary, depth=4000.0)
    inside = (f >= 0.02 - 1e-9) & (f <= 0.06 + 1e-9)  # stage F's band mask
    assert arrays["coherences"].shape == (3, len(f))
    assert arrays["measured_ratio"].shape == f.shape
    assert arrays["theoretical_ratio"].shape == (inside.sum(),)
    assert arrays["measured_ratio"][0] == pytest.approx(2.0e7 / (1028 * 4000))
    assert arrays["event_gains"].tolist() == [0.69, 0.71, 0.70]
    assert (arrays["gain"], arrays["depth"]) == (0.70, 4000.0)
    json.dumps({k: np.asarray(v).tolist() for k, v in arrays.items()})
