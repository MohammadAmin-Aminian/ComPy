"""Defects found by the audit of 20 and 21 November 2026 (ComPy 3.2.1).

Each test fails on 3.2.0. Stage A runs here on a small archive sampled at
2.5 Hz (one FIR stage to 1.25 Hz), so that whole days stay cheap.
"""

import importlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest
from obspy import Trace, UTCDateTime, read
from obspy.core.inventory import Channel, Inventory, Network, Site, Station
from obspy.core.inventory.response import Response

REPO = Path(__file__).resolve().parents[1]
T0 = UTCDateTime(2020, 1, 1)
EPOCH = T0 + 1.5 * 86400
RATE = 2.5
CHANNELS = ("BHZ", "BH1", "BH2", "BDH")
CONFIG = (
    'NETWORK = "XX"\nSTATIONS = ("OBS",)\nYEARS = (2020,)\nRAW_RATE = 2.5\n'
    'FACTORS = (2,)\nINVENTORY = ROOT / "metadata" / "XX.xml"\n'
)


def _channel(code, units, gain, start, end):
    response = Response.from_paz(
        zeros=[], poles=[], stage_gain=gain, input_units=units, output_units="COUNTS"
    )
    response.instrument_sensitivity.frequency = 0.05
    return Channel(code=code, location_code="00", latitude=1, longitude=2, elevation=-3000,
                   depth=0, sample_rate=RATE, start_date=start, end_date=end,
                   response=response)  # fmt: skip


def _write_day(root, day, rate=RATE, channels=CHANNELS):
    start = T0 + day * 86400
    n = int(86400 * rate)
    t = start.timestamp + np.arange(n) / rate
    noise = np.random.default_rng(day).standard_normal(n)
    data = (1000 * np.sin(2 * np.pi * 0.05 * t) + 20 * noise).astype(np.int32)
    for code in channels:
        stats = dict(network="XX", station="OBS", location="00", channel=code,
                     sampling_rate=rate, starttime=start)  # fmt: skip
        folder = root / "SDS" / "2020" / "XX" / "OBS" / f"{code}.D"
        folder.mkdir(parents=True, exist_ok=True)
        name = f"XX.OBS.00.{code}.D.2020.{start.julday:03d}"
        Trace(data, stats).write(str(folder / name), format="MSEED")


@pytest.fixture
def run_dir(tmp_path, monkeypatch):
    """A run folder with a StationXML whose pressure gain doubles at EPOCH."""
    channels = [
        _channel(c, "M/S", 4.0, T0 - 86400, T0 + 30 * 86400) for c in CHANNELS[:3]
    ]
    channels += [_channel("BDH", "PA", 4.0, T0 - 86400, EPOCH),
                 _channel("BDH", "PA", 8.0, EPOCH, T0 + 30 * 86400)]  # fmt: skip
    station = Station(code="OBS", latitude=1, longitude=2, elevation=-3000,
                      site=Site(name="synthetic"), channels=channels)  # fmt: skip
    (tmp_path / "metadata").mkdir()
    Inventory(networks=[Network(code="XX", stations=[station])], source="test").write(
        str(tmp_path / "metadata" / "XX.xml"), format="STATIONXML"
    )
    (tmp_path / "config.py").write_text(CONFIG)
    monkeypatch.setenv("COMPY_RUN_DIR", str(tmp_path))
    monkeypatch.setenv("COMPY_PIPELINE_CONFIG", str(tmp_path / "config.py"))
    monkeypatch.delenv("COMPY_SDS", raising=False)
    monkeypatch.delenv("COMPY_STATIONS", raising=False)
    monkeypatch.syspath_prepend(str(REPO / "pipeline"))
    return tmp_path


@pytest.fixture
def pipeline_module(run_dir):
    """Import a pipeline script against the test run folder."""
    names = ("config", "a_preprocess", "run_inversions")

    def load(name):
        for other in names:
            sys.modules.pop(other, None)
        return importlib.import_module(name)

    yield load
    for other in names:
        sys.modules.pop(other, None)


def _amplitude(stream, channel, start, end):
    trace = stream.select(channel=channel).slice(start, end).merge()[0]
    return float(np.sqrt(2) * np.std(np.asarray(trace.data, float)))


def test_default_archive_is_in_the_run_folder(run_dir, pipeline_module):
    # 3.2.0 pointed at a drive of the author's machine.
    config = pipeline_module("config")
    assert config.SDS == run_dir / "SDS"
    assert "/Volumes/" not in (REPO / "pipeline" / "config.py").read_text()


def test_stage_a_uses_the_response_of_each_epoch(run_dir, pipeline_module):
    for day in range(3):
        _write_day(run_dir, day)
    stage_a = pipeline_module("a_preprocess")
    days = [T0 + k * 86400 for k in range(3)]
    _, written, status = stage_a.process_chunk("OBS", days)
    assert status == "ok" and written
    phys = read(str(run_dir / "work" / "OBS" / "phys" / "*.mseed"))
    before = _amplitude(phys, "BDH", EPOCH - 20000, EPOCH - 4000)
    after = _amplitude(phys, "BDH", EPOCH + 4000, EPOCH + 20000)
    # 3.2.0 removed the response valid at the start of the segment from all of it.
    assert after / before == pytest.approx(0.5, rel=0.01)
    vertical = _amplitude(phys, "BHZ", EPOCH - 20000, EPOCH + 20000)
    assert vertical > 0  # the channel whose response does not change is not split
    qc = json.loads(next((run_dir / "work" / "OBS" / "qc").glob("*.json")).read_text())
    assert qc["response_epochs"] == [["BDH", str(EPOCH)]]


def test_epochs_with_the_same_response_do_not_split(run_dir, pipeline_module):
    stage_a = pipeline_module("a_preprocess")
    channels = [_channel("BHZ", "M/S", 4.0, T0 - 86400, EPOCH),
                _channel("BHZ", "M/S", 4.0, EPOCH, T0 + 30 * 86400)]  # fmt: skip
    station = Station(code="OBS", latitude=1, longitude=2, elevation=-3000,
                      site=Site(name="synthetic"), channels=channels)  # fmt: skip
    inventory = Inventory(networks=[Network(code="XX", stations=[station])], source="t")
    stats = dict(network="XX", station="OBS", location="00", channel="BHZ",
                 sampling_rate=1.25, starttime=T0)  # fmt: skip
    trace = Trace(np.zeros(int(3 * 86400 * 1.25)), stats)
    assert stage_a.epoch_boundaries(inventory, trace) == []


def test_split_at_keeps_every_sample_once(run_dir, pipeline_module):
    stage_a = pipeline_module("a_preprocess")
    trace = Trace(np.arange(100.0), dict(sampling_rate=1.25, starttime=T0))
    pieces = stage_a.split_at(trace, [T0 + 8.0, T0 + 40.3])
    assert [p.stats.npts for p in pieces] == [11, 40, 49]
    assert np.array_equal(np.concatenate([p.data for p in pieces]), trace.data)
    assert pieces[1].stats.starttime == T0 + 8.8  # the sample at 8.0 s stays before


def test_stage_a_leaves_out_a_day_at_another_rate(run_dir, pipeline_module, capfd):
    _write_day(run_dir, 0)
    _write_day(run_dir, 1, rate=5.0)
    _write_day(run_dir, 2)
    stage_a = pipeline_module("a_preprocess")
    _, written, status = stage_a.process_chunk(
        "OBS", [T0 + k * 86400 for k in range(3)]
    )
    assert status == "ok"
    qc = json.loads(next((run_dir / "work" / "OBS" / "qc").glob("*.json")).read_text())
    # 3.2.0 appended the 5 Hz day to the 2.5 Hz trace sample by sample.
    assert sorted(row[0] for row in qc["wrong_rate"]) == sorted(CHANNELS)
    assert qc["expected_rate_hz"] == 2.5
    assert "WRONG-RATE" in capfd.readouterr().err
    names = sorted(Path(p).name for p in written if "/phys/" in p)
    assert names == ["OBS.2020.001.mseed", "OBS.2020.003.mseed"]
    phys = read(str(run_dir / "work" / "OBS" / "phys" / "*.mseed"))
    assert {tr.stats.sampling_rate for tr in phys} == {1.25}
    # the 0.05 Hz line is still at 0.05 Hz on both days kept
    for day in (0, 2):
        tr = phys.select(channel="BHZ").slice(T0 + day * 86400 + 7200,
                                              T0 + day * 86400 + 7200 + 8000 * 0.8)[0]  # fmt: skip
        spectrum = np.abs(np.fft.rfft(tr.data[:8000] * np.hanning(8000)))
        peak = np.fft.rfftfreq(8000, 0.8)[int(np.argmax(spectrum[1:])) + 1]
        assert peak == pytest.approx(0.05, abs=2e-4)


def test_traces_at_different_rates_are_not_joined(run_dir, pipeline_module):
    stage_a = pipeline_module("a_preprocess")
    first = Trace(np.zeros(100), dict(sampling_rate=50.0, starttime=T0))
    second = Trace(np.ones(100), dict(sampling_rate=100.0, starttime=T0 + 2.0))
    assert len(stage_a.join_day_files([first, second])) == 2
    third = Trace(np.ones(100), dict(sampling_rate=50.0, starttime=T0 + 2.0))
    assert len(stage_a.join_day_files([first, third])) == 1


def test_stage_a_stops_when_a_station_has_no_day_file(
    run_dir, pipeline_module, monkeypatch
):
    stage_a = pipeline_module("a_preprocess")
    monkeypatch.setattr(sys, "argv", ["a_preprocess.py", "OBS", "--workers", "1"])
    # 3.2.0 printed "0 days, 0 to do" and ended with exit status 0.
    with pytest.raises(SystemExit) as error:
        stage_a.main()
    assert "no day file for OBS" in str(error.value)
    assert str(run_dir / "SDS") in str(error.value)


def test_inversion_queue_creates_its_log_folder(run_dir, pipeline_module):
    queue = pipeline_module("run_inversions")
    assert not (run_dir / "logs").exists()
    assert (
        queue.run(["-c", "print('ok')"], "probe.log") == 0
    )  # FileNotFoundError in 3.2.0
    assert (run_dir / "logs" / "probe.log").read_text().strip() == "ok"


@pytest.mark.parametrize("module", ["openpyxl", "cartopy"])
def test_commands_name_the_missing_extra(module, monkeypatch):
    from compy_pressure import _commands

    extra = _commands.EXTRAS[module]
    target = (
        "compy_pressure.inventory" if module == "openpyxl" else "compy_pressure.maps"
    )
    for name in [n for n in sys.modules if n == module or n.startswith(module + ".")]:
        monkeypatch.delitem(sys.modules, name)
    monkeypatch.delitem(sys.modules, target, raising=False)
    monkeypatch.setitem(sys.modules, module, None)  # as if it were not installed
    with pytest.raises(SystemExit) as error:
        _commands.run(target)
    assert f"seafloor-compy[{extra}]" in str(error.value)


def test_other_import_errors_are_not_hidden(monkeypatch):
    from compy_pressure import _commands

    with pytest.raises(ModuleNotFoundError):
        _commands.run("compy_pressure.no_such_module")


def test_every_pressure_command_goes_through_the_wrapper():
    text = (REPO / "pyproject.toml").read_text()
    scripts = re.findall(
        r'^(compy-pressure[\w-]*) = "([\w.]+):(\w+)"', text, flags=re.M
    )
    assert len(scripts) == 9
    from compy_pressure import _commands

    for name, module, function in scripts:
        assert module == "compy_pressure._commands", name
        assert callable(getattr(_commands, function)), name


def test_compy_invert_explains_a_csv_without_its_columns(tmp_path):
    from compy.inversion.cli import read_compliance

    bad = tmp_path / "bad.csv"
    bad.write_text("f,eta,err\n0.007,1e-10,1e-12\n0.008,1.1e-10,1e-12\n")
    with pytest.raises(SystemExit) as error:  # StopIteration in 3.2.0
        read_compliance(bad, depth=3000.0)
    assert "freq_hz, compliance and uncertainty" in str(error.value)
    one = tmp_path / "one.csv"
    one.write_text("freq_hz,compliance,uncertainty\n0.007,1e-10,1e-12\n")
    with pytest.raises(SystemExit) as error:  # ValueError about vectors in 3.2.0
        read_compliance(one, depth=3000.0)
    assert "at least three" in str(error.value)
    good = tmp_path / "good.csv"
    good.write_text("freq_hz,compliance_1_per_pa,uncertainty_1_per_pa\n"
                    "0.007,1e-10,1e-12\n0.008,1.1e-10,1e-12\n0.009,1.2e-10,1e-12\n")  # fmt: skip
    freqs, eta, cov, depth, _ = read_compliance(good, depth=3000.0)
    assert len(freqs) == 3 and cov.shape == (3, 3) and depth == 3000.0


def test_compy_invert_refuses_a_zero_gain_sigma(tmp_path, capsys):
    from compy.inversion.cli import main

    curve = tmp_path / "curve.csv"
    curve.write_text(
        "freq_hz,compliance,uncertainty\n0.007,1e-10,1e-12\n0.008,1.1e-10,1e-12\n"
    )
    # 3.2.0 took 0 as "not given" and used the default prior without a word.
    with pytest.raises(SystemExit):
        main([str(curve), "--depth", "3000", "--lat", "0", "--lon", "0",
              "--crust1", str(tmp_path), "--gain-sigma", "0"])  # fmt: skip
    assert "--gain-sigma must be above zero" in capsys.readouterr().err


def test_deprecation_message_does_not_name_a_version_already_released():
    import compy
    from compy import _deprecation

    major = int(compy.__version__.split(".")[0])
    announced = re.search(r"ComPy (\d+)", _deprecation.REMOVAL)
    # 3.2.0 still said "will be removed in ComPy 3.0".
    assert announced is None or int(announced.group(1)) > major
