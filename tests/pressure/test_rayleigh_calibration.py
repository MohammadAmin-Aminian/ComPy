"""Independent physical and gain checks plus real raw-waveform integration."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from obspy import Trace, UTCDateTime
from obspy.core.inventory import Channel, Inventory, Network, Site, Station
from obspy.core.inventory.response import Response
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation

from compy_pressure import rayleigh_calibration as ray
from compy_pressure.inventory import (
    load_master_dataset,
    synchronize_inventory,
    update_processing_results,
)
from compy_pressure.paired_spectra import KEYS


@pytest.mark.parametrize("phase_speed", [4000, 1300, 1500])
def test_water_model_against_acoustic_ode(phase_speed):
    depth, density, sound = 2918, 1028, 1500
    frequency = np.array([0, 0.02, 0.04, 0.06])
    expected = ray.water_ratio(frequency, phase_speed, depth, density, sound)
    for f, value in zip(frequency, expected, strict=True):
        q2 = (2 * np.pi * f) ** 2 * (sound**-2 - phase_speed**-2)
        # P(0)=0, arbitrary dP/dz(0)=1; acoustic ODE gives the floor ratio.
        solution = solve_ivp(
            lambda z, y: [y[1], -q2 * y[0]], [0, depth], [0, 1], rtol=1e-10, atol=1e-12
        )
        pressure, derivative = solution.y[:, -1]
        assert value == pytest.approx(density * pressure / derivative, rel=1e-9)
    assert expected[0] == density * depth
    with pytest.raises(ValueError, match="resonance"):
        ray.water_ratio([0.3], 4000, depth, density, sound)


def test_gain_direction_grid_and_no_zero_or_boundary_solution():
    theory = ray.water_ratio(np.linspace(0.02, 0.06, 9), 4000, 2918, 1028, 1500)
    fitted, grid, loss = ray.fit_multiplier(theory / 0.67, theory)
    assert fitted["pressure_multiplier"] == pytest.approx(0.67)
    assert fitted["sensitivity_counts_per_pa_multiplier"] == pytest.approx(1 / 0.67)
    assert grid[np.argmin(loss)] == pytest.approx(0.67)
    with pytest.raises(ValueError, match="boundary"):
        ray.fit_multiplier(theory / 4, theory)
    with pytest.raises(ValueError, match="positive"):
        ray.fit_multiplier(-theory, theory)
    for a in [np.array([1 + 1j] * 4), np.ma.array(np.ones(4), mask=[0, 1, 0, 0])]:
        with pytest.raises(ValueError):
            ray.fit_multiplier(a, np.ones(4))


def test_group_velocity_controls_arrivals_and_tilt_axis_is_independent():
    event = dict(latitude=0, longitude=90, origin_utc="2013-02-06")
    curve = np.array([[0.02, 4000, 3000], [0.06, 4000, 3000]])
    start, end, distance = ray.event_geometry(event, 0, 0, curve)
    assert (start + 900) - UTCDateTime(event["origin_utc"]) == pytest.approx(
        distance * 1000 / 3000
    )
    curve[:, 1] *= 2
    assert ray.event_geometry(event, 0, 0, curve)[:2] == (start, end)
    rotation = Rotation.from_rotvec([0.04, -0.03, 0]).as_matrix()
    cov = rotation @ np.diag([100, 80, 0.1]) @ rotation.T
    direction, angle = ray.tilt_direction(cov)
    np.testing.assert_allclose(direction, rotation[:, 2], atol=1e-12)
    assert angle == pytest.approx(np.degrees(0.05))
    with pytest.raises(ValueError, match="unresolved"):
        ray.tilt_direction(np.diag([0.1, 80, 100]))


@pytest.fixture
def synthetic_station(tmp_path, request):
    codes = getattr(request, "param", ("BDH", "BH1", "BH2", "BHZ"))
    rate, seconds = 4, 86400
    origin = UTCDateTime("2013-02-06")
    start = origin - 43200
    rng = np.random.default_rng(917)
    n = rate * seconds
    f = np.fft.rfftfreq(n, 1 / rate)
    horizontal = rng.normal(size=(2, n))
    zf = np.fft.rfft(rng.normal(size=n))
    band = (f >= 0.018) & (f <= 0.08)
    zf[~band] = 0
    velocity = np.fft.irfft(zf, n)
    pressure_f = np.zeros(len(f), complex)
    pressure_f[band] = (
        2j
        * np.pi
        * f[band]
        * zf[band]
        * ray.water_ratio(f[band], 4000, 1000, 1028, 1500)
        / 0.67
    )
    pressure = np.fft.irfft(pressure_f, n)
    rotation = Rotation.from_rotvec([0.04, -0.03, 0]).as_matrix()
    seismic = rotation @ np.vstack([horizontal, velocity])
    data, inventories, exports, results = [
        tmp_path / name for name in ["data", "inventories", "exports", "results"]
    ]
    data.mkdir()
    inventories.mkdir()
    arrays = [pressure, *seismic]
    channels = []
    for i, code in enumerate(codes):
        units, gain = ("PA", 4) if i == 0 else ("M/S", 8)
        response = Response.from_paz(
            zeros=[],
            poles=[],
            stage_gain=gain,
            input_units=units,
            output_units="COUNTS",
        )
        response.instrument_sensitivity.frequency = 0.04
        channels.append(
            Channel(
                code=code,
                location_code="00",
                latitude=0,
                longitude=0,
                elevation=-1000,
                depth=0,
                sample_rate=rate,
                dip=-90 if i == 3 else 0,
                response=response,
                start_date=start - 1,
                end_date=start + seconds + 1,
            )
        )
        Trace(
            arrays[i] * gain,
            header=dict(
                network="XX",
                station="TEST",
                location="00",
                channel=code,
                starttime=start,
                sampling_rate=rate,
            ),
        ).write(str(data / f"XX.TEST.00.{code}.D.2013.037"), format="MSEED")
    inv = Inventory(
        [
            Network(
                "XX",
                stations=[
                    Station(
                        "TEST",
                        latitude=0,
                        longitude=0,
                        elevation=-1000,
                        site=Site(name="Test"),
                        channels=channels,
                    )
                ],
            )
        ],
        "Test",
    )
    inv.write(str(inventories / "XX.xml"), format="STATIONXML")
    synchronize_inventory(inventories, exports)
    master = load_master_dataset(exports)
    patches = []
    for i, row in master.iterrows():
        folder = results / str(i)
        folder.mkdir(parents=True)
        inv.select(channel=row.Seed_id.split(".")[-1]).write(
            str(folder / "inventory.xml"), format="STATIONXML"
        )
        (folder / "ppsd.npz").write_bytes(b"existing spectrum")
        patches.append(
            {k: row[k] for k in KEYS}
            | dict(ppsd_status="complete", ppsd_npz_file=str(folder / "ppsd.npz"))
        )
    update_processing_results(exports, patches)
    event = dict(
        id="test_event",
        properties=dict(
            time=float(origin) * 1000,
            mag=8,
            magType="mww",
            url="https://example.invalid",
        ),
        geometry=dict(coordinates=[60, 0, 20]),
    )
    catalog = tmp_path / "catalog.json"
    catalog.write_text(json.dumps(dict(features=[event])))
    dispersion = tmp_path / "dispersion.csv"
    np.savetxt(
        dispersion,
        [[0.015, 4000, 3800], [0.085, 4000, 3800]],
        delimiter=",",
        header="frequency_hz,phase_velocity_m_s,group_velocity_m_s",
        comments="",
    )
    return data, exports, catalog, dispersion


def test_calibration_can_use_response_snapshots_without_ppsd(synthetic_station):
    data, exports, catalog, dispersion = synthetic_station
    master = load_master_dataset(exports)
    master = master[~master.Seed_id.str.endswith(("BH1", "BH2"))].copy()
    master["rayleigh_response_snapshot_file"] = master.ppsd_npz_file.map(
        lambda path: str(Path(path).parent / "inventory.xml")
    )
    for path in master.ppsd_npz_file:
        Path(path).unlink()
    master["ppsd_npz_file"] = ""
    master["ppsd_status"] = ""
    result = ray.run_calibration(
        catalog=catalog,
        dispersion_file=dispersion,
        master_dataset=master,
        data_dir=data,
        network="XX",
        station="TEST",
        tilt_method="none",
        pressure_sign=1,
        event_resample_rate=2,
    )
    assert result["events"][0]["status"] == "candidate"
    assert result["events"][0]["pressure_multiplier"] == pytest.approx(0.67, rel=0.03)
    assert result["settings"]["event_resample_rate_hz"] == 2
    assert Path(result["files"]["result"]).is_file()


def test_result_records_code_hash_from_import_not_write_time(
    synthetic_station, monkeypatch
):
    data, exports, catalog, dispersion = synthetic_station
    assert (
        ray.CODE_SHA256 == hashlib.sha256(Path(ray.__file__).read_bytes()).hexdigest()
    )
    # An edit while a run executes the loaded code must not change its record.
    monkeypatch.setattr(ray, "CODE_SHA256", "code-as-imported")
    result = ray.run_calibration(
        catalog=catalog,
        dispersion_file=dispersion,
        master_dataset=load_master_dataset(exports),
        data_dir=data,
        network="XX",
        station="TEST",
        tilt_method="none",
        pressure_sign=1,
        event_resample_rate=2,
    )
    assert result["code_sha256"] == "code-as-imported"
    saved = json.loads(Path(result["files"]["result"]).read_text())
    assert saved["code_sha256"] == "code-as-imported"


@pytest.mark.parametrize("downward", [False, True])
def test_no_tilt_uses_only_event_pressure_and_vertical(
    synthetic_station, monkeypatch, downward
):
    from obspy import read, read_inventory

    data, exports, catalog, dispersion = synthetic_station
    traces = {p.name.split(".")[3]: (p, read(str(p))[0]) for p in data.iterdir()}
    rotation = Rotation.from_rotvec([0.04, -0.03, 0]).as_matrix()
    # Recover the known original vertical, independent of the production projection.
    vertical = (
        rotation.T @ np.vstack([traces[c][1].data for c in ("BH1", "BH2", "BHZ")])
    )[2]
    traces["BHZ"][1].data = vertical * (-1 if downward else 1)
    origin = UTCDateTime("2013-02-06")
    for code, (path, trace) in traces.items():
        if code in {"BH1", "BH2"}:
            path.unlink()
        else:
            # No preceding quiet interval exists, but the event/overview are covered.
            trace.trim(starttime=origin)
            trace.write(str(path), format="MSEED")
    master = load_master_dataset(exports)
    master = master[~master.Seed_id.str.endswith(("BH1", "BH2"))]
    if downward:
        z = master[master.Seed_id.str.endswith("BHZ")].iloc[0]
        for path in (
            Path(z.inventory_file),
            Path(z.ppsd_npz_file).parent / "inventory.xml",
        ):
            inv = read_inventory(str(path))
            for n in inv:
                for st in n:
                    for c in st:
                        if c.code == "BHZ":
                            c.dip = 90
            inv.write(str(path), format="STATIONXML")
    protected = {
        p: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [*data.iterdir(), *exports.glob("master_inventory.*")]
    }
    monkeypatch.setattr(
        ray, "require_clean_rotator", lambda: pytest.fail("TiSKitPy was required")
    )
    monkeypatch.setattr(
        ray, "estimate_tilt", lambda *a: pytest.fail("Tilt was estimated")
    )
    result = ray.run_calibration(
        catalog=catalog,
        dispersion_file=dispersion,
        master_dir=exports,
        master_dataset=master,
        data_dir=data,
        network="XX",
        station="TEST",
        max_events=1,
        pressure_sign=1,
        tilt_method="none",
    )
    item = result["events"][0]
    assert item["status"] == "candidate"
    assert item["pressure_multiplier"] == pytest.approx(0.67, rel=0.02)
    assert item["quiet_start_utc"] is None
    assert item["tilt_angle_deg"] is None
    assert item["tilt"]["applied"] is False
    assert item["rotation_coefficients"] == [1.0]
    assert len(result["epoch_keys"]) == 2
    assert result["settings"]["quiet_window_seconds"] is None
    assert not result["settings"]["tilt_correction_applied"]
    assert not {f["id"] for f in result["review_flags"]} & {
        "tilt_stationarity_needs_testing",
        "quiet_screen_and_coda_heuristic",
        "sensor_axes_assumed_orthogonal",
    }
    with np.load(item["files"]["spectra"]) as spectra:
        np.testing.assert_allclose(spectra["h1_pa"], spectra["raw_h1_pa"])
        np.testing.assert_allclose(
            spectra["coherence_squared"], spectra["raw_coherence_squared"]
        )
    assert Path(item["files"]["arrival_review"]).is_file()
    assert Path(result["files"]["paper_comparison"]).is_file()
    assert protected == {
        p: hashlib.sha256(p.read_bytes()).hexdigest() for p in protected
    }


def test_raw_response_rotation_fit_preservation_and_phase_guard(
    synthetic_station, monkeypatch
):
    from obspy import read

    data, exports, catalog, dispersion = synthetic_station
    originals = {
        p: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in exports.glob("master_inventory.*")
    }
    raw = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in data.iterdir()}
    result = ray.run_calibration(
        catalog=catalog,
        dispersion_file=dispersion,
        master_dir=exports,
        data_dir=data,
        network="XX",
        station="TEST",
        max_events=1,
        pressure_sign=1,
    )
    item = result["events"][0]
    assert item["status"] == "candidate"
    assert item["tilt"]["method"] == "tiskitpy.CleanRotator"
    assert item["tilt"]["package_version"] == "2.3.1"
    assert item["tilt"]["training_rate_hz"] == 1
    assert "arrival_review" in item["files"]
    assert Path(result["files"]["paper_comparison"]).is_file()
    assert result["paper_comparison"]["event_ids"] == ["test_event"]
    assert result["paper_comparison"]["candidate_pressure_multiplier"] == pytest.approx(
        0.67, rel=0.02
    )
    assert item["pressure_multiplier"] == pytest.approx(0.67, rel=0.02)
    assert item["tilt_angle_deg"] == pytest.approx(np.degrees(0.05), abs=0.05)
    assert (
        result["status"] == "unresolved"
    )  # One event cannot establish station calibration.
    assert raw == {
        p: hashlib.sha256(p.read_bytes()).hexdigest() for p in data.iterdir()
    }
    original_spectra = ray.event_spectra

    def reject_display_band(*args, **kwargs):
        if args[6][0] == 0.005:
            raise ValueError("Synthetic response unavailable outside fit band")
        return original_spectra(*args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(ray, "event_spectra", reject_display_band)
        restricted = ray.run_calibration(
            catalog=catalog,
            dispersion_file=dispersion,
            master_dir=exports,
            data_dir=data,
            network="XX",
            station="TEST",
            max_events=1,
            pressure_sign=1,
        )
    assert restricted["events"][0]["status"] == "candidate"
    assert "paper_display_reason" in restricted["events"][0]
    assert restricted["paper_comparison"]["candidate_pressure_multiplier"] == (
        pytest.approx(0.67, rel=0.02)
    )
    pressure = next(p for p in data.iterdir() if "BDH" in p.name)
    stream = read(str(pressure))
    stream[0].data *= -1
    stream.write(str(pressure), format="MSEED")
    reversed_result = ray.run_calibration(
        catalog=catalog,
        dispersion_file=dispersion,
        master_dir=exports,
        data_dir=data,
        network="XX",
        station="TEST",
        max_events=1,
        pressure_sign=1,
    )
    reversed_item = reversed_result["events"][0]
    assert reversed_item["status"] == "polarity_review"
    assert "pressure_multiplier" not in reversed_item
    assert reversed_item["diagnostic_pressure_amplitude_multiplier"] == pytest.approx(
        0.67, rel=0.02
    )
    assert originals == {
        p: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in exports.glob("master_inventory.*")
    }
    with np.load(item["files"]["spectra"]) as saved:
        assert saved["paper_frequency_hz"][0] == pytest.approx(0.005)
        assert saved["paper_frequency_hz"][-1] == pytest.approx(0.1)
        assert saved["coherence_squared"].mean() > 0.95
        assert np.abs(saved["phase_deg"]).max() < 20
        np.testing.assert_allclose(
            abs(saved["cross_pa"]) / saved["acceleration_psd"], abs(saved["h1_pa"])
        )


def test_untrusted_catalog_identifier_cannot_escape_output(tmp_path):
    catalog = tmp_path / "catalog.json"
    catalog.write_text(json.dumps(dict(features=[dict(id="../../existing_result")])))
    with pytest.raises(ValueError, match="identifier"):
        ray.catalog_events(catalog)


def test_clean_rotator_projection_matches_actual_apply_and_preserves_input():
    from obspy import Stream
    from tiskitpy import CleanRotator

    rng = np.random.default_rng(882)
    earth = rng.normal(size=(3, 7200))
    earth[2] *= 0.002
    rotation = Rotation.from_rotvec([0.04, -0.03, 0]).as_matrix()
    instrument = rotation @ earth
    stream = Stream(
        [
            Trace(
                instrument[i].copy(),
                header=dict(
                    network="XX",
                    station="TEST",
                    location="00",
                    channel=code,
                    sampling_rate=1,
                ),
            )
            for i, code in enumerate(["BH1", "BH2", "BHZ"])
        ]
    )
    before = [t.data.copy() for t in stream]
    direction, angle, ratio, info = ray.clean_rotator_axis(stream)
    direct = CleanRotator(
        stream,
        remove_eqs=False,
        save_eq_file=False,
        plot=False,
        verbose=False,
        filt_band=(0.001, 0.01),
        H_over_Z=1,
        N_E_comps=("1", "2"),
    )
    applied = direct.apply(stream, horiz_too=False, rot_limit=10, set_dtype=False)
    np.testing.assert_allclose(
        direction @ instrument,
        applied.select(component="Z")[0].data,
        rtol=1e-12,
        atol=1e-12,
    )
    assert angle == pytest.approx(np.degrees(0.05), abs=0.05)
    assert ratio < 0.01
    assert info["method"] == "tiskitpy.CleanRotator"
    for old, new in zip(before, stream, strict=True):
        np.testing.assert_array_equal(old, new.data)


def test_local_method_is_explicit_and_flagged(synthetic_station):
    data, exports, catalog, dispersion = synthetic_station
    result = ray.run_calibration(
        catalog=catalog,
        dispersion_file=dispersion,
        master_dir=exports,
        data_dir=data,
        network="XX",
        station="TEST",
        max_events=1,
        pressure_sign=1,
        tilt_method="local_eigen",
    )
    item = result["events"][0]
    assert item["tilt"]["validation_status"] == "needs_to_be_tested"
    assert item["tilt"]["method"] == "local_eigen"
    assert item["pressure_multiplier"] == pytest.approx(0.67, rel=0.02)


def test_explicit_assumed_sign_handles_opposite_polarity(synthetic_station):
    from obspy import read

    data, exports, catalog, dispersion = synthetic_station
    pressure = next(p for p in data.iterdir() if "BDH" in p.name)
    stream = read(str(pressure))
    stream[0].data *= -1
    stream.write(str(pressure), format="MSEED")
    result = ray.run_calibration(
        catalog=catalog,
        dispersion_file=dispersion,
        master_dir=exports,
        data_dir=data,
        network="XX",
        station="TEST",
        max_events=1,
        pressure_sign=-1,
    )
    item = result["events"][0]
    assert item["status"] == "candidate"
    assert item["signed_pressure_multiplier"] == pytest.approx(-0.67, rel=0.02)
    assert result["pressure_sign_assumption"]["status"] == "assumed_needs_confirmation"
    with np.load(item["files"]["spectra"]) as saved:
        assert np.abs(saved["phase_deg"]).max() < 20
        assert np.abs(saved["nominal_phase_deg"]).min() > 160


def test_timing_error_reports_actual_offsets(synthetic_station):
    from obspy import read

    from compy_pressure.archive import discover_waveforms
    from compy_pressure.ppsd import select_epoch

    data, exports, _, _ = synthetic_station
    master = load_master_dataset(exports)
    rows = [
        master[master.Seed_id.str.endswith(code)].iloc[0]
        for code in ["BDH", "BH1", "BH2", "BHZ"]
    ]
    channels = [select_epoch(r)[1] for r in rows]
    path = next(p for p in data.iterdir() if "BH2" in p.name)
    stream = read(str(path))
    stream[0].stats.starttime += 0.001
    stream.write(str(path), format="MSEED")
    with pytest.raises(ray.SampleAlignmentError) as err:
        ray.load_aligned(
            discover_waveforms(data), rows, channels, UTCDateTime("2013-02-06"), 1800
        )
    assert err.value.span_seconds == pytest.approx(0.001, abs=1e-6)
    assert max(t["offset_samples"] for t in err.value.timing) == pytest.approx(
        0.004, abs=1e-6
    )


def test_manual_window_decisions_do_not_calibrate(synthetic_station, tmp_path):
    data, exports, catalog, dispersion = synthetic_station
    override = tmp_path / "windows.json"
    override.write_text(
        json.dumps(
            {"test_event": {"decision": "exclude", "note": "overlapping arrival"}}
        )
    )
    result = ray.run_calibration(
        catalog=catalog,
        dispersion_file=dispersion,
        master_dir=exports,
        data_dir=data,
        network="XX",
        station="TEST",
        pressure_sign=1,
        window_overrides=override,
    )
    assert result["events"][0]["status"] == "excluded_manual"
    assert result["candidate_events"] == 0
    start = UTCDateTime("2013-02-06") + 1000
    override.write_text(
        json.dumps(
            {
                "test_event": {
                    "decision": "use",
                    "start_utc": str(start),
                    "end_utc": str(start + 2000),
                    "note": "Synthetic known-transfer window with ten segments",
                }
            }
        )
    )
    reviewed = ray.run_calibration(
        catalog=catalog,
        dispersion_file=dispersion,
        master_dir=exports,
        data_dir=data,
        network="XX",
        station="TEST",
        pressure_sign=1,
        window_overrides=override,
    )
    item = reviewed["events"][0]
    assert item["window_review"] == "reviewed_override"
    assert item["start_utc"] == str(start)
    assert item["segment_count"] == 10
    assert item["pressure_multiplier"] == pytest.approx(0.67, rel=0.02)
    assert reviewed["status"] == "unresolved"
    assert (
        reviewed["inputs"][str(override)]
        == hashlib.sha256(override.read_bytes()).hexdigest()
    )


def test_default_sampling_exhausts_catalogue_after_three_successes_and_attempt_cap():
    origin = UTCDateTime("2013-01-01")
    events = [
        dict(id=f"large{i}", magnitude=7.5, origin_utc=str(origin + i * 172800))
        for i in range(10)
    ] + [
        dict(id="medium", magnitude=6.5, origin_utc=str(origin + 10 * 172800)),
        dict(id="six", magnitude=6.0, origin_utc=str(origin + 11 * 172800)),
        dict(id="cluster", magnitude=6.1, origin_utc=str(origin + 11 * 172800 + 3600)),
    ]
    events.append(events[0].copy())
    summary = dict(events=[], magnitude_search=[])
    attempted = []
    for e in ray.magnitude_schedule(events, summary, 1, 3, True, None):
        attempted.append(e["id"])
        summary["events"].append(e | dict(status="candidate"))
    assert attempted == [f"large{i}" for i in range(10)] + ["medium", "six"]
    assert len([e for e in summary["events"] if e["status"] == "candidate"]) == 12
    assert [s["outcome"] for s in summary["magnitude_search"]] == [
        "eligible_events_exhausted",
        "eligible_events_exhausted",
        "eligible_events_exhausted",
    ]
    assert (
        next(e for e in summary["events"] if e["id"] == "cluster")["status"]
        == "excluded_cluster"
    )


def test_all_sampling_can_restrict_to_large_events_without_early_stop():
    events = [
        dict(
            id=f"large{i}",
            magnitude=7.5,
            origin_utc=str(UTCDateTime("2013-01-01") + i * 172800),
        )
        for i in range(5)
    ] + [dict(id="small", magnitude=6.5, origin_utc="2013-02-01")]
    summary = dict(events=[], magnitude_search=[])
    attempted = []
    for e in ray.magnitude_schedule(events, summary, 1, 3, False, None):
        attempted.append(e["id"])
        summary["events"].append(e | dict(status="candidate"))
    assert attempted == [f"large{i}" for i in range(5)]
    assert len(summary["magnitude_search"]) == 1


def test_magnitude_fallback_uses_qc_count_and_stops_at_target():
    events = [
        dict(id="large1", magnitude=7.5, origin_utc="2013-01-01"),
        dict(id="large2", magnitude=7.2, origin_utc="2013-01-03"),
        dict(id="large3", magnitude=7.1, origin_utc="2013-01-05"),
        dict(id="medium", magnitude=6.5, origin_utc="2013-01-07"),
        dict(id="small", magnitude=6.0, origin_utc="2013-01-09"),
    ]
    summary = dict(events=[], magnitude_search=[])
    attempted = []
    events.append(events[0].copy())  # Duplicate catalogue entries are not retries.
    for e in ray.magnitude_schedule(
        events, summary, 8, 3, True, None, sampling="bounded"
    ):
        attempted.append(e["id"])
        # Plenty of large earthquakes exist, but one fails scientific QC.
        summary["events"].append(
            e | dict(status="rejected_qc" if e["id"] == "large2" else "candidate")
        )
    assert attempted == ["large1", "large2", "large3", "medium"]
    assert [s["minimum_magnitude"] for s in summary["magnitude_search"]] == [7, 6.5]
    assert summary["magnitude_search"][-1]["outcome"] == "target_reached"


def test_magnitude_fallback_reaches_six_without_repeating_or_using_clusters():
    events = [
        dict(id="medium", magnitude=6.5, origin_utc="2013-01-01"),
        dict(id="cluster", magnitude=6.4, origin_utc="2013-01-01T12:00:00"),
        dict(id="six1", magnitude=6.0, origin_utc="2013-01-03"),
        dict(id="six2", magnitude=6.2, origin_utc="2013-01-05"),
        dict(id="six3", magnitude=6.1, origin_utc="2013-01-07"),
    ]
    summary = dict(events=[], magnitude_search=[])
    attempted = []
    for e in ray.magnitude_schedule(
        events, summary, 8, 3, True, None, sampling="bounded"
    ):
        attempted.append(e["id"])
        summary["events"].append(
            e | dict(status="rejected_qc" if e["id"] == "medium" else "candidate")
        )
    assert attempted == ["medium", "six1", "six2", "six3"]
    assert [s["minimum_magnitude"] for s in summary["magnitude_search"]] == [7, 6.5, 6]
    assert (
        next(e for e in summary["events"] if e["id"] == "cluster")["status"]
        == "excluded_cluster"
    )


def test_fallback_disabled_or_explicit_small_event_does_not_expand_search():
    events = [dict(id="small", magnitude=6.0, origin_utc="2013-01-01")]
    summary = dict(events=[], magnitude_search=[])
    assert list(ray.magnitude_schedule(events, summary, 8, 3, False, None)) == []
    assert len(summary["magnitude_search"]) == 1
    summary = dict(events=[], magnitude_search=[])
    selected = list(ray.magnitude_schedule(events, summary, 1, 3, True, "small"))
    assert [e["id"] for e in selected] == ["small"]
    assert len(summary["magnitude_search"]) == 1


@pytest.mark.parametrize(
    "synthetic_station", [("LDH", "LH1", "LH2", "LHZ")], indirect=True
)
def test_tilt_horizontals_follow_the_paired_vertical_band(synthetic_station):
    data, exports, catalog, dispersion = synthetic_station
    master = load_master_dataset(exports)
    rows, _ = ray.select_calibration_channels(
        master, "XX", "TEST", tilt_method="local_eigen"
    )
    assert [r.Seed_id.split(".")[-1] for r in rows] == ["LDH", "LH1", "LH2", "LHZ"]
    result = ray.run_calibration(
        catalog=catalog,
        dispersion_file=dispersion,
        master_dir=exports,
        data_dir=data,
        network="XX",
        station="TEST",
        max_events=1,
        pressure_sign=1,
    )
    item = result["events"][0]
    assert item["status"] == "candidate"
    assert item["tilt"]["method"] == "tiskitpy.CleanRotator"
    assert item["pressure_multiplier"] == pytest.approx(0.67, rel=0.02)


@pytest.mark.parametrize("polarity", [1, -1])
def test_boundary_fit_keeps_spectra_and_actual_qc_outcome(synthetic_station, polarity):
    from obspy import read

    data, exports, catalog, dispersion = synthetic_station
    pressure = next(p for p in data.iterdir() if "BDH" in p.name)
    stream = read(str(pressure))
    # g = 0.67 * 5 lies beyond the 0.1-3 grid; -1 also fails the phase guard.
    stream[0].data = np.round(stream[0].data * polarity / 5).astype(np.int32)
    stream.write(str(pressure), format="MSEED")
    result = ray.run_calibration(
        catalog=catalog,
        dispersion_file=dispersion,
        master_dir=exports,
        data_dir=data,
        network="XX",
        station="TEST",
        pressure_sign=1,
        tilt_method="none",
    )
    item = result["events"][0]
    assert item["median_coherence_squared"] > 0.95
    assert "pressure_multiplier" not in item
    assert "diagnostic_pressure_amplitude_multiplier" not in item
    assert "boundary" in item["fit_reason"]
    assert item["status"] == ("unresolved" if polarity == 1 else "polarity_review")
    assert result["candidate_events"] == 0
    with np.load(item["files"]["spectra"]) as saved:
        observed = saved["psd_ratio_pa"] / saved["expected_pa"]
        assert np.median(observed) == pytest.approx(1 / (0.67 * 5), rel=0.02)
        if polarity == -1:
            assert np.all(abs(saved["phase_deg"]) > 160)
    assert Path(item["files"]["plot"]).is_file()


@pytest.mark.parametrize(
    "diagnostic,call_number", [("h1", 2), ("raw", 3), ("paper", 4)]
)
def test_diagnostic_fit_failure_preserves_primary_candidate(
    synthetic_station, monkeypatch, diagnostic, call_number
):
    data, exports, catalog, dispersion = synthetic_station
    original = ray.fit_multiplier
    calls = 0

    def fail_diagnostic(observed, expected):
        nonlocal calls
        calls += 1
        if calls == call_number:
            raise ValueError("Diagnostic minimum lies on a boundary")
        return original(observed, expected)

    monkeypatch.setattr(ray, "fit_multiplier", fail_diagnostic)
    result = ray.run_calibration(
        catalog=catalog,
        dispersion_file=dispersion,
        master_dir=exports,
        data_dir=data,
        network="XX",
        station="TEST",
        pressure_sign=1,
        tilt_method="none",
    )
    item = result["events"][0]
    assert item["status"] == "candidate"
    assert 0.1 < item["pressure_multiplier"] < 3
    assert Path(item["files"]["spectra"]).is_file()
    if diagnostic == "paper":
        assert "boundary" in result["paper_comparison"]["fit_reason"]
        assert result["paper_comparison"]["candidate_pressure_multiplier"] is None
        assert Path(result["files"]["paper_comparison"]).is_file()
    else:
        assert "boundary" in item[f"diagnostic_{diagnostic}_fit_reason"]


def test_event_spectra_phase_and_ratio_from_time_domain_derivative(synthetic_station):
    _, exports, _, _ = synthetic_station
    master = load_master_dataset(exports)
    rows, channels = ray.select_calibration_channels(
        master, "XX", "TEST", tilt_method="none"
    )
    rate = 4
    t = np.arange(rate * 1800) / rate
    # Whole cycles per 200 s segment, every other FFT bin: Hann leakage and the
    # taper-derivative term vanish at these bins, isolating the conventions.
    f0 = np.arange(4, 15, 2) / 200
    phase = np.random.default_rng(31).uniform(0, 2 * np.pi, f0.size)
    velocity = np.cos(2 * np.pi * f0 * t[:, None] + phase).sum(axis=1)
    # Analytic time derivative, independent of the code's i*omega convention.
    acceleration = (-2 * np.pi * f0 * np.sin(2 * np.pi * f0 * t[:, None] + phase)).sum(
        axis=1
    )
    transfer = 2.5e6  # Positive Pa per upward m/s^2
    arrays = np.array([4 * transfer * acceleration, 8 * velocity])
    spectra = ray.event_spectra(
        arrays, rows, channels, rate, np.array([1.0]), 200, (0.02, 0.07)
    )
    at = np.isin(np.round(spectra["frequency_hz"] * 200), np.round(f0 * 200))
    assert at.sum() == f0.size
    np.testing.assert_allclose(spectra["psd_ratio_pa"][at], transfer, rtol=1e-9)
    np.testing.assert_allclose(spectra["phase_deg"][at], 0, atol=1e-6)
    np.testing.assert_allclose(spectra["coherence_squared"][at], 1, atol=1e-9)
