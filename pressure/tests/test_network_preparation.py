"""Header-based views preserve samples; APG selection is explicit."""

import hashlib
from importlib import import_module

import numpy as np
import pandas as pd
import pytest
from obspy import Trace


@pytest.mark.parametrize("module", ["prepare_7d_run", "prepare_xo_run"])
def test_location_alias_uses_header_without_changing_source(tmp_path, module):
    helper = import_module(f"scripts.{module}")
    net = "7D" if "7d" in module else "XO"
    trace = Trace(np.arange(1000, dtype=np.int32))
    trace.stats.update(
        dict(
            network=net,
            station="TEST",
            location="",
            channel="LDH",
            sampling_rate=1,
            starttime=0,
        )
    )
    suffix = ".mseed" if net == "XO" else ""
    source = tmp_path / f"{net}.TEST.00.LDH.D.1970.001{suffix}"
    trace.write(str(source), format="MSEED")
    before = hashlib.sha256(source.read_bytes()).hexdigest()
    item = helper.source_header(source)
    assert item["filename_seed_id"] == f"{net}.TEST.00.LDH"
    assert item["seed_id"] == f"{net}.TEST..LDH"
    assert item["sample_rate_hz"] == 1
    assert hashlib.sha256(source.read_bytes()).hexdigest() == before
    bad = tmp_path / f"{net}.OTHER.00.LDH.D.1970.001{suffix}"
    bad.write_bytes(source.read_bytes())
    with pytest.raises(ValueError, match="Network/station/channel mismatch"):
        helper.source_header(bad)


@pytest.mark.parametrize("module", ["prepare_7d_run", "prepare_xo_run"])
def test_mixed_native_rates_remain_rejected(tmp_path, module):
    helper = import_module(f"scripts.{module}")
    net = "7D" if "7d" in module else "XO"
    from obspy import Stream

    traces = []
    for start, rate in [(0, 1), (2000, 2)]:
        trace = Trace(np.arange(1000, dtype=np.int32))
        trace.stats.update(
            dict(
                network=net,
                station="TEST",
                location="",
                channel="LDH",
                sampling_rate=rate,
                starttime=start,
            )
        )
        traces.append(trace)
    source = tmp_path / f"{net}.TEST..LDH.D.1970.001"
    Stream(traces).write(str(source), format="MSEED")
    with pytest.raises(ValueError, match="Mixed identifiers or sample rates"):
        helper.source_header(source)


def test_xo_equal_rate_pair_does_not_resample_125_hz_apg():
    from scripts.prepare_xo_run import choose_pairs

    rows, available = [], {}
    for code, rate, units, dip in [
        ("HDH", 125, "Pa", 0),
        ("EDH", 100, "Pa", 90),
        ("HHZ", 100, "m/s", -90),
    ]:
        seed = f"XO.LA33..{code}"
        rows.append(
            dict(
                Seed_id=seed,
                input_units=units,
                sample_rate_hz=str(rate),
                dip_deg=str(dip),
                channel_start_time="2018-05-01",
                channel_end_time="2019-01-01",
                inventory_file="source.xml",
                Profondeur="2000",
                sensor_description="Test sensor",
            )
        )
        available[seed] = [
            dict(
                sample_rate_hz=rate,
                start_utc="2018-06-01",
                end_utc_exclusive="2018-07-01",
            )
        ]
    selected, plans = choose_pairs(pd.DataFrame(rows), available)
    assert set(selected.Seed_id) == {"XO.LA33..EDH", "XO.LA33..HHZ"}
    assert plans[0]["sample_rate_hz"] == 100
    only_apg = pd.DataFrame(rows).query("Seed_id != 'XO.LA33..EDH'")
    with pytest.raises(ValueError, match="No equal-rate"):
        choose_pairs(only_apg, available)

    apg = pd.DataFrame(rows)
    apg.loc[apg.Seed_id.str.endswith("HDH"), "sensor_description"] = "Parosci APG"
    selected, plans = choose_pairs(apg, available, "apg", 2)
    assert set(selected.Seed_id) == {"XO.LA33..HDH", "XO.LA33..HHZ"}
    assert list(selected.sample_rate_hz) == ["125", "100"]
    assert plans[0]["event_resample_rate_hz"] == 2
    # HDH is also used by high-rate DPGs: hardware labels control APG selection.
    dpg = []
    for code, rate, units, dip in [
        ("HDH", 100, "Pa", 0),
        ("HHZ", 100, "m/s", -90),
        ("LDH", 1, "Pa", 0),
        ("LHZ", 1, "m/s", -90),
    ]:
        row = rows[0].copy()
        row.update(
            Seed_id=f"XO.WD66..{code}",
            sample_rate_hz=str(rate),
            input_units=units,
            dip_deg=str(dip),
            sensor_description="SIO Differential Pressure Gauge",
        )
        dpg.append(row)
        available[row["Seed_id"]] = [
            dict(
                sample_rate_hz=rate,
                start_utc="2018-06-01",
                end_utc_exclusive="2018-07-01",
            )
        ]
    gauges, plans = choose_pairs(
        pd.concat([apg, pd.DataFrame(dpg)]), available, "gauges", 1
    )
    assert set(gauges.Seed_id) == {
        "XO.LA33..HDH",
        "XO.LA33..HHZ",
        "XO.WD66..LDH",
        "XO.WD66..LHZ",
    }
    assert [p["pressure_kind"] for p in plans] == ["apg", "dpg"]
    with pytest.raises(ValueError, match="not upsample"):
        choose_pairs(apg, available, "apg", 126)


def test_apg_mixed_rate_roundoff_is_recorded_without_rewriting(tmp_path, monkeypatch):
    from obspy import Stream

    from scripts.prepare_xo_run import source_header

    traces = []
    for start, rate in [(0, 125), (10, 125.00000762939453)]:
        traces.append(
            Trace(
                np.arange(1250, dtype=np.int32),
                header=dict(
                    network="XO",
                    station="TEST",
                    channel="HDH",
                    location="",
                    sampling_rate=rate,
                    starttime=start,
                ),
            )
        )
    path = tmp_path / "XO.TEST..HDH.D.1970.001"
    Stream(traces).write(str(path), format="MSEED")
    # libmseed may choose integer rate factors when writing this near-integer
    # rate; emulate the mixed blockette-100 headers found in the downloaded APGs.
    monkeypatch.setattr("scripts.prepare_xo_run.read", lambda *a, **k: Stream(traces))
    before = path.read_bytes()
    with pytest.raises(ValueError, match="Mixed identifiers or sample rates"):
        source_header(path)
    result = source_header(path, 0.1)
    assert result["sample_rates_hz"] == [125, 125.00000762939453]
    assert path.read_bytes() == before
