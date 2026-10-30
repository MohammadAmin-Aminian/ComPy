import numpy as np
import pandas as pd
import pytest
from obspy import Stream, Trace, UTCDateTime

from compy_pressure import tilt


def row(
    seed, start="2012-01-01T00:00:00", end="2013-01-01T00:00:00", status="complete"
):
    return {
        "Seed_id": seed,
        "channel_start_time": start,
        "channel_end_time": end,
        "ppsd_status": status,
        "ppsd_npz_file": f"/r/{seed}/ppsd.npz",
    }


def test_choose_triplets_prefers_one_hertz_band_and_needs_both_horizontals():
    master = pd.DataFrame(
        [
            row("XX.A..BHZ"),
            row("XX.A..BH1"),
            row("XX.A..BH2"),
            row("XX.A..LHZ"),
            row("XX.A..LHN"),
            row("XX.A..LHE"),
            row("XX.B..BHZ"),
            row("XX.B..BH1"),
            row("XX.C..HHZ"),
            row("XX.C..HH1", status="failed"),
            row("XX.C..HH2"),
            row("XX.D..BDH"),
            row("XX.E..EH3"),
            row("XX.E..EH1"),
            row("XX.E..EH2"),
        ]
    )
    chosen = [[r["Seed_id"] for r in t] for t in tilt.choose_triplets(master)]
    assert chosen == [
        ["XX.A..LHZ", "XX.A..LHN", "XX.A..LHE"],
        ["XX.E..EH3", "XX.E..EH1", "XX.E..EH2"],
    ]


def test_to_target_rate_resamples_onto_whole_seconds():
    start = UTCDateTime(2012, 1, 1, 0, 0, 0, 300000)
    trace = Trace(
        np.random.default_rng(0).normal(size=20 * 3600), {"sampling_rate": 20.0}
    )
    trace.stats.starttime = start
    out = tilt.to_target_rate(Stream([trace]))
    assert len(out) == 1
    assert out[0].stats.sampling_rate == 1.0
    assert out[0].stats.starttime == UTCDateTime(2012, 1, 1, 0, 0, 1)


def test_common_intervals_renames_and_aligns_components():
    def trace(code, offset, seconds):
        t = Trace(
            np.zeros(seconds), {"sampling_rate": 1.0, "channel": code, "station": "A"}
        )
        t.stats.starttime = UTCDateTime(2012, 1, 1) + offset
        return Stream([t])

    hours = tilt.WINDOW_HOURS * 3600
    triples = tilt.common_intervals(
        trace("BHZ", 0, hours + 600),
        trace("BHN", 100, hours + 600),
        trace("BHE", 0, hours + 600),
    )
    assert len(triples) == 1
    assert [t.stats.channel for t in triples[0]] == ["BHZ", "BH1", "BH2"]
    assert len({len(t) for t in triples[0]}) == 1
    assert triples[0][0].stats.starttime == UTCDateTime(2012, 1, 1) + 100


def test_band_differences_reports_each_band():
    class Fake:
        period_bin_centers = 1 / np.array([0.005, 0.02, 0.05, 0.2])

        def __init__(self, value):
            self.psd_values = [np.full(4, value)] * 3

    result = tilt.band_differences(Fake(-150.0), Fake(-160.0))
    assert result == pytest.approx({k: -10.0 for k in result})
