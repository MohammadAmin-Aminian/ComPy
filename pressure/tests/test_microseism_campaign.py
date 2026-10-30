"""Checks for the one-station-at-a-time microseism campaign runner."""

import csv
from collections import Counter
from pathlib import Path

import pytest

from compy_pressure import microseism_campaign as campaign


def test_manifest_matches_all_saved_first_method_gains():
    rows = campaign.load_manifest()
    assert len(rows) == 46
    assert Counter(row["network"] for row in rows) == {
        "YV": 16,
        "YL": 22,
        "7D": 4,
        "XO": 4,
    }
    with (campaign.ROOT / "docs/socrates_all_station_comparison.csv").open(
        newline=""
    ) as handle:
        comparison = list(csv.DictReader(handle))
    expected = {
        f"{row['network']}.{row['station']}": float(row["previous_g"])
        for row in comparison
        if row["previous_g"]
    }
    actual = {
        f"{row['network']}.{row['station']}": float(row["fixed_pressure_multiplier"])
        for row in rows
    }
    assert actual == expected


def test_manifest_channel_and_paper_directions():
    rows = campaign.load_manifest()
    rr52 = campaign.select_station(rows, "YV.RR52")
    assert rr52["pressure_seed_id"] == "YV.RR52.00.BDH"
    assert rr52["vertical_seed_id"] == "YV.RR52.00.BHZ"
    assert rr52["comparison_direction"] == "g"
    assert rr52["published_reference"] == "Aminian2025"
    yl = campaign.select_station(rows, "YL.A02W")
    context = campaign.comparison_context(yl, campaign.DEFAULT_MANIFEST)
    assert context["first_method"]["comparison_coefficient"] == pytest.approx(
        1 / 0.6962199558732434
    )
    assert context["publication"]["coefficient"] == 1.5
    la33 = campaign.select_station(rows, "XO.LA33")
    assert la33["event_resample_rate_hz"] == "1"
    assert la33["native_rate_tolerance_ppm"] == "0.1"


def test_run_station_reuses_designated_comparison_and_fixed_gain(tmp_path, monkeypatch):
    comparison = tmp_path / "comparison.json"
    comparison.write_text(
        '{"pressure_epoch":{"Seed_id":"XX.ONE..LDH"},"vertical_seed_id":"XX.ONE..LHZ"}'
    )
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "network,station,pressure_seed_id,vertical_seed_id,"
        "fixed_pressure_multiplier,comparison_direction,published_coefficient,"
        "published_reference,published_comparison_status,master_dir,data_dir,"
        "event_resample_rate_hz,native_rate_tolerance_ppm,saved_comparison\n"
        f"XX,ONE,XX.ONE..LDH,XX.ONE..LHZ,0.75,g,,,none,unused,unused,,,{comparison}\n"
    )
    captured = {}

    def fake_diagnostic(source, output, *, pressure_multiplier, comparison_context):
        captured.update(
            source=source,
            output=output,
            gain=pressure_multiplier,
            context=comparison_context,
        )
        return Path(output) / "diagnostic.json"

    monkeypatch.setattr(campaign, "run_diagnostic", fake_diagnostic)
    result = campaign.run_station("XX.ONE", manifest_path=manifest)
    assert captured["source"] == comparison
    assert captured["gain"] == 0.75
    assert captured["context"]["microseism_gain_estimated"] is False
    assert result["pressure_seed_id"] == "XX.ONE..LDH"


def test_comparison_binding_refuses_another_station(tmp_path):
    comparison = tmp_path / "comparison.json"
    comparison.write_text(
        '{"pressure_epoch":{"Seed_id":"XX.TWO..LDH"},"vertical_seed_id":"XX.TWO..LHZ"}'
    )
    row = campaign.select_station(campaign.load_manifest(), "YV.RR52")
    with pytest.raises(ValueError, match="do not match"):
        campaign.validate_comparison_binding(comparison, row)
