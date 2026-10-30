import hashlib
import json
from pathlib import Path

import pytest


def week(shift, g1=None, g2=None, label="nominal", sign=-1):
    return {
        "shift": shift,
        "start": None,
        "label": label,
        "sign": sign,
        "g1": g1,
        "g2": g2,
    }


def test_aggregate_classes():
    from scripts.combine_microseism_weeks import aggregate

    stable = aggregate([week(0, 0.69, 0.71), week(6, 0.70, 0.72), week(-6, 0.68, 0.70)])
    assert stable["stability"] == "stable"
    assert stable["median_g"] == pytest.approx(0.70)
    assert stable["spread_pct"] == pytest.approx(100 * 0.02 / 0.70)
    assert (stable["pooled_h1"], stable["pooled_h2"]) == (0.68, 0.72)

    variable = aggregate([week(0, 0.60, 0.62), week(6, 0.70, 0.72)])
    assert variable["stability"] == "variable"

    single = aggregate([week(0, 0.60, 0.62), week(6), week(12, label=None)])
    assert single["stability"] == "single_week"
    assert (single["weeks_checked"], single["weeks_informative"]) == (2, 1)

    none = aggregate([week(0), week(6, label="unexplained_phase", sign=None)])
    assert none["stability"] == "not_informative"
    assert none["sign_stable"] is True

    flipped = aggregate([week(0), week(6, label="reversed_nominal", sign=1)])
    assert flipped["sign_stable"] is False


def test_week_runs_finds_only_runs_made_from_the_base(tmp_path):
    from scripts.combine_microseism_weeks import week_runs

    base = tmp_path / "run_base"
    for name, settings in (
        ("run_base", {"base_run": None}),
        ("run_plus6", {"base_run": str(base), "shift_weeks": 6}),
        ("run_minus6", {"base_run": str(base), "shift_weeks": -6}),
        ("run_other", {"base_run": str(tmp_path / "run_x"), "shift_weeks": 6}),
    ):
        (tmp_path / name).mkdir()
        (tmp_path / name / "settings.json").write_text(json.dumps(settings))
    for name in ("run_plus6", "run_minus6", "run_other"):
        (tmp_path / name / "summary.json").write_text("{}")
        (tmp_path / name / "method_comparison.json").write_text("{}")
    (tmp_path / "run_plus6_again").mkdir()
    (tmp_path / "run_plus6_again/settings.json").write_text(
        json.dumps({"base_run": str(base), "shift_weeks": 6})
    )
    (tmp_path / "run_plus6_again/summary.json").write_text("{}")
    (tmp_path / "run_plus6_again/method_comparison.json").write_text("{}")
    (tmp_path / "run_plus6_unfinished").mkdir()
    (tmp_path / "run_plus6_unfinished/settings.json").write_text(
        json.dumps({"base_run": str(base), "shift_weeks": 6})
    )
    skipped = []
    found = week_runs(base, skipped)
    assert [(s, p.name) for s, p in found] == [
        (-6, "run_minus6"),
        (0, "run_base"),
        (6, "run_plus6_again"),
    ]
    assert [Path(p).name for p in skipped] == ["run_plus6"]
    (tmp_path / "run_plus6_again/method_comparison.json").write_text(
        json.dumps({"rows": [{"station": "XX.S", "claude_g_h1": 999}]})
    )
    with pytest.raises(ValueError, match="Conflicting completed repeats"):
        week_runs(base)


def test_shifted_start_moves_whole_weeks(tmp_path):
    from scripts.run_microseism_check import shifted_start

    (tmp_path / "stations").mkdir()
    (tmp_path / "stations/YV.RR52.json").write_text(
        json.dumps({"start_utc": "2013-06-10T00:00:00.000000Z"})
    )
    (tmp_path / "stations/YV.RR12.json").write_text(json.dumps({"status": "x"}))
    assert shifted_start(tmp_path, "YV", "RR52", -2).startswith("2013-05-27T00:00:00")
    assert shifted_start(tmp_path, "YV", "RR12", 2) is None
    assert shifted_start(tmp_path, "YV", "RR99", 2) is None


def test_published_values_are_frequency_matched():
    from scripts.compare_microseism_methods import published

    def xo(station, band):
        return {"network": "XO", "station": station, "best_sub_band_hz": band}

    assert published(xo("WD66", [0.034, 0.043]))[3] is None  # partial overlap only
    assert published(xo("WD66", [0.032, 0.039]))[3] == 0.92
    assert published(xo("WD66", [0.03, 0.065]))[3] is None  # spans both bands
    assert published(xo("WD48", [0.115, 0.171]))[3] is None  # outside their band
    assert published(xo("WD70", [0.044, 0.058]))[3] is None  # mapping not given
    low, high, _, matched = published(xo("WD70", [0.044, 0.058]))
    assert (low, high, matched) == (0.91, 1.71, None)
    yl = {"network": "YL", "station": "A02W", "published_g": 0.667}
    assert published(yl)[3] == 0.667


def test_combine_sums_counts_and_rejects_duplicate_stations(tmp_path):
    from scripts.compare_microseism_methods import combine

    def run(name, stations, counts):
        folder = tmp_path / name
        folder.mkdir()
        rows = [
            {"station": s, "network": s[:2], "status": "unresolved_pairing"}
            for s in stations
        ]
        (folder / "method_comparison.json").write_text(
            json.dumps(
                {
                    "rows": rows,
                    "run_summary_sha256": name,
                    "columns": {},
                    "phase_label_counts": counts,
                    "informative_stations": [],
                }
            )
        )
        return folder

    docs = tmp_path / "docs"
    (docs / "figures").mkdir(parents=True)
    a = run("a", ["YV.RR01"], {"YV:unresolved_pairing": 1})
    b = run("b", ["YV.RR02"], {"YV:unresolved_pairing": 1})
    _, summary = combine([a, b], docs)
    assert summary["phase_label_counts"] == {"YV:unresolved_pairing": 2}
    with pytest.raises(ValueError, match="same station"):
        combine([a, run("c", ["YV.RR01"], {})], docs)


def test_codex_column_reuses_only_a_matching_saved_diagnostic(tmp_path):
    from scripts.compare_microseism_methods import codex_column

    paired = tmp_path / "comparison.json"
    paired.write_text("{}")
    spectra = tmp_path / "spectra.npz"
    spectra.write_bytes(b"saved input")
    output = tmp_path / "codex"
    output.mkdir()
    saved = {
        "source_comparison": str(paired.resolve()),
        "provenance": {
            "saved_input_artifacts": {
                name: {
                    "path": str(path),
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                }
                for name, path in (("comparison", paired), ("spectra", spectra))
            }
        },
        "results": {
            "fixed_pressure_multiplier": 0.7,
            "median_coherence_squared": 0.9,
            "median_nominal_transfer_phase_deg": 179.0,
            "median_nominal_transfer_amplitude_m_pa": 1e-7,
            "ruan_high_coherence_reference": {"fraction_of_band": 0.5},
            "saikia_daily_reference": {"selected_days": 3, "available_days": 7},
        },
    }
    (output / "diagnostic.json").write_text(json.dumps(saved))
    row = codex_column({}, paired, 0.7, output)
    assert row["codex_phase_family"] == "near_180_deg"
    assert row["codex_selected_days"] == "3/7"
    with pytest.raises(ValueError, match="other inputs"):
        codex_column({}, paired, 0.8, output)
    spectra.write_bytes(b"different data at the same path")
    with pytest.raises(ValueError, match="stale diagnostic"):
        codex_column({}, paired, 0.7, output)
