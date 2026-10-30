"""Physical and provenance regressions for the frozen handoff rule."""

import hashlib
import json

import pytest

from compy_pressure.microseism_validation import RULE, frequency_station, frequency_week


def band(low, high, gain=0.7, width=0.01, conditioning=1.2):
    g1, g2 = gain - width / 2, gain + width / 2
    return {
        "frequency_hz": [low, high],
        "mode0_water_h1": g1,
        "mode0_water_h2": g2,
        "upper_bound_h2": g1 * conditioning,
        "median_coherence_squared": 0.98,
    }


def check(*bands):
    return {"long_wave_sub_bands": list(bands)}


def test_high_frequency_model_artifact_cannot_be_selected():
    result = frequency_week(
        check(band(0.03, 0.05), band(0.05, 0.075), band(0.14, 0.20, gain=1.4)),
        "nominal",
    )
    assert result["passes"]
    assert result["conditional_week_g"] == pytest.approx(0.7)
    assert len(result["eligible_sub_bands"]) == 2
    assert result["g_h1"] == pytest.approx(0.695)
    assert result["g_h2"] == pytest.approx(0.705)


def test_a_single_good_band_cannot_pass():
    result = frequency_week(
        check(band(0.03, 0.05), band(0.05, 0.075, conditioning=1.5001)), "nominal"
    )
    assert not result["passes"]
    assert result["conditional_week_g"] is None
    assert "fewer_than_two_eligible_sub_bands" in result["reasons"]


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), -1, 0])
def test_nonphysical_or_missing_values_are_ineligible(value):
    bad = band(0.05, 0.075)
    bad["upper_bound_h2"] = value
    assert not frequency_week(check(band(0.03, 0.05), bad), "nominal")["passes"]


def test_frozen_thresholds_include_the_boundary_and_reject_above_it():
    assert RULE["maximum_sub_band_edge_hz"] == 0.10
    assert RULE["maximum_conditioning_ratio"] == 1.5
    result = frequency_week(
        check(
            band(0.03, 0.05, gain=1.0, conditioning=1.5),
            band(0.05, 0.10, gain=1.1, conditioning=1.5),
        ),
        "nominal",
    )
    assert result["passes"]
    above = frequency_week(
        check(band(0.03, 0.05, gain=1.0), band(0.05, 0.10, gain=1.10001)),
        "nominal",
    )
    assert not above["passes"] and "sub_bands_disagree" in above["reasons"]


@pytest.mark.parametrize(
    "label", ["approximate_nominal", "unexplained_phase", "insufficient_coherent_bins"]
)
def test_phase_failure_is_not_rescued_by_amplitude_agreement(label):
    result = frequency_week(check(band(0.03, 0.05), band(0.05, 0.075)), label)
    assert not result["passes"]
    assert "phase_not_explained_exactly" in result["reasons"]


def test_the_reference_gain_cannot_influence_selection():
    data = check(band(0.03, 0.05), band(0.05, 0.075))
    data["reference"] = {"pressure_multiplier": 0.7}
    first = frequency_week(data, "nominal")
    data["reference"]["pressure_multiplier"] = 100
    assert frequency_week(data, "nominal") == first


def test_repeated_bands_cannot_fake_frequency_consistency():
    data = band(0.03, 0.05)
    with pytest.raises(ValueError, match="overlap or repeat"):
        frequency_week(check(data, data), "nominal")


def week(gain):
    return {
        "sign": -1,
        "frequency_check": frequency_week(
            check(band(0.03, 0.05, gain=gain), band(0.05, 0.075, gain=gain)),
            "nominal",
        ),
    }


def test_three_weeks_are_required_and_twofold_states_remain_unresolved():
    assert not frequency_station([week(0.7), week(0.71)])["conditional_reportable"]
    stable = frequency_station([week(0.7), week(0.71), week(0.69)])
    assert stable["conditional_reportable"] and stable["status"] == "stable"
    assert stable["median_g"] == pytest.approx(0.7)
    variable = frequency_station([week(0.28), week(0.29), week(0.56)])
    assert not variable["conditional_reportable"]
    assert variable["status"] == "variable"
    assert variable["week_spread_pct"] > 90


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return path


@pytest.fixture
def saved_weeks(tmp_path):
    """Three independent saved weeks, including exact input-hash bindings."""
    base = tmp_path / "run_base"
    paths = []
    for shift, start, end in [(0, "01", "08"), (1, "08", "15"), (2, "15", "22")]:
        run = base if shift == 0 else tmp_path / f"run_shift{shift}"
        write_json(
            run / "settings.json",
            {
                "base_run": str(base) if shift else None,
                "shift_weeks": shift,
                "networks": ["XX"],
                "days": 7,
            },
        )
        paired = write_json(
            run / "paired.json",
            {"pressure_epoch": {"id": "P"}, "vertical_epoch": {"id": "Z"}},
        )
        data = check(band(0.03, 0.05), band(0.05, 0.075)) | {
            "station": "XX.S",
            "water": {"depth_m": 3000},
            "model": {"layers_sha256": "same-model"},
            "interval_utc": [f"2020-01-{start}T00:00:00Z", f"2020-01-{end}T00:00:00Z"],
            "inputs": {
                str(paired): {"sha256": hashlib.sha256(paired.read_bytes()).hexdigest()}
            },
        }
        check_path = write_json(run / "check.json", data)
        summary = write_json(
            run / "summary.json",
            {
                "stations": [
                    {
                        "network": "XX",
                        "station": "S",
                        "status": "checked",
                        "check_summary": str(check_path),
                        "paired_result": str(paired),
                    },
                ]
            },
        )
        write_json(
            run / "method_comparison.json",
            {
                "run_summary_sha256": hashlib.sha256(summary.read_bytes()).hexdigest(),
                "rows": [
                    {
                        "station": "XX.S",
                        "network": "XX",
                        "status": "checked",
                        "interval_start": f"2020-01-{start}",
                        "claude_phase_label": "nominal",
                        "claude_implied_sign": -1,
                        "claude_g_h1": 0.695,
                        "claude_g_h2": 0.705,
                    }
                ],
            },
        )
        paths.append(run)
    return base, paths


def test_saved_review_checks_provenance_and_temporal_independence(saved_weeks):
    from scripts.combine_microseism_weeks import InputGuard, load_weeks, summarize

    base, runs = saved_weeks
    guard = InputGuard()
    stations, _ = load_weeks([base], guard=guard)
    result = summarize(stations)[0]["frequency_consistency"]
    assert result["conditional_reportable"] and result["independent_weeks"]
    guard.verify()
    (runs[1] / "paired.json").write_text("changed")
    with pytest.raises(ValueError, match="Input changed"):
        guard.verify()
    with pytest.raises(ValueError, match="Recorded input hash mismatch"):
        load_weeks([base])


def test_overlapping_weeks_cannot_establish_repeatability(saved_weeks):
    from scripts.combine_microseism_weeks import load_weeks, summarize

    base, _ = saved_weeks
    stations, _ = load_weeks([base])
    stations["XX.S"]["weeks"][1]["interval_utc"] = [
        "2020-01-02T00:00:00Z",
        "2020-01-09T00:00:00Z",
    ]
    result = summarize(stations)[0]["frequency_consistency"]
    assert not result["conditional_reportable"]
    assert result["status"] == "overlapping_or_invalid_intervals"


def test_wrong_network_shift_is_excluded_and_recorded(saved_weeks):
    from scripts.combine_microseism_weeks import week_runs

    base, runs = saved_weeks
    path = runs[1] / "settings.json"
    data = json.loads(path.read_text())
    data["networks"] = ["OTHER"]
    write_json(path, data)
    excluded = []
    selected = week_runs(base, excluded=excluded)
    assert len(selected) == 2
    assert excluded == [
        {"run": str(runs[1]), "reason": "networks_or_duration_differ_from_base"}
    ]


def test_different_earth_models_are_not_pooled(saved_weeks):
    from scripts.combine_microseism_weeks import load_weeks

    base, runs = saved_weeks
    path = runs[1] / "check.json"
    data = json.loads(path.read_text())
    data["model"]["layers_sha256"] = "other-model"
    write_json(path, data)
    with pytest.raises(ValueError, match="epochs vary"):
        load_weeks([base])


@pytest.mark.parametrize("arguments", [["--network", "YV"], ["--days", "14"]])
def test_shift_runner_refuses_different_network_or_duration(
    saved_weeks, monkeypatch, arguments
):
    from scripts.run_microseism_check import main

    base, _ = saved_weeks
    monkeypatch.setattr(
        "sys.argv",
        ["run_microseism_check.py", "--base-run", str(base), "--shift-weeks", "1"]
        + arguments,
    )
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 2


def test_review_outputs_are_fresh_and_preserve_legacy_results(saved_weeks, monkeypatch):
    from scripts.combine_microseism_weeks import main

    base, _ = saved_weeks
    sentinel = write_json(
        base.parent / "weeks_combined/microseism_weeks.json", {"old": True}
    )
    out = base.parent / "new_review"
    monkeypatch.setattr(
        "sys.argv", ["combine_microseism_weeks.py", str(base), "--output-dir", str(out)]
    )
    main()
    result = json.loads((out / "microseism_weeks.json").read_text())
    assert result["frequency_status_counts"] == {"stable": 1}
    assert result["protected_inputs_unchanged"]
    assert result["input_sha256"]
    assert json.loads(sentinel.read_text()) == {"old": True}
    assert (out / "microseism_frequency_table.md").is_file()
    with pytest.raises(FileExistsError):
        main()
