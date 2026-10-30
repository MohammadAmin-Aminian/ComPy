"""Review saved microseism weeks with legacy and frozen frequency rules.

For each base run, finds the runs made from it with --base-run/--shift-weeks in
the same output root and reads their method_comparison.json; nothing is
recomputed. The old best-band calculation is retained for comparison. The new
rule requires consistent sub-bands <=0.10 Hz and >=3 passing independent weeks.
Neither repeatability nor H1-H2 spread establishes absolute calibration or a
confidence interval. All derived outputs use a fresh directory.
"""

import argparse
import csv
import hashlib
import json
import shutil
import statistics
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np

from compy_pressure.microseism_validation import RULE, frequency_station, frequency_week
from compy_pressure.ppsd import atomic_json

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

STABLE_PCT, MODERATE_PCT = 5.0, 10.0  # fixed before the repeated weeks were run
REFERENCES = (
    "rayleigh_baseline_g",
    "rayleigh_trial_g",
    "published_g_low",
    "published_g_high",
    "published_g_frequency_matched",
    "published_source",
)


class InputGuard:
    """Hash exact saved inputs, checking recorded hashes and concurrent edits."""

    def __init__(self):
        self.hashes = {}

    @staticmethod
    def digest(path):
        digest = hashlib.sha256()
        with Path(path).open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    def protect(self, path, expected=None):
        path = Path(path).resolve()
        if str(path) not in self.hashes:
            self.hashes[str(path)] = self.digest(path)
        if expected is not None and self.hashes[str(path)] != expected:
            raise ValueError(f"Recorded input hash mismatch: {path}")
        return path

    def json(self, path):
        return json.loads(self.protect(path).read_text())

    def verify(self):
        for path, expected in self.hashes.items():
            if self.digest(path) != expected:
                raise ValueError(f"Input changed during review: {path}")


def week_runs(base, skipped=None, excluded=None):
    """[(shift_weeks, run)] for the base run (shift 0) and runs made from it.

    Identical repeated shifts keep their latest completed run (folder names
    sort by time); the others are recorded in ``skipped``. Conflicting repeats
    require review rather than assuming the newest is scientifically equivalent.
    """
    base = Path(base).resolve()
    original = json.loads((base / "settings.json").read_text())
    latest = {}
    for settings in sorted(base.parent.glob("run_*/settings.json")):
        saved = json.loads(settings.read_text())
        run = settings.parent
        if not saved.get("base_run") or Path(saved["base_run"]).resolve() != base:
            continue
        if not (run / "summary.json").exists():
            continue
        reason = None
        if not (run / "method_comparison.json").exists():
            reason = "comparison_not_complete"
        elif any(
            key in original and saved.get(key) != original[key]
            for key in ("networks", "days")
        ):
            reason = "networks_or_duration_differ_from_base"
        if reason:
            if excluded is not None:
                excluded.append({"run": str(run), "reason": reason})
            continue
        shift = saved["shift_weeks"]
        if shift in latest:

            def scientific_rows(folder):
                rows = json.loads((folder / "method_comparison.json").read_text())
                return sorted(
                    json.dumps(
                        {k: v for k, v in row.items() if k != "codex_diagnostic"},
                        sort_keys=True,
                    )
                    for row in rows.get("rows", [])
                )

            if scientific_rows(latest[shift]) != scientific_rows(run):
                raise ValueError(
                    f"Conflicting completed repeats for shift {shift}: "
                    f"review {latest[shift]} and {run}"
                )
            if skipped is not None:
                skipped.append(str(latest[shift]))
        latest[shift] = run
    return sorted([(0, base), *latest.items()])


def aggregate(weeks):
    """Summarise one station's weeks (dicts: shift, start, label, sign, g1, g2)."""
    checked = [w for w in weeks if w["label"] is not None]
    good = [w for w in checked if w["g1"] is not None]
    mids = [(w["g1"] + w["g2"]) / 2 for w in good]
    signs = {w["sign"] for w in checked if w["sign"] is not None}
    out = {
        "weeks_checked": len(checked),
        "weeks_informative": len(good),
        "phase_labels": sorted({w["label"] for w in checked}),
        "sign_stable": (len(signs) == 1) if signs else None,
    }
    if not mids:
        out["stability"] = "not_informative"
        return out
    median = statistics.median(mids)
    spread = 100 * (max(mids) - min(mids)) / median
    out.update(
        median_g=median,
        min_week_g=min(mids),
        max_week_g=max(mids),
        spread_pct=spread,
        pooled_h1=min(w["g1"] for w in good),
        pooled_h2=max(w["g2"] for w in good),
    )
    if len(good) < 2:
        out["stability"] = "single_week"
    elif spread <= STABLE_PCT:
        out["stability"] = "stable"
    elif spread <= MODERATE_PCT:
        out["stability"] = "moderate"
    else:
        out["stability"] = "variable"
    return out


def load_weeks(bases, skipped=None, excluded=None, guard=None):
    guard = guard or InputGuard()
    stations, runs = {}, []
    seen_bases = set()
    for base in bases:
        base = Path(base).resolve()
        if base in seen_bases:
            raise ValueError("A base run was supplied more than once")
        seen_bases.add(base)
        for shift, run in week_runs(base, skipped, excluded):
            path = run / "method_comparison.json"
            if not path.exists():
                raise FileNotFoundError(f"Run compare_microseism_methods.py on {run}")
            guard.json(run / "settings.json")
            comparison = guard.json(path)
            guard.protect(run / "summary.json", comparison["run_summary_sha256"])
            records = {
                f"{r['network']}.{r['station']}": r
                for r in guard.json(run / "summary.json")["stations"]
            }
            seen = set()
            runs.append({"run": str(run), "shift_weeks": shift})
            for r in comparison["rows"]:
                key = r["station"]
                if key in seen:
                    raise ValueError(f"Duplicate station in {path}: {key}")
                seen.add(key)
                entry = stations.setdefault(key, {"weeks": [], "base_run": str(base)})
                if entry["base_run"] != str(base):
                    raise ValueError(f"Station occurs in multiple base runs: {key}")
                if shift == 0:
                    entry["base"] = r
                week = {
                    "shift": shift,
                    "start": r.get("interval_start"),
                    "label": r.get("claude_phase_label"),
                    "sign": r.get("claude_implied_sign"),
                    "g1": r.get("claude_g_h1"),
                    "g2": r.get("claude_g_h2"),
                    "run": str(run),
                }
                if r["status"] == "checked":
                    record = records[key]
                    check = guard.json(record["check_summary"])
                    if check["station"] != key:
                        raise ValueError(
                            f"Check summary belongs to another station: {key}"
                        )
                    for name, metadata in check["inputs"].items():
                        guard.protect(name, metadata["sha256"])
                    paired = guard.json(record["paired_result"])
                    identity = {
                        "water": check["water"],
                        "layers_sha256": check["model"]["layers_sha256"],
                        "pressure_epoch": paired["pressure_epoch"],
                        "vertical_epoch": paired["vertical_epoch"],
                    }
                    if entry.setdefault("science_inputs", identity) != identity:
                        raise ValueError(
                            f"Water/model/channel epochs vary across {key}"
                        )
                    week["interval_utc"] = check["interval_utc"]
                    week["frequency_check"] = frequency_week(check, week["label"])
                    week["codex_amplitude_m_pa"] = r.get("codex_nominal_amplitude_m_pa")
                entry["weeks"].append(week)
    return stations, runs


def percent(value, reference):
    return None if value is None or not reference else 100 * (value / reference - 1)


def summarize(stations):
    rows = []
    for key, entry in stations.items():
        base = entry.get("base", {})
        row = {"station": key, "network": key.split(".")[0]}
        row |= aggregate(sorted(entry["weeks"], key=lambda w: w["shift"]))
        row |= {k: base.get(k) for k in REFERENCES}
        median = row.get("median_g")
        for name, reference in (
            ("baseline", base.get("rayleigh_baseline_g")),
            ("trial", base.get("rayleigh_trial_g")),
            ("published", base.get("published_g_frequency_matched")),
        ):
            row[f"median_vs_{name}_pct"] = percent(median, reference)
            if median is not None and reference:
                row[f"{name}_within_week_range"] = bool(
                    row["min_week_g"] <= reference <= row["max_week_g"]
                )
        row["weeks"] = sorted(entry["weeks"], key=lambda w: w["shift"])
        row["frequency_consistency"] = frequency_station(row["weeks"])
        passing = [
            w for w in row["weeks"] if w.get("frequency_check", {}).get("passes")
        ]
        intervals = sorted(
            [
                datetime.fromisoformat(v.replace("Z", "+00:00"))
                for v in w["interval_utc"]
            ]
            for w in passing
        )
        independent = all(start < end for start, end in intervals) and all(
            left[1] <= right[0] for left, right in zip(intervals, intervals[1:])
        )
        row["frequency_consistency"]["independent_weeks"] = independent
        if not independent:
            row["frequency_consistency"].update(
                conditional_reportable=False, status="overlapping_or_invalid_intervals"
            )
        rows.append(row)
    return rows


def write_frequency_markdown(rows, destination):
    lines = [
        "# Conditional microseism estimates: frequency and week consistency",
        "",
        f"Rule `{RULE['id']}` was designed on YL; all values assume the vertical",
        "response and historical RR52 model. No calibration is accepted/applied.",
        "At least two eligible bands <=0.10 Hz must agree within 10%, and at least",
        "three independent passing weeks must agree within 10% (stable <=5%).",
        "Legacy best-band values remain in microseism_weeks_table.md.",
        "H1-H2 envelopes and temporal spread are descriptive, not confidence limits.",
        "",
        "| Station | Passing / checked weeks | Status | Conditional g | Week spread"
        " | Rayleigh baseline | Rayleigh trial | Published scalar g |",
        "|---|---:|---|---:|---:|---:|---:|---:|",
    ]
    for row in sorted(rows, key=lambda r: r["station"]):
        review = row["frequency_consistency"]
        value = review["median_g"] if review["conditional_reportable"] else None
        # Deng's frequency-dependent values are not matched to this multi-band rule.
        published = (
            row.get("published_g_low")
            if row["network"] in {"YV", "YL"}
            and row.get("published_g_low") == row.get("published_g_high")
            else None
        )
        cells = [
            row["station"],
            f"{review['passing_weeks']}/{row['weeks_checked']}",
            review["status"],
            fmt(value),
            fmt(review["week_spread_pct"], 1),
            fmt(row.get("rayleigh_baseline_g")),
            fmt(row.get("rayleigh_trial_g")),
            fmt(published),
        ]
        lines.append("| " + " | ".join(cells) + " |")
    destination.write_text("\n".join(lines) + "\n")


def fmt(value, digits=3):
    return "—" if value is None else f"{value:.{digits}f}"


def write_markdown(rows, destination):
    lines = [
        "# Microseism check: week-to-week spread (generated)",
        "",
        "Legacy best-band diagnostic; use microseism_frequency_table.md for the",
        "frozen frequency/week rule. These values are not accepted calibrations.",
        "",
        "Median and range of informative weekly mode-0 midpoints. Stability:",
        f"stable ≤ {STABLE_PCT:g}%, moderate ≤ {MODERATE_PCT:g}%, otherwise",
        "variable (range / median). Published-band matches are inherited from",
        "the historical base-week comparison and may not match every later week.",
        "",
        "| Station | weeks (inf./checked) | median g | week range | spread"
        " | stability | Rayleigh baseline | vs baseline | published | vs published"
        " | phase labels |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        if r["weeks_checked"] == 0:
            continue
        week_range = (
            f"{r['min_week_g']:.3f}–{r['max_week_g']:.3f}"
            if r.get("median_g") is not None
            else "—"
        )
        spread = r.get("spread_pct")
        cells = [
            r["station"],
            f"{r['weeks_informative']}/{r['weeks_checked']}",
            fmt(r.get("median_g")),
            week_range,
            "—" if spread is None else f"{spread:.1f}%",
            r["stability"],
            fmt(r.get("rayleigh_baseline_g")),
            fmt(r.get("median_vs_baseline_pct"), 1),
            fmt(r.get("published_g_frequency_matched")),
            fmt(r.get("median_vs_published_pct"), 1),
            ", ".join(r["phase_labels"]),
        ]
        lines.append("| " + " | ".join(cells) + " |")
    destination.write_text("\n".join(lines) + "\n")


def plot(rows, destination):
    shown = [r for r in rows if r.get("median_g") is not None]
    figure, axis = plt.subplots(figsize=(max(12, 0.3 * len(shown)), 6))
    for i, r in enumerate(shown):
        good = [w for w in r["weeks"] if w["g1"] is not None]
        offsets = np.linspace(-0.25, 0.25, len(good)) if len(good) > 1 else [0.0]
        for dx, w in zip(offsets, good, strict=True):
            axis.vlines(i + dx, w["g1"], w["g2"], lw=3, alpha=0.6)
        axis.plot(i, r["median_g"], "b_", ms=12)
        if r.get("rayleigh_baseline_g") is not None:
            axis.plot(i, r["rayleigh_baseline_g"], "ko", ms=4)
        if r.get("published_g_frequency_matched") is not None:
            axis.plot(i, r["published_g_frequency_matched"], "g^", ms=4)
    axis.plot([], [], lw=3, alpha=0.6, label="weekly H1–H2 (informative)")
    axis.plot([], [], "b_", ms=12, label="median of weekly midpoints")
    axis.plot([], [], "ko", label="Rayleigh baseline")
    axis.plot([], [], "g^", label="published (historical base-week reference)")
    axis.set(
        ylabel="pressure multiplier g",
        xticks=range(len(shown)),
        title="Legacy best-band estimates: superseded by frequency/week review",
    )
    axis.set_xticklabels([r["station"] for r in shown], rotation=90, fontsize=7)
    axis.grid(True, alpha=0.3)
    axis.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(destination, dpi=130)
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base", type=Path, nargs="+", help="Base (shift 0) runs")
    parser.add_argument("--docs", type=Path, help="Copy curated outputs into docs/")
    parser.add_argument(
        "--output-dir", type=Path, help="Fresh derived-review directory"
    )
    args = parser.parse_args()
    guard = InputGuard()
    guard.protect(__file__)
    guard.protect(
        Path(__file__).parents[1] / "../compy_pressure/microseism_validation.py"
    )
    skipped, excluded = [], []
    stations, runs = load_weeks(args.base, skipped, excluded, guard)
    rows = summarize(stations)
    out = args.output_dir or (
        Path(__file__).resolve().parents[1]
        / "results/REVIEWS/MICROSEISM_INTEGRATION"
        / f"run_{datetime.now(UTC):%Y%m%dT%H%M%SZ}_{uuid4().hex[:8]}"
    )
    out.mkdir(parents=True, exist_ok=False)
    flat = [
        {k: v for k, v in r.items() if k not in {"weeks", "frequency_consistency"}}
        | {f"frequency_{k}": v for k, v in r["frequency_consistency"].items()}
        for r in rows
    ]
    fields = sorted({k for r in flat for k in r}, key=lambda k: (k != "station", k))
    with (out / "microseism_weeks.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(flat)
    counts = {}
    for r in rows:
        counts[r["stability"]] = counts.get(r["stability"], 0) + 1
    frequency_counts = {}
    for row in rows:
        label = row["frequency_consistency"]["status"]
        frequency_counts[label] = frequency_counts.get(label, 0) + 1
    guard.verify()
    atomic_json(
        out / "microseism_weeks.json",
        {
            "runs": runs,
            "duplicate_runs_skipped": skipped,
            "excluded_runs": excluded,
            "input_sha256": guard.hashes,
            "protected_inputs_unchanged": True,
            "frequency_rule": RULE,
            "frequency_status_counts": frequency_counts,
            "thresholds_pct": {"stable": STABLE_PCT, "moderate": MODERATE_PCT},
            "stability_counts": counts,
            "rows": rows,
            "status": "diagnostic_only; no gain fitted or applied",
        },
    )
    write_markdown(rows, out / "microseism_weeks_table.md")
    write_frequency_markdown(rows, out / "microseism_frequency_table.md")
    plt.close(plot(rows, out / "microseism_weeks.png"))
    if args.docs:
        for name, target in (
            ("microseism_weeks.csv", "microseism_weeks.csv"),
            ("microseism_weeks.json", "microseism_weeks.json"),
            ("microseism_weeks_table.md", "microseism_weeks_table.md"),
            ("microseism_frequency_table.md", "microseism_frequency_table.md"),
            ("microseism_weeks.png", "figures/microseism_weeks.png"),
        ):
            shutil.copyfile(out / name, args.docs / target)
    print(
        f"{len(runs)} runs; legacy {counts}; frequency rule {frequency_counts}: {out}"
    )


if __name__ == "__main__":
    main()
