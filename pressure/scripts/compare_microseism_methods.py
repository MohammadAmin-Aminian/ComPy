"""Four-column microseism comparison for a completed all-network run.

1. Claude microseism check (this branch): phase screen, implied pressure sign,
   model-conditional gain bracket where informative.
2. Rayleigh teleseismic calibration: saved baseline and pole-trial gains.
3. Codex microseism diagnostic: ``compy_pressure.microseism_diagnostic.run_diagnostic``
   run here on the *same* paired spectra as column 1 (coherence, phase,
   transfer amplitude; Codex's method fits no gain).
4. Published or reported values, converted to the g direction.

Reads the run folder and tracked metadata; writes only inside the run folder
(and curated copies with --docs). No gain is fitted, applied or stored.
"""

import argparse
import csv
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from obspy import read_inventory
from obspy.core.inventory.response import PolesZerosResponseStage

from compy_pressure.microseism_check import sha256
from compy_pressure.microseism_compare import phase_family
from compy_pressure.microseism_diagnostic import run_diagnostic
from compy_pressure.microseism_validation import RULE_ID
from compy_pressure.ppsd import atomic_json
from scripts.summarize_yv_microseism import (
    EXPLAINED,
    INFORMATIVE_CONDITIONING,
    load_bins,
    phase_screen,
)

ROOT = Path(__file__).resolve().parents[1]
CODEX_CAMPAIGN = (
    Path.home() / "Projects/bruit-fm/docs/microseism_baseline_comparison.csv"
)
SIGN_BY_PHASE = {"nominal": -1, "reflected_pole": -1}
SIGN_BY_PHASE |= {f"reversed_{k}": 1 for k in SIGN_BY_PHASE}
# Deng et al. (2022) XO values, Cp multiplies recorded pressure (g direction).
# Values are per band within their ~0.02-0.07 Hz analysis; a comparison is
# frequency matched only when the check's sub-band overlaps exactly one band.
DENG_2022 = {
    "WD48": (
        [((0.02, 0.07), 0.95)],
        "Deng 2022 text: Cp ~0.95, nearly frequency independent (0.02-0.07 Hz)",
    ),
    "WD70": (
        [((0.02, 0.07), (0.91, 1.71))],
        "Deng 2022 text: Cp 0.91-1.71, frequency dependent (band mapping not given)",
    ),
    "WD66": (
        [((0.02, 0.04), 0.92), ((0.06, 0.07), 0.80)],
        "Deng 2022 Fig. 4g: ~0.92 (0.02-0.04 Hz) to ~0.80 (0.06-0.07 Hz)",
    ),
}
NO_PUBLISHED = {
    "7D": "Ruan et al. 2014 fit DPG scale factors but did not tabulate them",
    "XO": "Deng et al. 2022 values for this station are only in its supplement",
    "YV": "not in Aminian et al. 2025 Table 1",
    "YL": "not in Zha & Webb 2016 Table 1",
}


def positive_real_poles(snapshot, seed_id):
    """Positive-real Laplace poles (rad/s) of one channel in a response snapshot."""
    net, sta, loc, cha = seed_id.split(".")
    selected = read_inventory(str(snapshot)).select(
        network=net, station=sta, location=loc, channel=cha
    )
    poles = []
    for channel in (c for n in selected for s in n for c in s):
        for stage in channel.response.response_stages:
            if not isinstance(stage, PolesZerosResponseStage):
                continue
            scale = 2 * np.pi if "HERTZ" in stage.pz_transfer_function_type else 1.0
            for pole in stage.poles:
                pole = complex(pole)
                if pole.real > 0 and abs(pole.imag) <= 1e-9 * abs(pole):
                    poles.append(pole.real * scale)
    return poles


def channel_poles(paired):
    """Pressure and vertical positive-real poles from the paired run's snapshots."""
    result = {}
    for key, seed in (
        ("pressure", paired["pressure_epoch"]["Seed_id"]),
        ("vertical", paired["vertical_seed_id"]),
    ):
        found = [
            positive_real_poles(path, seed) for path in paired["response_snapshots"]
        ]
        result[key] = sorted({p for group in found for p in group})
    return result


def codex_column(record, paired_path, multiplier, output):
    """Run Codex's diagnostic on the same paired spectra; never fits a gain.

    run_diagnostic refuses an existing folder, so a rerun reuses the saved
    diagnostic only if it was made from the same spectra and multiplier.
    """
    path = Path(output) / "diagnostic.json"
    if path.exists():
        diagnostic = json.loads(path.read_text())
        same = (
            Path(diagnostic["source_comparison"]) == Path(paired_path).resolve()
            and diagnostic["results"]["fixed_pressure_multiplier"] == multiplier
        )
        if not same:
            raise ValueError(f"{output} holds a diagnostic for other inputs")
        artifacts = diagnostic.get("provenance", {}).get("saved_input_artifacts", {})
        if not {"comparison", "spectra"}.issubset(artifacts):
            raise ValueError(f"{output} lacks saved spectral input hashes")
        if (
            Path(artifacts["comparison"]["path"]).resolve()
            != Path(paired_path).resolve()
        ):
            raise ValueError(f"{output} records a different paired comparison")
        for artifact in artifacts.values():
            if sha256(artifact["path"]) != artifact["sha256"]:
                raise ValueError(f"{output} has a stale diagnostic: input hash changed")
    else:
        path = Path(run_diagnostic(paired_path, output, pressure_multiplier=multiplier))
        diagnostic = json.loads(path.read_text())
    results = diagnostic["results"]
    saikia = results["saikia_daily_reference"]
    phase = results["median_nominal_transfer_phase_deg"]
    return {
        "codex_diagnostic": str(path),
        "codex_display_multiplier": multiplier,
        "codex_median_coherence_squared": results["median_coherence_squared"],
        "codex_median_phase_deg": phase,
        "codex_phase_family": phase_family(phase)[0],
        "codex_ruan_high_coherence_fraction": results["ruan_high_coherence_reference"][
            "fraction_of_band"
        ],
        "codex_selected_days": f"{saikia['selected_days']}/{saikia['available_days']}",
        "codex_nominal_amplitude_m_pa": results[
            "median_nominal_transfer_amplitude_m_pa"
        ],
        "codex_gain_estimate": "none (method reports no gain)",
    }


def codex_campaign():
    """The user's completed 46-station Codex campaign, if available."""
    if not CODEX_CAMPAIGN.exists():
        return {}
    with CODEX_CAMPAIGN.open() as stream:
        return {r["station"]: r for r in csv.DictReader(stream)}


def published(record):
    """(low, high, source, frequency-matched value or None).

    Table values (Aminian 2025, Zha & Webb 2016) are single teleseismic gains and
    count as matched. Deng values are matched only when the check's whole
    sub-band lies within exactly one of their bands with a single value.
    """
    network, station = record["network"], record["station"]
    if record.get("published_g") is not None:
        source = (
            "Aminian et al. 2025 Table 1 (same Rayleigh P/a method family as column 2)"
            if network == "YV"
            else ("Zha & Webb 2016 Table 1 (1/g inverted)")
        )
        value = record["published_g"]
        return value, value, source, value
    if network == "XO" and station in DENG_2022:
        bands, source = DENG_2022[station]
        values = [v for _, v in bands for v in (v if isinstance(v, tuple) else (v,))]
        sub = record.get("best_sub_band_hz")
        hits = [v for (lo, hi), v in bands if sub and lo <= sub[0] < sub[1] <= hi]
        matched = hits[0] if len(hits) == 1 and not isinstance(hits[0], tuple) else None
        return min(values), max(values), source, matched
    return None, None, NO_PUBLISHED[network], None


def percent(value, reference):
    if value is None or reference is None:
        return None
    return 100 * (value / reference - 1)


def compare(run, docs=None, codex=True, show=False):
    run = Path(run)
    stations = json.loads((run / "summary.json").read_text())["stations"]
    campaign = codex_campaign()
    rows = []
    for record in stations:
        key = f"{record['network']}.{record['station']}"
        row = {"station": key, "network": record["network"], "status": record["status"]}
        low, high, source, matched = published(record)
        row.update(
            published_g_low=low,
            published_g_high=high,
            published_g_frequency_matched=matched,
            published_source=source,
            rayleigh_baseline_g=record.get("baseline_g"),
            rayleigh_trial_g=record.get("trial_g"),
            rayleigh_assumed_sign=record.get("rayleigh_pressure_sign"),
        )
        if low is not None and low == high and record.get("baseline_g") is not None:
            row["rayleigh_baseline_vs_published_pct"] = percent(
                record["baseline_g"], low
            )
        if record["status"] != "checked":
            row["reason"] = record.get("paired_reason") or record.get("error")
            rows.append(row)
            continue
        paired_path = Path(record["paired_result"])
        paired = json.loads(paired_path.read_text())
        poles = channel_poles(paired)
        check = json.loads(Path(record["check_summary"]).read_text())
        bins = load_bins(check["files"]["bins"])
        screen = phase_screen(
            bins["frequency_hz"],
            bins["phase_deg"],
            bins["coherence_squared"],
            poles["pressure"],
            record["first_resonance_hz"],
            poles["vertical"],
        )
        sign = SIGN_BY_PHASE.get(screen["label"].replace("approximate_", ""))
        conditioning = record["best_sub_band_conditioning_ratio"]
        informative = bool(
            conditioning <= INFORMATIVE_CONDITIONING and screen["label"] in EXPLAINED
        )
        g1, g2 = (
            record["best_sub_band_mode0_water_h1"],
            record["best_sub_band_mode0_water_h2"],
        )
        row.update(
            depth_m=record["depth_m"],
            interval_start=record["start_utc"][:10],
            retained_windows=record["retained_windows"],
            pressure_positive_poles_rad_s=";".join(f"{p:g}" for p in poles["pressure"]),
            vertical_positive_poles_rad_s=";".join(f"{p:g}" for p in poles["vertical"]),
            claude_phase_label=screen["label"],
            claude_phase_mean_deg=screen.get("circular_mean_phase_deg"),
            claude_phase_rms_deg=(
                min(screen["rms_residual_deg"].values())
                if "rms_residual_deg" in screen
                else None
            ),
            claude_implied_sign=sign,
            sign_consistent_with_rayleigh=(
                None if sign is None else sign == record.get("rayleigh_pressure_sign")
            ),
            claude_conditioning=conditioning,
            claude_gain_informative=informative,
            claude_gain_status="legacy_single_week_candidate",
            claude_sub_band_hz="-".join(f"{v:.4f}" for v in record["best_sub_band_hz"]),
            claude_g_h1=g1 if informative else None,
            claude_g_h2=g2 if informative else None,
            # Structure-free limit: vertical incidence with H2; mode-0 values
            # above depend on the Earth model up to this ratio.
            claude_structure_free_max=(
                record["best_sub_band_upper_bound_h2"] if informative else None
            ),
        )
        if informative:
            for name, reference in (
                ("baseline", record.get("baseline_g")),
                ("trial", record.get("trial_g")),
                ("published", matched),
            ):
                if reference is not None:
                    row[f"claude_vs_{name}_pct"] = (
                        f"{percent(g1, reference):+.1f}..{percent(g2, reference):+.1f}"
                    )
                    row[f"{name}_in_claude_bracket"] = bool(g1 <= reference <= g2)
        if codex:
            multiplier = record.get("baseline_g") or record.get("trial_g") or 1.0
            output = run / "codex_diagnostic" / record["network"] / record["station"]
            try:
                row.update(codex_column(record, paired_path, multiplier, output))
            except Exception as exc:  # one station must not abort the comparison
                row["codex_error"] = repr(exc)
        prior = campaign.get(key)
        if prior and Path(prior["source_comparison"]).exists():
            prior_paired = json.loads(Path(prior["source_comparison"]).read_text())
            saved_phase = prior.get("median_nominal_transfer_phase_deg")
            row["codex_campaign_phase_deg"] = (
                None if saved_phase in (None, "") else float(saved_phase)
            )
            row["codex_campaign_same_interval"] = prior_paired.get(
                "start_utc"
            ) == paired.get("start_utc")
        rows.append(row)

    fields = sorted({k for r in rows for k in r}, key=lambda k: (k != "station", k))
    with (run / "method_comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    counts = {}
    for r in rows:
        label = r.get("claude_phase_label", r["status"])
        counts[f"{r['network']}:{label}"] = counts.get(f"{r['network']}:{label}", 0) + 1
    summary = {
        "run": str(run),
        "run_summary_sha256": sha256(run / "summary.json"),
        "codex_campaign_csv": str(CODEX_CAMPAIGN) if campaign else None,
        "columns": {
            "1": "Claude microseism check (phase screen; gain bracket if informative)",
            "2": "Rayleigh teleseismic calibration, saved baseline/trial g",
            "3": "Codex microseism diagnostic on the same paired spectra (no gain)",
            "4": "Published/reported values in the g direction",
        },
        "phase_label_counts": counts,
        "informative_stations": [
            r["station"] for r in rows if r.get("claude_gain_informative")
        ],
        "status": "diagnostic_only; no gain fitted or applied",
        "gain_policy": f"Legacy best-band candidates; use {RULE_ID} for week review",
    }
    atomic_json(run / "method_comparison.json", summary | {"rows": rows})
    write_markdown(rows, summary, run / "method_comparison.md")
    figure = plot(rows, run / "method_comparison.png")
    if show:
        plt.show()
    plt.close(figure)
    if docs:
        docs = Path(docs)
        for name, target in (
            ("method_comparison.csv", docs / "microseism_method_comparison.csv"),
            ("method_comparison.json", docs / "microseism_method_comparison.json"),
            ("method_comparison.md", docs / "microseism_method_comparison_table.md"),
            (
                "method_comparison.png",
                docs / "figures/microseism_method_comparison.png",
            ),
        ):
            shutil.copyfile(run / name, target)
    return rows, summary


def combine(runs, docs):
    """Merge saved method_comparison.json rows from several runs into docs/."""
    rows, sources, counts, informative = [], [], {}, []
    for run in map(Path, runs):
        data = json.loads((run / "method_comparison.json").read_text())
        rows += data["rows"]
        sources.append(
            {"run": str(run), "run_summary_sha256": data["run_summary_sha256"]}
        )
        for label, count in data["phase_label_counts"].items():
            counts[label] = counts.get(label, 0) + count
        informative += data["informative_stations"]
    keys = [r["station"] for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("--combine got the same station from more than one run")
    summary = {
        "runs": sources,
        "columns": data["columns"],
        "phase_label_counts": counts,
        "informative_stations": informative,
        "status": "diagnostic_only; no gain fitted or applied",
    }
    docs = Path(docs)
    fields = sorted({k for r in rows for k in r}, key=lambda k: (k != "station", k))
    with (docs / "microseism_method_comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    atomic_json(docs / "microseism_method_comparison.json", summary | {"rows": rows})
    write_markdown(rows, summary, docs / "microseism_method_comparison_table.md")
    plt.close(plot(rows, docs / "figures/microseism_method_comparison.png"))
    return rows, summary


def fmt(value, digits=3):
    if value is None or value == "":
        return "—"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def write_markdown(rows, summary, destination):
    lines = [
        "# Microseism method comparison (generated)",
        "",
        "Columns: (1) Claude microseism check, (2) Rayleigh teleseismic g,",
        "(3) Codex microseism diagnostic on the same spectra, (4) published g.",
        f"These are legacy single-week candidates, superseded by `{RULE_ID}`.",
        "Gain brackets appear only where column 1 is informative. The H1–H2",
        "bracket is a one-week noise bracket for the mode-0 RR52-model value;",
        "the structure-free max is the vertical-incidence limit (H2), so the",
        "Earth model can move g anywhere up to it. † = published value not",
        "frequency matched to the check's sub-band (not a test).",
        "",
        "| Station | (1) phase label | (1) sign | (1) g H1–H2 | (1) structure-free"
        " max | (2) baseline g | (2) trial g | (3) Codex phase | (3) coh² 0.1–0.2"
        " | (4) published g | (4) source |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        bracket = (
            f"{r['claude_g_h1']:.3f}–{r['claude_g_h2']:.3f}"
            if r.get("claude_g_h1") is not None
            else "—"
        )
        pub = r.get("published_g_low")
        pub_text = (
            "—"
            if pub is None
            else (
                f"{pub:.3f}"
                if pub == r["published_g_high"]
                else f"{pub:.2f}–{r['published_g_high']:.2f}"
            )
        )
        if pub is not None and r.get("published_g_frequency_matched") is None:
            pub_text += " †"
        sign = r.get("claude_implied_sign")
        sign_text = (
            "—"
            if sign is None
            else (
                f"{sign:+d}"
                + ("" if r.get("sign_consistent_with_rayleigh") else " ≠ Rayleigh")
            )
        )
        label = r.get("claude_phase_label", r["status"])
        cells = [
            r["station"],
            label,
            sign_text,
            bracket,
            fmt(r.get("claude_structure_free_max")),
            fmt(r.get("rayleigh_baseline_g")),
            fmt(r.get("rayleigh_trial_g")),
            fmt(r.get("codex_median_phase_deg"), 1),
            fmt(r.get("codex_median_coherence_squared"), 2),
            pub_text,
            r.get("published_source", ""),
        ]
        lines.append("| " + " | ".join(cells) + " |")
    lines += ["", f"Phase labels: `{json.dumps(summary['phase_label_counts'])}`", ""]
    destination.write_text("\n".join(lines))


def plot(rows, destination):
    done = [r for r in rows if r["status"] == "checked"]
    x = np.arange(len(done))
    figure, axes = plt.subplots(
        2, 1, figsize=(max(15, len(done) * 0.22), 10), sharex=True
    )
    for i, r in enumerate(done):
        if r.get("claude_phase_mean_deg") is not None:
            axes[0].plot(i, r["claude_phase_mean_deg"], "bo", ms=4)
        if r.get("codex_median_phase_deg") is not None:
            axes[0].plot(i, r["codex_median_phase_deg"], "rx", ms=5)
    axes[0].plot([], [], "bo", label="(1) Claude coherent-bin circular mean")
    axes[0].plot([], [], "rx", label="(3) Codex median 0.1–0.2 Hz")
    for value in (0, 180, -180):
        axes[0].axhline(value, color="gray", lw=0.8)
    axes[0].set(ylabel="P/Z phase (deg)", ylim=(-200, 200), title="Phase: same spectra")
    axes[0].legend(fontsize=8)
    for i, r in enumerate(done):
        if r.get("claude_g_h1") is not None:
            axes[1].vlines(i, r["claude_g_h1"], r["claude_g_h2"], lw=5, alpha=0.5)
        for key, style in (("rayleigh_baseline_g", "ko"), ("rayleigh_trial_g", "C1x")):
            if r.get(key) is not None:
                axes[1].plot(i, r[key], style, ms=4)
        if r.get("published_g_frequency_matched") is not None:
            axes[1].plot(i, r["published_g_frequency_matched"], "g^", ms=4)
    axes[1].plot(
        [], [], lw=5, alpha=0.5, label="(1) Claude microseism H1–H2 (informative)"
    )
    axes[1].plot([], [], "ko", label="(2) Rayleigh baseline")
    axes[1].plot([], [], "C1x", label="(2) Rayleigh trial")
    axes[1].plot([], [], "g^", label="(4) published")
    axes[1].set(ylabel="pressure multiplier g", ylim=(0, 2), xticks=x)
    axes[1].set_xticklabels([r["station"] for r in done], rotation=90, fontsize=6)
    axes[1].legend(fontsize=8, ncol=4)
    for axis in axes:
        axis.grid(True, alpha=0.3)
    figure.tight_layout()
    figure.savefig(destination, dpi=130)
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "run", type=Path, nargs="+", help="Completed run_microseism_check.py folder"
    )
    parser.add_argument("--docs", type=Path, help="Copy curated outputs into docs/")
    parser.add_argument("--no-codex", action="store_true", help="Skip column 3 runs")
    parser.add_argument(
        "--combine",
        action="store_true",
        help="Merge already compared runs into --docs without recomputing",
    )
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    if args.combine:
        if not args.docs:
            parser.error("--combine needs --docs")
        _, summary = combine(args.run, args.docs)
    elif len(args.run) == 1:
        _, summary = compare(args.run[0], args.docs, not args.no_codex, args.show)
    else:
        parser.error("compare one run at a time, or merge several with --combine")
    print(json.dumps(summary["phase_label_counts"]))
    print(f"informative: {summary['informative_stations']}")


if __name__ == "__main__":
    main()
