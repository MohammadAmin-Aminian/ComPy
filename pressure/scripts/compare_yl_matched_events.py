"""Read-only matched-event/response review for YL C07W, C08W and C09W."""

import argparse
import csv
import hashlib
import json
import re
import webbrowser
from html import escape
from pathlib import Path
from uuid import uuid4

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime, read_inventory

STATIONS = ("C07W", "C08W", "C09W")
COLORS = ("tab:red", "tab:blue", "tab:green")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workflow",
        type=Path,
        default=Path("results/NETWORKS/YL_ALL_NO_TILT/_workflow/rayleigh.json"),
    )
    parser.add_argument("--event-id", action="append", help="Optional matched event ID")
    parser.add_argument("--min-magnitude", type=float, default=0)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("results/VALIDATION/YL_MATCHED_EVENTS"),
    )
    parser.add_argument(
        "--show", action="store_true", help="Open the saved HTML review"
    )
    args = parser.parse_args()
    if not np.isfinite(args.min_magnitude):
        raise ValueError("Minimum magnitude must be finite")
    workflow = json.loads(args.workflow.read_text())
    if workflow.get("execution_status") != "complete":
        raise ValueError("Workflow must be complete")
    protected = {str(args.workflow): sha(args.workflow)}
    for p in Path("exports").glob("master_inventory.*"):
        protected[str(p)] = sha(p)
    results, events = {}, {}
    for station in STATIONS:
        entries = [
            e
            for e in workflow["entries"].values()
            if e["seed_id"].split(".")[1] == station
        ]
        if len(entries) != 1:
            raise ValueError(f"Select one saved pressure epoch for {station}")
        path = Path(entries[0]["update"]["rayleigh_result_file"])
        protected[str(path)] = sha(path)
        result = json.loads(path.read_text())
        for name, digest in (
            result["source_inventory_sha256"] | result["response_snapshot_sha256"]
        ).items():
            if sha(name) != digest:
                raise ValueError(f"Saved response/inventory changed: {name}")
            protected[name] = digest
        results[station] = result
        events[station] = {
            e["id"]: e for e in result["events"] if e.get("files", {}).get("spectra")
        }
    shared = set.intersection(*(set(events[s]) for s in STATIONS))
    shared = {
        eid
        for eid in shared
        if events[STATIONS[0]][eid]["magnitude"] >= args.min_magnitude
    }
    if args.event_id:
        missing = set(args.event_id) - shared
        if missing:
            raise ValueError(f"No shared saved spectra for requested events: {missing}")
        shared = set(args.event_id)
    if not shared:
        raise ValueError("No matched saved spectral events satisfy the selection")
    if any(not re.fullmatch(r"[A-Za-z0-9_-]+", eid) for eid in shared):
        raise ValueError("Unsafe saved event identifier")
    for eid in shared:
        origins = {events[s][eid]["origin_utc"] for s in STATIONS}
        if len(origins) != 1:
            raise ValueError(f"Event origins disagree: {eid}")
        for station in STATIONS:
            for path in events[station][eid]["files"].values():
                protected[path] = sha(path)
    output = args.output_root.resolve() / f"run_{uuid4().hex[:12]}"
    output.mkdir(parents=True, exist_ok=False)
    metadata, response_curves = [], {}
    band = results[STATIONS[0]]["settings"]["band_hz"]
    for station in STATIONS:
        if results[station]["settings"]["band_hz"] != band:
            raise ValueError("Saved fit bands differ; no implicit rebinning")
        curves = []
        for key in results[station]["epoch_keys"]:
            seed = key["Seed_id"]
            net, sta, loc, code = seed.split(".")
            inventory = read_inventory(key["inventory_file"]).select(
                network=net,
                station=sta,
                location=loc,
                channel=code,
            )
            channels = [
                c
                for n in inventory
                for st in n
                for c in st
                if str(c.start_date or "") == key["channel_start_time"]
                and str(c.end_date or "") == key["channel_end_time"]
            ]
            if len(channels) != 1:
                raise ValueError(f"Ambiguous exact response epoch: {seed}")
            channel = channels[0]
            sensitivity = channel.response.instrument_sensitivity
            kind = "pressure" if code == "LDH" else "vertical"
            f = np.linspace(*band, 161)
            response = channel.response.get_evalresp_response_for_frequencies(
                f,
                output="DEF" if kind == "pressure" else "VEL",
            )
            if not np.all(np.isfinite(response) & (abs(response) > 0)):
                raise ValueError(f"Invalid response in fit band: {seed}")
            curves.append((kind, response))
            metadata.append(
                dict(
                    station=station,
                    seed_id=seed,
                    kind=kind,
                    input_units=sensitivity.input_units,
                    output_units=sensitivity.output_units,
                    sensitivity=sensitivity.value,
                    sensitivity_frequency_hz=sensitivity.frequency,
                    sample_rate_hz=float(channel.sample_rate),
                    dip_deg=channel.dip,
                    stage_gains=[
                        float(s.stage_gain) for s in channel.response.response_stages
                    ],
                    inventory_file=key["inventory_file"],
                )
            )
        if {kind for kind, _ in curves} != {"pressure", "vertical"}:
            raise ValueError("Expected saved no-tilt P/Z epochs only")
        response_curves[station] = dict(curves)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for station, color in zip(STATIONS, COLORS, strict=True):
        for j, kind in enumerate(("pressure", "vertical")):
            response = response_curves[station][kind]
            reference = response_curves["C08W"][kind]
            axes[0, j].plot(f, abs(response), label=station, color=color)
            axes[1, j].plot(
                f,
                np.angle(response / reference, deg=True),
                color=color,
                label=station,
            )
    axes[0, 0].set_ylabel("|Response| [counts/Pa]")
    axes[0, 1].set_ylabel("|Response| [counts/(m/s)]")
    for j, kind in enumerate(("pressure", "vertical")):
        axes[0, j].set_title(f"Nominal {kind} response")
        axes[1, j].set_ylabel("Phase relative to C08W [degrees]")
        axes[1, j].set_xlabel("Frequency [Hz]")
        extent = max(
            1.0,
            max(
                np.max(
                    abs(
                        np.angle(
                            response_curves[s][kind] / response_curves["C08W"][kind],
                            deg=True,
                        )
                    )
                )
                for s in STATIONS
            )
            * 1.05,
        )
        axes[1, j].set_ylim(-extent, extent)
    for ax in axes.flat:
        ax.grid(alpha=0.25)
        ax.legend()
    fig.suptitle("Supplied responses, unchanged — not hardware calibration validation")
    fig.tight_layout()
    fig.savefig(output / "responses.png", dpi=150)
    plt.close(fig)
    rows, cards = [], []
    all_candidate = 0
    for eid in sorted(shared, key=lambda e: events[STATIONS[0]][e]["origin_utc"]):
        fig, axes = plt.subplots(5, 1, figsize=(11, 14), sharex=True)
        candidate = all(events[s][eid]["status"] == "candidate" for s in STATIONS)
        all_candidate += int(candidate)
        links = []
        for station, color in zip(STATIONS, COLORS, strict=True):
            event = events[station][eid]
            with np.load(event["files"]["spectra"], allow_pickle=False) as a:
                freq = a["frequency_hz"]
                ratio = a["psd_ratio_pa"] / a["expected_pa"]
                label = f"{station}: {event['status']}"
                axes[0].semilogy(freq, a["pressure_psd"], color=color, label=label)
                axes[1].semilogy(freq, a["acceleration_psd"], color=color)
                axes[2].plot(freq, ratio, color=color)
                axes[3].plot(freq, a["coherence_squared"], color=color)
                axes[4].plot(freq, a["phase_deg"], color=color)
                valid = a["accepted_bins"].astype(bool)
                axes[2].scatter(freq[valid], ratio[valid], color=color, s=25)
                med = float(np.median(ratio))
            origin = UTCDateTime(event["origin_utc"])
            left, right = UTCDateTime(event["start_utc"]), UTCDateTime(event["end_utc"])
            body_overlap = any(
                event.get(k) and left <= UTCDateTime(event[k]) <= right
                for k in ("predicted_first_p_utc", "predicted_first_s_utc")
            )
            rows.append(
                dict(
                    event_id=eid,
                    station=station,
                    magnitude=event["magnitude"],
                    origin_utc=event["origin_utc"],
                    status=event["status"],
                    all_three_candidates=candidate,
                    saved_pressure_multiplier=event.get("pressure_multiplier"),
                    median_measured_over_expected=med,
                    start_seconds_after_origin=float(left - origin),
                    end_seconds_after_origin=float(right - origin),
                    first_body_arrival_in_window=body_overlap,
                    saved_other_arrival_ids=";".join(
                        event.get("other_predicted_arrival_ids", [])
                    ),
                )
            )
            for kind in ("arrival_review", "waveform", "plot"):
                if path := event["files"].get(kind):
                    url = escape(Path(path).resolve().as_uri(), quote=True)
                    links.append(f'<a href="{url}">{station} {kind}</a>')
        axes[0].legend(fontsize=8)
        for ax, label in zip(
            axes,
            (
                "Pressure PSD [Pa²/Hz]",
                "Upward acceleration PSD [(m/s²)²/Hz]",
                "Measured P/a / saved theory",
                "Squared coherence",
                "P/a phase [degrees]",
            ),
            strict=True,
        ):
            ax.set_ylabel(label)
            ax.grid(alpha=0.25)
        axes[2].axhline(1, ls="--", color="black")
        axes[3].set_ylim(0, 1.02)
        axes[4].axhspan(-20, 20, color="gray", alpha=0.1)
        axes[4].set_ylim(-180, 180)
        axes[-1].set_xlabel("Frequency [Hz]")
        first = events[STATIONS[0]][eid]
        fig.suptitle(
            f"{eid} — M {first['magnitude']} — {first['origin_utc']}\n"
            "Same earthquake; station-specific windows. No gain applied or refitted."
        )
        fig.tight_layout()
        filename = f"{eid}.png"
        fig.savefig(output / filename, dpi=130)
        plt.close(fig)
        cards.append(
            f"<section><h2>{escape(eid)} — all three candidates: {candidate}</h2>"
            f'<p>{" | ".join(links)}</p><img src="{filename}" width="850"></section>'
        )
    changed = [p for p, h in protected.items() if sha(p) != h]
    if changed:
        raise ValueError(f"Protected inputs changed: {changed}")
    write_csv(output / "events.csv", rows)
    write_csv(output / "responses.csv", metadata)
    audit = dict(
        stations=STATIONS,
        workflow=str(args.workflow),
        matched_saved_events=len(shared),
        all_three_candidate_events=all_candidate,
        metadata=metadata,
        protected_sha256=protected,
        protected_inputs_preserved=True,
        gains_applied=False,
        spectra_recomputed=False,
        caveats=[
            "Same earthquakes; different station windows and propagation paths",
            "Measured/theory curves use the saved provisional model at each station",
            "Response agreement does not prove absolute hardware calibration",
            "Missing spectra are not filled; rejected events remain visibly labelled",
            "Body-arrival flags do not prove contamination or mode isolation",
        ],
    )
    (output / "review.json").write_text(json.dumps(audit, indent=2) + "\n")
    html = (
        '<!doctype html><meta charset="utf-8"><title>YL matched-event review</title>'
        "<style>body{font:16px sans-serif;max-width:1100px;margin:30px auto}"
        "img{max-width:100%}section{border-top:1px solid #aaa;margin-top:30px}</style>"
        "<h1>C07W / C08W / C09W matched-event review</h1>"
        f"<p>{len(shared)} matched saved events; "
        f"{all_candidate} candidates at all three."
        " No refitting, raw processing or calibration updates.</p>"
        "<p>Different station windows and paths remain a limitation. Open the saved "
        "arrival/waveform links to review packet isolation.</p>"
        '<h2>Nominal responses</h2><img src="responses.png" width="1000">'
        + "".join(cards)
    )
    index = output / "index.html"
    index.write_text(html)
    print(f"Matched saved events: {len(shared)}; all-three candidates: {all_candidate}")
    print(f"Review: {index}")
    print(f"Event table: {output / 'events.csv'}")
    if args.show:
        webbrowser.open(index.as_uri())


if __name__ == "__main__":
    main()
