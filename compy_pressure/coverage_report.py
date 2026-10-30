"""Descriptive coverage tables and pressure probability heat maps.

Every saved product is opened read-only and hash-checked before and after use.
Heat maps are drawn per network and pressure channel code, so different
pressure sensors (DPG, APG, hydrophone) and two streams of one sensor are
never pooled. Two separate weightings are written for each group:

* window-pooled: every saved PSD window counts once; the denominator at each
  frequency is the number of windows covering that frequency;
* station-balanced: each station's windows (all its epochs of that channel
  code) form one distribution normalized per frequency; stations are averaged
  with equal weight; the denominator is the number of stations covering that
  frequency.

No mean, median or percentile curve is computed or drawn. Values are kept,
including numerical zero-power floors and values outside the saved histogram
range; they are counted separately. Nothing here is a validated noise model.
"""

import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from obspy.imaging.cm import pqlx
from obspy.signal import PPSD

from compy_pressure.coverage import KEYS, digest, fresh_name, public, write_csv
from compy_pressure.coverage_run import is_done, reference_pool
from compy_pressure.plot_network_ppsd import (
    accumulate_histogram,
    add_references,
    spectrum,
)

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

# ObsPy stores 10*log10 of the smallest positive double for zero power.
ZERO_POWER_DB = 10 * np.log10(np.finfo(float).tiny)
SAVED_HISTOGRAM_DB = (-250.0, 100.0)
DETAIL_DB = (-100.0, 120.0)
# Keyword rule on the declared StationXML sensor label, first match wins.
# It only separates declared instrument families for plotting; it does not
# verify hardware. Exact labels are listed in the heat-map manifest.
SENSOR_FAMILIES = (
    ("apg", ("APG", "PAROSCI")),
    ("dpg", ("DIFFERENTIAL PRESSURE", "DPG", "COX-WEBB")),
    ("hydrophone", ("HYDROPHONE", "HTI")),
)


def sensor_label(channel):
    return str(channel.get("sensor_description") or channel.get("sensor_type") or "")


def sensor_family(label):
    upper = str(label).upper()
    for family, keywords in SENSOR_FAMILIES:
        if any(word in upper for word in keywords):
            return family
    return "unclassified"


def campaign_npz(campaign_dir, task_id, tasks):
    task = tasks.get(task_id)
    if not task:
        return None
    complete, state = is_done(campaign_dir, task)
    if not complete:
        return None
    paths = [o["path"] for o in state.get("outputs", []) if o["path"].endswith(".npz")]
    return paths[0] if len(paths) == 1 else None


def pressure_sources(plan, campaign_dir):
    """The product each pressure epoch would use now, with its origin."""
    tasks = {t["task_id"]: t for t in plan["tasks"]}
    rows = []
    for channel in plan["channels"]:
        if channel["channel_class"] != "pressure":
            continue
        target = channel.get("ppsd_target") or ""
        path, origin = None, "none"
        if target.startswith(("ppsd-", "ppsd_merge-")):
            path = campaign_npz(campaign_dir, target, tasks)
            origin = "campaign" if path else "pending"
            if not path and channel.get("ppsd_product"):
                path, origin = channel["ppsd_product"], "historical_partial"
        elif target:
            path = target
            origin = (
                "historical"
                if channel["ppsd_decision"] == "reuse_existing"
                else "historical_partial"
            )
        rows.append(
            {
                **{k: channel[k] for k in KEYS},
                "network": channel["network"],
                "station": channel["station"],
                "channel": channel["channel"],
                "sensor_label": sensor_label(channel),
                "sensor_family": sensor_family(sensor_label(channel)),
                "ppsd_decision": channel["ppsd_decision"],
                "npz": path or "",
                "origin": origin,
            }
        )
    return rows


def coverage_tables(plan, campaign_dir, sources):
    tasks = {t["task_id"]: t for t in plan["tasks"]}
    by_key = {tuple(s[k] for k in KEYS): s for s in sources}
    channel_rows, stations = [], defaultdict(Counter)
    for c in plan["channels"]:
        product = next(
            (p for p in c["products"] if p["npz"] == c.get("ppsd_product")), None
        )
        uncovered = product["uncovered_seconds"] if product else c["available_seconds"]
        states = []
        usable_states = []
        for task_id in c.get("ppsd_tasks", []):
            complete, state = (
                is_done(campaign_dir, tasks[task_id])
                if task_id in tasks
                else (
                    False,
                    None,
                )
            )
            states.append("complete" if complete else "pending")
            usable_states.append(
                complete
                and bool((state or {}).get("outputs"))
                and (state or {}).get("result_status")
                not in {"no_data", "insufficient_data"}
            )
        if c.get("ppsd_tasks") and all(usable_states):
            uncovered = 0.0
        covered = c["available_seconds"] - uncovered if (product or states) else 0.0
        if c["ppsd_decision"] in {"blocked", "review_required"} and not product:
            covered = 0.0
        source = by_key.get(tuple(c[k] for k in KEYS), {})
        channel_rows.append(
            {
                **{
                    k: c[k]
                    for k in (
                        "network",
                        "station",
                        "channel",
                        "Seed_id",
                        "channel_start_time",
                        "channel_end_time",
                        "inventory_file",
                        "channel_class",
                        "sample_rate_hz",
                        "response_status",
                        "files",
                        "available_start",
                        "available_end",
                        "ppsd_decision",
                        "ppsd_reason",
                    )
                },
                "available_days": round(c["available_seconds"] / 86400, 3),
                "product_covered_days": round(max(covered, 0.0) / 86400, 3),
                "coverage_basis": "processed_input_span_not_continuous_PSD_coverage",
                "product": c.get("ppsd_product") or "",
                "heatmap_source": source.get("npz", ""),
                "heatmap_origin": source.get("origin", ""),
                "campaign_task_states": ",".join(states),
                **{
                    f"check_{name}": entry["decision"]
                    for name, entry in c.get("checks", {}).items()
                },
            }
        )
        key = (c["network"], c["station"])
        stations[key][f"{c['channel_class']}_epochs"] += 1
        stations[key][f"{c['channel_class']}_available_days"] += (
            c["available_seconds"] / 86400
        )
        stations[key][f"{c['channel_class']}_covered_days"] += max(covered, 0) / 86400
        stations[key][f"ppsd_{c['ppsd_decision']}"] += 1
        for name, entry in c.get("checks", {}).items():
            stations[key][f"{name}_{entry['decision']}"] += 1
    station_rows = [
        {"network": net, "station": sta}
        | {k: round(v, 3) if isinstance(v, float) else v for k, v in sorted(c.items())}
        for (net, sta), c in sorted(stations.items())
    ]
    return channel_rows, station_rows


def map_counts(state, edges, lo, hi):
    """Station counts on the union frequency grid, without interpolation."""
    centers = np.sqrt(edges[:-1] * edges[1:])
    index = np.searchsorted(state["edges"], centers, side="right") - 1
    keep = (index >= 0) & (index < state["counts"].shape[1])
    block = np.zeros((state["hi"] - state["lo"], len(centers)))
    block[:, keep] = state["counts"][:, index[keep]]
    out = np.zeros((hi - lo, len(centers)))
    out[state["lo"] - lo : state["hi"] - lo] = block
    return out


def draw(path_stub, edges, lo, hi, occurrence, title, note, max_percentage):
    figure, axis = plt.subplots(figsize=(14, 8))
    figure.subplots_adjust(left=0.08, right=0.97, bottom=0.2, top=0.86)
    db_edges = np.arange(lo, hi + 1)
    mesh = axis.pcolormesh(
        edges,
        db_edges,
        np.ma.masked_equal(occurrence, 0),
        cmap=pqlx,
        shading="flat",
        rasterized=True,
        vmin=0,
        vmax=max_percentage,
    )
    figure.colorbar(mesh, ax=axis, label="Occurrence [% per 1-dB bin]", pad=0.02)
    add_references(axis, True, edges[0], edges[-1])
    axis.set(
        xscale="log",
        xlim=(edges[0], edges[-1]),
        xlabel="Frequency [Hz]",
        ylabel=r"Pressure PSD [dB re 1 Pa$^2$/Hz]",
    )
    axis.grid(which="both", alpha=0.2)
    axis.legend(loc="lower left", fontsize=8, framealpha=0.9)
    axis.set_title(title, fontsize=10)
    figure.text(0.08, 0.03, note, fontsize=8)
    full = Path(f"{path_stub}_full_range.png")
    figure.savefig(full, dpi=130)
    axis.set_ylim(*DETAIL_DB)
    axis.set_title(
        title + f"\nDetail view {DETAIL_DB[0]:g} to {DETAIL_DB[1]:g} dB; "
        "counts outside it are retained (full-range PNG and NPZ)",
        fontsize=10,
    )
    detail = Path(f"{path_stub}.png")
    figure.savefig(detail, dpi=130)
    plt.close(figure)
    return [str(detail), str(full)]


def heatmaps(sources, output, max_percentage=30.0):
    output.mkdir(parents=True, exist_ok=False)
    groups = defaultdict(list)
    for source in sources:
        if source["npz"]:
            key = (source["network"], source["channel"], source["sensor_family"])
            groups[key].append(source)
    manifest = {"groups": [], "inputs": []}
    for (network, code, family), items in sorted(groups.items()):
        pooled, stations, stats = {}, defaultdict(dict), defaultdict(Counter)
        for item in items:
            path = Path(item["npz"])
            before = digest(path)
            ppsd = PPSD.load_npz(str(path))
            f, values, _ = spectrum(
                ppsd, {"Seed_id": item["Seed_id"], "ppsd_units": "Pa^2/Hz"}
            )
            edges = 1 / np.asarray(ppsd.period_xedges)[::-1]
            accumulate_histogram(pooled, f, values, edges)
            accumulate_histogram(stations[item["station"]], f, values, edges)
            floor = values <= ZERO_POWER_DB + 0.01
            counts = {
                "windows": len(values),
                "values": int(values.size),
                "zero_power_floor_values": int(floor.sum()),
                "windows_with_zero_power_floor": int(floor.any(axis=1).sum()),
                "below_saved_histogram_values": int(
                    ((values < SAVED_HISTOGRAM_DB[0]) & ~floor).sum()
                ),
                "above_saved_histogram_values": int(
                    (values > SAVED_HISTOGRAM_DB[1]).sum()
                ),
                "outside_detail_view_values": int(
                    ((values < DETAIL_DB[0]) | (values > DETAIL_DB[1])).sum()
                ),
            }
            stats[item["station"]].update(counts)
            if digest(path) != before:
                raise RuntimeError(f"Saved PPSD changed while reading: {path}")
            manifest["inputs"].append(
                item
                | counts
                | {
                    "sha256": before,
                    "first_window": str(ppsd.times_processed[0]),
                    "last_window": str(ppsd.times_processed[-1]),
                    "frequency_min_hz": float(f[0]),
                    "frequency_max_hz": float(f[-1]),
                }
            )
        edges, lo, hi = pooled["edges"], pooled["lo"], pooled["hi"]
        coverage = pooled["counts"].sum(axis=0)
        pooled_occurrence = 100.0 * pooled["counts"] / np.maximum(coverage, 1)
        balanced = np.zeros((hi - lo, len(edges) - 1))
        covering = np.zeros(len(edges) - 1)
        for state in stations.values():
            counts = map_counts(state, edges, lo, hi)
            total = counts.sum(axis=0)
            present = total > 0
            balanced[:, present] += counts[:, present] / total[present]
            covering += present
        balanced_occurrence = 100.0 * balanced / np.maximum(covering, 1)
        names = {k for counts in stats.values() for k in counts}
        totals = {k: sum(c.get(k, 0) for c in stats.values()) for k in names}
        floors = totals["zero_power_floor_values"]
        stub = output / f"{network}_{code}_{family}"
        label = f"{network} {code} ({family})"
        windows = int(pooled["windows"])
        common = (
            f"Zero-power floor values (ObsPy {ZERO_POWER_DB:.2f} dB): {floors:,} "
            f"in {totals['windows_with_zero_power_floor']:,} windows; "
            f"below/above saved histogram {totals['below_saved_histogram_values']:,}/"
            f"{totals['above_saved_histogram_values']:,} values. All retained.\n"
            "Saved nominal-response spectra; no gain, acceptance filter, "
            "usable-band censoring or outlier removal. Not a validated noise model. "
            "Overlapping windows are correlated; counts are not independent samples."
        )
        figures = draw(
            f"{stub}_window_pooled",
            edges,
            lo,
            hi,
            pooled_occurrence,
            f"{label}: {len(stations)} stations, {windows:,} saved PSD "
            "windows - window-pooled occurrence\nDenominator per frequency: "
            "windows covering it; longer records weigh more",
            common,
            max_percentage,
        )
        figures += draw(
            f"{stub}_station_balanced",
            edges,
            lo,
            hi,
            balanced_occurrence,
            f"{label}: {len(stations)} stations - station-balanced "
            "occurrence\nEach station normalized per frequency, then equal "
            "weight; denominator: stations covering it",
            common,
            max_percentage,
        )
        np.savez_compressed(
            f"{stub}_heatmap.npz",
            pooled_counts=pooled["counts"],
            windows_per_frequency=coverage,
            station_balanced_probability=balanced_occurrence / 100.0,
            stations_per_frequency=covering,
            frequency_edges_hz=edges,
            db_edges=np.arange(lo, hi + 1),
        )
        manifest["groups"].append(
            {
                "network": network,
                "channel": code,
                "sensor_family": family,
                "declared_sensor_labels": dict(
                    Counter(item["sensor_label"] for item in items)
                ),
                "stations": len(stations),
                "windows": windows,
                "db_range_retained": [int(lo), int(hi)],
                "figures": figures,
                "per_station": {s: dict(c) for s, c in sorted(stats.items())},
                **dict(totals),
            }
        )
    manifest["weightings"] = {
        "window_pooled": "each saved window once; denominator = windows covering "
        "the frequency",
        "station_balanced": "each station normalized per frequency; equal station "
        "weight; denominator = stations covering the frequency",
    }
    manifest["grouping"] = (
        "network, channel code and sensor family from the declared StationXML "
        f"label (keyword rule {SENSOR_FAMILIES}); labels are not hardware checks"
    )
    manifest["colormap"] = "ObsPy pqlx"
    manifest["summary_curves"] = "none (no mean, median or percentile)"
    with (output / "heatmap_manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, allow_nan=False, default=str)
    return manifest


def task_state(complete, state):
    """A recorded completion whose outputs or nearby pool changed is stale."""
    if complete:
        return "complete"
    status = (state or {}).get("status", "not_started")
    return "stale" if status == "complete" else status


def write_report(plan, campaign_dir, output=None, *, max_percentage=30.0):
    campaign_dir = Path(campaign_dir)
    output = Path(output or campaign_dir / "reports" / fresh_name(plan["plan_id"]))
    output.mkdir(parents=True, exist_ok=False)
    sources = pressure_sources(plan, campaign_dir)
    channels, stations = coverage_tables(plan, campaign_dir, sources)
    write_csv(output / "coverage_channels.csv", channels, list(channels[0]))
    fields = sorted({k for row in stations for k in row} - {"network", "station"})
    write_csv(
        output / "coverage_stations.csv", stations, ["network", "station"] + fields
    )
    tasks = []
    # A nearby result is current only for the present reference pool.
    pool = reference_pool(
        campaign_dir, plan, {tuple(c[k] for k in KEYS): c for c in plan["channels"]}
    )[1]
    for task in plan["tasks"]:
        complete, state = is_done(
            campaign_dir, task, pool if task["kind"] == "nearby" else None
        )
        tasks.append(
            {
                "task_id": task["task_id"],
                "kind": task["kind"],
                "seed_id": task["seed_id"],
                "interval": task.get("interval", ""),
                "reason": task["reason"],
                "dependencies": task["dependencies"],
                "state": task_state(complete, state),
                "result_status": (state or {}).get("result_status", ""),
                "error": (state or {}).get("error", ""),
            }
        )
    write_csv(
        output / "tasks_status.csv", tasks, list(tasks[0]) if tasks else ["task_id"]
    )
    manifest = heatmaps(sources, output / "heatmaps", max_percentage)
    summary = {
        "created_utc": datetime.now(UTC).isoformat(),
        "plan_id": plan["plan_id"],
        "plan_summary": plan["summary"],
        "task_states": dict(Counter(t["state"] for t in tasks)),
        "heatmap_sources": dict(Counter(s["origin"] for s in sources)),
        "heatmap_groups": [
            {
                k: g[k]
                for k in (
                    "network",
                    "channel",
                    "sensor_family",
                    "stations",
                    "windows",
                    "zero_power_floor_values",
                    "below_saved_histogram_values",
                    "above_saved_histogram_values",
                )
            }
            for g in manifest["groups"]
        ],
        "gain_applied": False,
        "scores_or_thresholds_assigned": False,
    }
    (output / "summary.json").write_text(
        json.dumps(public(summary), indent=2, allow_nan=False) + "\n"
    )
    return output


def shorten(value, plan):
    """Portable paths for curated, shareable tables."""
    text = str(value)
    for prefix, label in (
        (plan["config"]["data_directory"], "<DATA>"),
        (plan["config"]["source_root"], "<SOURCE>"),
        (plan["campaign_dir"], "<CAMPAIGN>"),
        (str(Path.home()), "~"),
    ):
        text = text.replace(prefix, label)
    return text


def write_curated(plan, report_dir, curated_dir):
    """Small tables and detail heat maps for review; full products stay local."""
    import shutil

    curated = Path(curated_dir)
    curated.mkdir(parents=True, exist_ok=True)
    rows = []
    for c in plan["channels"]:
        if c["channel_class"] != "pressure":
            continue
        rows.append(
            {
                **{k: c[k] for k in ("network", "station", "channel", "Seed_id")},
                "channel_start_time": c["channel_start_time"],
                "channel_end_time": c["channel_end_time"],
                "inventory_file": shorten(c["inventory_file"], plan),
                "sample_rate_hz": c["sample_rate_hz"],
                "response_status": c.get("response_status", ""),
                "available_start": c["available_start"],
                "available_end": c["available_end"],
                "available_days": round(c["available_seconds"] / 86400, 3),
                "available_intervals": len(c["available_intervals"]),
                "files": c["files"],
                "file_layouts": c["file_layouts"],
                "ppsd_decision": c["ppsd_decision"],
                "ppsd_reason": shorten(c["ppsd_reason"], plan),
                "existing_product": shorten(c.get("ppsd_product") or "", plan),
                "planned_ppsd_tasks": len(c.get("ppsd_tasks", [])),
                **{
                    f"check_{name}": entry["decision"]
                    for name, entry in c.get("checks", {}).items()
                },
                **{
                    f"history_{name}": ",".join(
                        sorted({e["state"] for e in entry.get("evidence", [])})
                    )
                    for name, entry in c.get("checks", {}).items()
                    if name in ("rayleigh", "microseism", "local_p")
                },
                **{
                    f"colocated_{name}": ",".join(entry.get("colocated_evidence", []))
                    for name, entry in c.get("checks", {}).items()
                    if name in ("rayleigh", "microseism", "local_p")
                },
            }
        )
    write_csv(curated / "pressure_epochs.csv", rows, list(rows[0]))
    intervals = {
        f"{c['Seed_id']}|{c['channel_start_time']}": c["available_intervals"]
        for c in plan["channels"]
        if c["channel_class"] == "pressure"
    }
    (curated / "pressure_available_intervals.json").write_text(
        json.dumps(intervals, indent=1) + "\n"
    )
    groups = defaultdict(Counter)
    days = defaultdict(float)
    for c in plan["channels"]:
        key = (c["network"], c["channel"], c["channel_class"])
        groups[key]["epochs"] += 1
        groups[key][c["ppsd_decision"]] += 1
        days[key] += c["available_seconds"] / 86400
    decisions = sorted({d for g in groups.values() for d in g} - {"epochs"})
    summary_rows = [
        {
            "network": net,
            "channel": cha,
            "class": cls,
            "available_days": round(days[(net, cha, cls)], 1),
            **{k: groups[(net, cha, cls)].get(k, 0) for k in ["epochs", *decisions]},
        }
        for net, cha, cls in sorted(groups)
    ]
    write_csv(curated / "channel_summary.csv", summary_rows, list(summary_rows[0]))
    blocked = Counter()
    for c in plan["channels"]:
        if c["ppsd_decision"] in {
            "blocked",
            "review_required",
            "reuse_existing_partial",
        }:
            reason = shorten(c["ppsd_reason"], plan)
            blocked[(c["network"], c["channel_class"], c["ppsd_decision"], reason)] += 1
    write_csv(
        curated / "blocked_prerequisites.csv",
        [
            {"network": n, "class": k, "decision": d, "reason": r, "epochs": v}
            for (n, k, d, r), v in sorted(blocked.items())
        ],
        ["network", "class", "decision", "reason", "epochs"],
    )
    archive = plan["archive"]
    issues = Counter(
        (f["filename_seed_id"], f["issue"].split(":")[0])
        for f in archive["files_with_issues"]
    )
    duplicates = Counter(d["seed_id"] for d in archive["excluded_duplicate_copies"])
    inventory_only = Counter(
        (e["Seed_id"].split(".")[0], e["channel_class"])
        for e in plan["inventory_only_epochs"]
    )
    (curated / "archive_and_inventory.json").write_text(
        json.dumps(
            {
                "data_directory": "<DATA>",
                "recognized_files": archive["recognized_files"],
                "ignored_files": archive["ignored_files"],
                "layouts": archive["layouts"],
                "header_streams": archive["streams"],
                "files_with_issues": [
                    {"filename_seed_id": s, "issue": i, "files": n}
                    for (s, i), n in sorted(issues.items())
                ],
                "excluded_duplicate_copies_by_stream": dict(sorted(duplicates.items())),
                "unmatched_streams": archive["unmatched_streams"],
                "inventory_sources": [
                    s | {"path": shorten(s["path"], plan)}
                    for s in plan["inventory_sources"]
                ],
                "inventory_epochs_without_local_data": {
                    f"{n}:{k}": v for (n, k), v in sorted(inventory_only.items())
                },
                "evidence_unjoined": plan["evidence_unjoined"],
                "evidence_errors": plan["evidence_errors"],
            },
            indent=1,
            default=str,
        )
        + "\n"
    )
    summary = json.loads((Path(report_dir) / "summary.json").read_text())
    summary["config_sha256"] = plan["config"]["config_sha256"]
    summary["scope"] = plan["scope"]
    (curated / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    figures = curated / "heatmaps"
    figures.mkdir(exist_ok=True)
    for png in sorted((Path(report_dir) / "heatmaps").glob("*.png")):
        if not png.stem.endswith("_full_range"):
            shutil.copy2(png, figures / png.name)
    return curated
