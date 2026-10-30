"""Inventory the whole local archive and plan only missing, reproducible work.

``plan`` reads waveform headers, StationXML, registered PPSD products and saved
check evidence. It writes only below one campaign folder (a header cache and a
fresh plan folder); every input stays unchanged. ``run`` (dry run unless
``--execute``) and ``report`` live in ``coverage_run`` and ``coverage_report``.

Completeness is judged from available waveform intervals, channel epochs,
response identity, product parameters and recorded source files, never from
file existence or code fingerprints alone. Products and evidence come only from
explicitly configured registries and roots; nothing is selected by recency.
No gain is applied, no score/threshold/rating is assigned or adopted.
"""

import argparse
import contextlib
import csv
import glob
import hashlib
import json
import math
import os
import sys
import warnings
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import numpy as np
from obspy import UTCDateTime, read, read_inventory
from obspy.io.mseed.util import get_record_information

from compy_pressure import quality_review
from compy_pressure.archive import waveform_seed_id
from compy_pressure.inventory import EPOCH_COLUMNS
from compy_pressure.ppsd import atomic_json, is_pressure_unit, response_handling

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

ROOT = Path(__file__).resolve().parents[1] / "pressure"
CONFIG = ROOT / "docs/coverage_sources.json"
SCHEMA = "bruitfm-coverage-sources-v1"
PLAN_VERSION = "coverage-plan-v1"
KEYS = EPOCH_COLUMNS + ["inventory_file"]
# Equal to the `compy_pressure.network ppsd` defaults used for every saved product.
PPSD_DEFAULTS = {
    "ppsd_length": 10800.0,
    "overlap": 0.5,
    "minimum_frequency": 0.001,
    "chunk_seconds": 86400.0,
}
PRODUCT_PARAMETERS = {
    "ppsd_length": "ppsd_length_seconds",
    "overlap": "overlap",
    "minimum_frequency": "minimum_frequency_hz",
}
FIXED_PRODUCT_CONFIG = {"db_bins": [-250.0, 100.0, 1.0], "skip_on_gaps": True}
# Check name -> (network stage, master/result prefix) of the existing estimator.
SPECTRAL_CHECKS = {
    "expected_pressure": ("pressure", "pressure_expected"),
    "nearby": ("nearby", "nearby"),
    "infragravity": ("infragravity", "infragravity"),
}
CALIBRATION_CHECKS = ("rayleigh", "microseism", "local_p")
SEISMIC_UNITS = {"M", "M/S", "M/S**2", "M/S^2", "M/S/S"}
SIGNATURE_FREQUENCIES = np.logspace(-3, 1, 41)
SCOPES = ("pressure", "all-channels")


def digest(path):
    with Path(path).open("rb") as stream:
        # hashlib.file_digest needs Python 3.11; ComPy supports 3.10
        checksum = hashlib.sha256()
        for block in iter(lambda: stream.read(1 << 20), b""):
            checksum.update(block)
        return checksum.hexdigest()


def json_hash(value):
    encoded = json.dumps(value, sort_keys=True, allow_nan=False, default=str)
    return hashlib.sha256(encoded.encode()).hexdigest()


def utc_now():
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


def fresh_name(plan_id):
    """Unique folder name: timestamp, plan identifier and a random suffix."""
    return f"{utc_now()}_{plan_id[:8]}_{uuid4().hex[:6]}"


def timestamp(value):
    return float(UTCDateTime(value).timestamp)


def iso(seconds):
    return str(UTCDateTime(seconds))


def text(value):
    """Format like the CSV master: blanks for missing values, str() otherwise."""
    return "" if value is None else str(value)


# --------------------------------------------------------------------------
# Intervals: half-open [start, end) POSIX seconds.


def merge_intervals(intervals, tolerance=0.0):
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if merged and start <= merged[-1][1] + tolerance:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(span) for span in merged]


def subtract_intervals(intervals, removed):
    result = []
    removed = merge_intervals(removed)
    for start, end in merge_intervals(intervals):
        cursor = start
        for cut_start, cut_end in removed:
            if cut_end <= cursor or cut_start >= end:
                continue
            if cut_start > cursor:
                result.append((cursor, cut_start))
            cursor = max(cursor, cut_end)
        if cursor < end:
            result.append((cursor, end))
    return result


def clip_intervals(intervals, start=None, end=None):
    low = -math.inf if start is None else start
    high = math.inf if end is None else end
    return [(max(s, low), min(e, high)) for s, e in intervals if e > low and s < high]


def total_seconds(intervals):
    return float(sum(end - start for start, end in intervals))


def interval_text(intervals):
    return [[iso(start), iso(end)] for start, end in intervals]


# --------------------------------------------------------------------------
# Configuration


def load_config(path=CONFIG, source_root=ROOT):
    """Resolve explicit sources; relative paths are below the source checkout."""
    path = Path(path).expanduser().resolve()
    config = json.loads(path.read_text())
    if config.get("schema") != SCHEMA:
        raise ValueError(f"Unsupported coverage source schema in {path}")
    root = Path(source_root).expanduser().resolve()

    def where(value):
        candidate = Path(value).expanduser()
        return str(candidate if candidate.is_absolute() else root / candidate)

    evidence = config.get("evidence", {})
    required = [
        where(p)
        for p in config.get("product_registries", [])
        + config.get("product_roots", [])
        + evidence.get("masters", [])
        + [evidence[k] for k in CALIBRATION_CHECKS if evidence.get(k)]
    ]
    missing = [p for p in required if not Path(p).exists()]
    if missing:
        raise FileNotFoundError(
            "Configured saved sources are missing (wrong --source-root?): "
            + ", ".join(missing)
        )
    parameters = PPSD_DEFAULTS | config.get("ppsd_parameters", {})
    if set(parameters) != set(PPSD_DEFAULTS):
        raise ValueError("Unknown PPSD parameter in coverage sources")
    return {
        "config_file": str(path),
        "config_sha256": digest(path),
        "source_root": str(root),
        "networks": list(config["networks"]),
        "data_directory": str(Path(where(config["data_directory"])).resolve()),
        "inventory_sources": {
            net: [where(p) for p in patterns]
            for net, patterns in config["inventory_sources"].items()
            if net in config["networks"]
        },
        "layout_preference": list(config.get("layout_preference", ["sds"])),
        "allowed_location_aliases": [
            tuple(pair) for pair in config.get("allowed_location_aliases", [])
        ],
        "gap_tolerance_seconds": float(config.get("gap_tolerance_seconds", 1.0)),
        "ppsd_parameters": {k: float(v) for k, v in parameters.items()},
        "product_registries": [where(p) for p in config.get("product_registries", [])],
        "product_roots": [where(p) for p in config.get("product_roots", [])],
        "evidence": {
            "masters": [where(p) for p in evidence.get("masters", [])],
            **{
                kind: where(evidence[kind])
                for kind in CALIBRATION_CHECKS
                if evidence.get(kind)
            },
        },
    }


def expand_sources(patterns):
    paths = []
    for pattern in patterns:
        found = sorted(glob.glob(pattern)) if glob.has_magic(pattern) else [pattern]
        paths.extend(str(Path(p).resolve()) for p in found if Path(p).is_file())
    return list(dict.fromkeys(paths))


# --------------------------------------------------------------------------
# Archive discovery


def filename_identity(name):
    """Return (layout, filename SEED ID) for the supported daily file names."""
    name = Path(name).name
    seed = waveform_seed_id(name)
    if seed is not None:
        return ("daily" if name.lower().endswith(".mseed") else "sds"), seed
    if name.lower().endswith(".mseed") and len(name[:-6].split(".")) == 7:
        seed = waveform_seed_id(name[:-6])
        if seed is not None:
            return "sds_mseed", seed
    return None, None


def probe_header(path, full=False):
    """Read the first/last MiniSEED record headers; fall back to all headers.

    The fast path gives the file span only. Internal gaps are found by the
    full-header option or later by PPSD processing, which skips gapped windows.
    """
    path = str(path)
    method = "full_header" if full else "first_last_record"
    try:
        if not full:
            first = get_record_information(path)
            size, length = first["filesize"], first["record_length"]
            if length and size % length == 0:
                last = (
                    first
                    if size == length
                    else get_record_information(path, offset=size - length)
                )
                ids = [
                    tuple(
                        str(r[k]).strip()
                        for k in ("network", "station", "location", "channel")
                    )
                    for r in (first, last)
                ]
                rate = float(first["samp_rate"])
                if (
                    ids[0] == ids[1]
                    and rate > 0
                    and float(last["samp_rate"]) == rate
                    and last["endtime"] >= first["starttime"]
                ):
                    return {
                        "seed_id": ".".join(ids[0]),
                        "sample_rate_hz": rate,
                        "start": str(first["starttime"]),
                        "end": str(last["endtime"] + 1 / rate),
                        "method": method,
                    }
            method = "full_header_fallback"
        traces = read(path, format="MSEED", headonly=True)
        ids = sorted({t.id for t in traces})
        rates = sorted({float(t.stats.sampling_rate) for t in traces})
        if len(ids) != 1 or len(rates) != 1:
            return {"error": f"mixed identifiers/rates {ids} {rates}", "method": method}
        spans = merge_intervals(
            [
                (
                    t.stats.starttime.timestamp,
                    (t.stats.endtime + t.stats.delta).timestamp,
                )
                for t in traces
            ]
        )
        return {
            "seed_id": ids[0],
            "sample_rate_hz": rates[0],
            "start": iso(spans[0][0]),
            "end": iso(spans[-1][1]),
            "method": method,
            "traces": len(traces),
            "intervals": interval_text(spans) if len(spans) > 1 else None,
        }
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}", "method": method}


class HeaderCache:
    """Header probes keyed by path, reused only for an unchanged size/mtime."""

    def __init__(self, path):
        self.path = Path(path)
        self.entries = json.loads(self.path.read_text()) if self.path.is_file() else {}
        self.hits = self.misses = 0

    def get(self, path, full=False):
        stat = Path(path).stat()
        cached = self.entries.get(str(path))
        if (
            cached
            and cached["size"] == stat.st_size
            and cached["mtime_ns"] == stat.st_mtime_ns
            and (not full or cached["info"].get("method", "").startswith("full"))
        ):
            self.hits += 1
            return cached["info"], stat
        self.misses += 1
        info = probe_header(path, full)
        self.entries[str(path)] = {
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "info": info,
        }
        return info, stat

    def save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        # A concurrent read-only planner may share this cache. Its snapshot
        # can safely replace ours, but it must not steal our temporary file.
        temporary = self.path.with_name(f".{self.path.name}.{uuid4().hex}.tmp")
        temporary.write_text(json.dumps(self.entries))
        temporary.replace(self.path)


def identity_issue(filename_seed, header_seed, aliases):
    if header_seed == filename_seed:
        return None
    a, b = filename_seed.split("."), header_seed.split(".")
    if len(b) == 4 and a[:2] == b[:2] and a[3] == b[3] and (a[2], b[2]) in aliases:
        return None
    return f"filename_header_mismatch: {filename_seed} versus {header_seed}"


def scan_archive(config, cache, *, full_headers=False, progress=True):
    """Every recognized waveform file below the data directory, header-checked."""
    root = Path(config["data_directory"])
    if not root.is_dir():
        raise FileNotFoundError(f"Waveform directory does not exist: {root}")
    networks = set(config["networks"])
    aliases = set(config["allowed_location_aliases"])
    files, ignored = [], Counter()
    for directory, folders, names in os.walk(root):
        folders.sort()
        for name in sorted(names):
            path = Path(directory) / name
            layout, filename_seed = filename_identity(name)
            if layout is None:
                ignored[Path(name).suffix.lower() or name] += 1
                continue
            if filename_seed.split(".")[0] not in networks:
                ignored["out_of_scope_network"] += 1
                continue
            info, stat = cache.get(path, full_headers)
            record = {
                "path": str(path),
                "layout": layout,
                "filename_seed_id": filename_seed,
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                **info,
            }
            if info.get("error"):
                record["issue"] = f"unreadable_header: {info['error']}"
            else:
                issue = identity_issue(filename_seed, info["seed_id"], aliases)
                if issue:
                    record["issue"] = issue
            files.append(record)
            if progress and len(files) % 10000 == 0:
                print(f"Indexed {len(files)} waveform files", flush=True)
                cache.save()
    cache.save()
    return files, dict(ignored)


def select_copies(files, preference, tolerance):
    """Keep the preferred layout; add other copies only where it has no data.

    Different copies of one channel are never merged where they overlap.
    Overlaps within one layout (records crossing midnight) are normal.
    """
    rank = {layout: index for index, layout in enumerate(preference)}
    groups = defaultdict(list)
    for record in files:
        if not record.get("issue"):
            groups[record["seed_id"]].append(record)
    selected, excluded = {}, []
    for seed, group in sorted(groups.items()):
        layouts = sorted({r["layout"] for r in group}, key=lambda x: rank.get(x, 99))
        kept, spans = [], []
        for layout in layouts:
            accepted = []
            for record in sorted(group, key=lambda r: (r["start"], r["path"])):
                if record["layout"] != layout:
                    continue
                start, end = timestamp(record["start"]), timestamp(record["end"])
                if any(
                    start < s_end - tolerance and end > s_start + tolerance
                    for s_start, s_end in spans
                ):
                    excluded.append(
                        record | {"exclusion": f"overlaps_preferred_copy:{layouts[0]}"}
                    )
                    continue
                accepted.append(record)
            kept.extend(accepted)
            spans.extend((timestamp(r["start"]), timestamp(r["end"])) for r in accepted)
        selected[seed] = sorted(kept, key=lambda r: (r["start"], r["path"]))
    return selected, excluded


def file_intervals(record):
    if record.get("intervals"):
        return [(timestamp(a), timestamp(b)) for a, b in record["intervals"]]
    return [(timestamp(record["start"]), timestamp(record["end"]))]


# --------------------------------------------------------------------------
# StationXML epochs and response identity


def channel_class(units, code):
    units = str(units or "").upper().replace(" ", "")
    if is_pressure_unit(units):
        return "pressure"
    if units in SEISMIC_UNITS:
        if code.endswith("Z"):
            return "vertical"
        if code[-1] in "12NE":
            return "horizontal"
        return "other_seismic"
    return "other"


@contextlib.contextmanager
def quiet_c_stderr():
    """evalresp prints FIR normalization notes from C; keep plan output legible."""
    sys.stderr.flush()
    saved = os.dup(2)
    try:
        with open(os.devnull, "w") as null:
            os.dup2(null.fileno(), 2)
            yield
    finally:
        os.dup2(saved, 2)
        os.close(saved)


def response_signature(channel):
    """Amplitude and phase on a fixed grid below Nyquist, plus sensitivity."""
    handling, units = response_handling(channel)
    rate = float(channel.sample_rate)
    frequencies = SIGNATURE_FREQUENCIES[SIGNATURE_FREQUENCIES < 0.45 * rate]
    with warnings.catch_warnings(), quiet_c_stderr():
        warnings.simplefilter("ignore")
        values = channel.response.get_evalresp_response_for_frequencies(
            frequencies, output="DEF" if handling else "VEL"
        )
    sensitivity = channel.response.instrument_sensitivity
    payload = {
        "psd_units": units,
        "input_units": str(sensitivity.input_units).upper(),
        "sensitivity": float(sensitivity.value),
        "frequencies_hz": frequencies.tolist(),
        "amplitude": np.abs(values).tolist(),
        "phase_deg": np.angle(values, deg=True).tolist(),
    }
    rounded = {
        **payload,
        "sensitivity": float(f"{payload['sensitivity']:.6g}"),
        "amplitude": [float(f"{a:.6g}") for a in payload["amplitude"]],
        "phase_deg": [round(p, 2) for p in payload["phase_deg"]],
    }
    return json_hash(rounded)[:16], payload


def same_response(a, b, rtol=1e-6, phase_tolerance_deg=0.01):
    if a is None or b is None:
        return None
    if a["psd_units"] != b["psd_units"] or a["frequencies_hz"] != b["frequencies_hz"]:
        return False
    if not math.isclose(a["sensitivity"], b["sensitivity"], rel_tol=rtol):
        return False
    phase = np.angle(
        np.exp(1j * np.radians(np.subtract(a["phase_deg"], b["phase_deg"]))), deg=True
    )
    return bool(
        np.allclose(a["amplitude"], b["amplitude"], rtol=rtol, atol=0)
        and np.all(np.abs(phase) <= phase_tolerance_deg)
    )


def epoch_row(network, station, channel, path, sha256):
    sensitivity = channel.response.instrument_sensitivity if channel.response else None
    sensor = channel.sensor
    depth = (
        -station.elevation
        if station.elevation is not None and station.elevation < 0
        else None
    )
    return {
        "Seed_id": f"{network.code}.{station.code}.{channel.location_code}."
        f"{channel.code}",
        "channel_start_time": text(channel.start_date),
        "channel_end_time": text(channel.end_date),
        "inventory_file": str(path),
        "inventory_sha256": sha256,
        "sample_rate_hz": text(channel.sample_rate),
        "Lat": text(channel.latitude),
        "Lon": text(channel.longitude),
        "Profondeur": text(depth),
        "azimuth_deg": text(channel.azimuth),
        "dip_deg": text(channel.dip),
        "station_name": text(station.site.name if station.site else None),
        "sensor_type": text(sensor.type if sensor else None),
        "sensor_model": text(sensor.model if sensor else None),
        "sensor_description": text(sensor.description if sensor else None),
        "sensitivity": text(sensitivity.value if sensitivity else None),
        "input_units": text(sensitivity.input_units if sensitivity else None),
        "output_units": text(sensitivity.output_units if sensitivity else None),
        "response_stage_count": text(
            len(channel.response.response_stages) if channel.response else 0
        ),
    }


class ResponseLibrary:
    """Exact channel epochs of any StationXML file, parsed once per file."""

    def __init__(self):
        self.files = {}

    def epochs(self, path):
        path = str(Path(path).resolve())
        if path not in self.files:
            sha = digest(path)
            inventory = read_inventory(path, format="STATIONXML")
            entries = []
            for network in inventory:
                for station in network:
                    for channel in station:
                        entries.append((network, station, channel))
            self.files[path] = (sha, entries)
        return self.files[path]

    def signature(self, path, seed, start, end):
        """Response signature of one exact epoch, or (None, reason)."""
        if not Path(path).is_file():
            return None, "missing_inventory_file"
        _, entries = self.epochs(path)
        found = [
            c
            for n, s, c in entries
            if f"{n.code}.{s.code}.{c.location_code}.{c.code}" == seed
            and text(c.start_date) == start
            and text(c.end_date) == end
        ]
        if len(found) != 1:
            return None, f"expected_one_epoch_found_{len(found)}"
        try:
            return response_signature(found[0])[1], None
        except Exception as exc:
            return None, f"unusable_response: {exc}"


def read_epochs(config, library, seeds_with_data):
    """Every channel epoch of the configured StationXML, primary source first."""
    epochs, sources = [], []
    for network, patterns in config["inventory_sources"].items():
        for order, path in enumerate(expand_sources(patterns)):
            sha, entries = library.epochs(path)
            count = 0
            for net, sta, cha in entries:
                if net.code != network:
                    continue
                row = epoch_row(net, sta, cha, path, sha)
                row["network"] = network
                row["source_order"] = order
                row["channel_class"] = channel_class(row["input_units"], cha.code)
                if row["Seed_id"] in seeds_with_data:
                    try:
                        handling, units = response_handling(cha)
                        signature, payload = response_signature(cha)
                        row.update(
                            response_status="usable",
                            psd_units=units,
                            response_signature=signature,
                            _response=payload,
                        )
                    except Exception as exc:
                        row.update(
                            response_status="unusable",
                            response_error=f"{type(exc).__name__}: {exc}",
                        )
                epochs.append(row)
                count += 1
            sources.append(
                {"network": network, "path": path, "sha256": sha, "epochs": count}
            )
    return epochs, sources


def assign_primary(epochs):
    """Epochs of the first source defining a SEED ID are primary.

    Later sources never add competing epochs for that SEED ID; their overlapping
    epochs are listed with an explicit response comparison instead.
    """
    groups = defaultdict(list)
    for row in epochs:
        groups[row["Seed_id"]].append(row)
    primary = []
    for rows in groups.values():
        first = min(r["source_order"] for r in rows)
        heads = [r for r in rows if r["source_order"] == first]
        others = [r for r in rows if r["source_order"] != first]
        for head in heads:
            head["alternative_sources"] = [
                {
                    "inventory_file": other["inventory_file"],
                    "epoch": [other["channel_start_time"], other["channel_end_time"]],
                    "same_epoch": (
                        other["channel_start_time"],
                        other["channel_end_time"],
                    )
                    == (head["channel_start_time"], head["channel_end_time"]),
                    "same_response": same_response(
                        head.get("_response"), other.get("_response")
                    ),
                }
                for other in others
                if overlaps(*epoch_bounds(head), *epoch_bounds(other))
            ]
            primary.append(head)
    return sorted(primary, key=lambda r: (r["Seed_id"], r["channel_start_time"]))


# --------------------------------------------------------------------------
# Existing products


def read_registry(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def discover_products(config, library):
    """Products pointed to by explicit registries or inside explicit roots."""
    registered = defaultdict(list)
    for registry in config["product_registries"]:
        for row in read_registry(registry):
            value = row.get("ppsd_npz_file", "")
            if not value or row.get("ppsd_status") not in {"complete", "partial"}:
                continue
            registered[str(Path(value).expanduser().resolve())].append(
                {"registry": registry, "row": row}
            )
    candidates = set(registered)
    for root in config["product_roots"]:
        if Path(root).is_dir():
            candidates.update(
                str(p.with_name("ppsd.npz").resolve())
                for p in Path(root).rglob("processing.json")
            )
    return [
        describe_product(p, registered.get(p, []), library) for p in sorted(candidates)
    ]


def describe_product(npz, registry_rows, library):
    npz = Path(npz)
    product = {
        "npz": str(npz),
        "registered_in": sorted({r["registry"] for r in registry_rows}),
        "registry_keys": sorted(
            {tuple(r["row"].get(k, "") for k in KEYS) for r in registry_rows}
        ),
        "issues": [],
    }
    product["_registry_rows"] = registry_rows
    manifest = npz.with_name("processing.json")
    if not npz.is_file():
        product["issues"].append("missing_npz")
    if not manifest.is_file():
        product["issues"].append("no_provenance_manifest")
        return product
    data = json.loads(manifest.read_text())
    config, summary = data.get("config", {}), data.get("summary", {})
    product.update(
        provenance=str(manifest),
        provenance_sha256=digest(manifest),
        npz_sha256=digest(npz) if npz.is_file() else "",
        seed_id=config.get("seed_id", ""),
        epoch_start=config.get("epoch_start", ""),
        epoch_end=config.get("epoch_end", ""),
        sample_rate_hz=config.get("sample_rate_hz"),
        parameters={k: config.get(v) for k, v in PRODUCT_PARAMETERS.items()},
        fixed={k: config.get(k) for k in FIXED_PRODUCT_CONFIG},
        special_handling=config.get("special_handling"),
        psd_units=config.get("psd_units", summary.get("ppsd_units", "")),
        inventory_sha256=config.get("inventory_sha256", ""),
        source_inventory=data.get("source_inventory", ""),
        obspy_version=config.get("obspy_version", ""),
        chunk_strategy=config.get("chunk_strategy", ""),
        chunk_seconds=config.get("chunk_seconds"),
        status=summary.get("ppsd_status", "unknown"),
        segments=summary.get("ppsd_segments"),
        processed=[summary.get("ppsd_start_time"), summary.get("ppsd_end_time")],
        requested=[
            summary.get("ppsd_requested_start_time"),
            summary.get("ppsd_requested_end_time"),
        ],
        errors=dict(data.get("errors", {})),
        interrupted="interrupted" in str(summary.get("ppsd_error") or "").lower(),
        completed_chunks=data.get("completed_chunks", {}),
        out_of_range_values=summary.get("ppsd_histogram_out_of_range_values"),
    )
    sources = []
    for item in data.get("source_files", []):
        real = os.path.realpath(item["path"])
        state = "missing"
        if Path(real).is_file():
            stat = Path(real).stat()
            unchanged = (stat.st_size, stat.st_mtime_ns) == (
                item["size"],
                item["mtime_ns"],
            )
            state = "unchanged" if unchanged else "modified"
        sources.append({"path": item["path"], "realpath": real, "state": state})
    product["_sources"] = sources
    product["source_counts"] = dict(Counter(s["state"] for s in sources))
    snapshot = npz.with_name("inventory.xml")
    payload, reason = (None, "missing_snapshot")
    if snapshot.is_file():
        payload, reason = library.signature(
            snapshot, product["seed_id"], product["epoch_start"], product["epoch_end"]
        )
        product["snapshot_sha256"] = digest(snapshot)
    product["_response"] = payload
    if reason:
        product["issues"].append(f"snapshot_response: {reason}")
    return product


def product_coverage(product, parameters):
    """Intervals the product's PSD windows could have used, minus error chunks."""
    start, end = product["requested"]
    if not start or not end:
        return []
    covered = [(timestamp(start), timestamp(end))]
    length = parameters["ppsd_length"]
    step = length * (1 - parameters["overlap"])
    chunk = product.get("chunk_seconds") or parameters["chunk_seconds"]
    failed = [
        (timestamp(key), timestamp(key) + chunk + length - step)
        for key in product.get("errors", {})
    ]
    if product.get("interrupted") and product.get("completed_chunks"):
        resume = max(
            timestamp(v["next_start"]) for v in product["completed_chunks"].values()
        )
        failed.append((resume, timestamp(end)))
    return subtract_intervals(covered, failed)


# --------------------------------------------------------------------------
# Historical evidence for the six checks (existing six-check review loaders)


def load_evidence(config, networks):
    """Join saved evidence with the reviewed exact-key loaders; never by name."""
    inputs = quality_review.Inputs(config["source_root"])
    networks = set(networks)
    errors, unjoined = [], []
    try:
        records = quality_review.load_masters(
            inputs, config["evidence"]["masters"], networks
        )
    except Exception as exc:
        records = {}
        errors.append(f"masters: {type(exc).__name__}: {exc}")
    items = {
        key: quality_review.base_review(key, row)
        for key, row in records.items()
        if row.get("ppsd_npz_file")
    }
    for kind in ("microseism", "local_p"):
        if kind not in config["evidence"]:
            continue
        try:
            quality_review.attach_staging(
                inputs,
                items,
                records,
                config["evidence"][kind],
                kind,
                networks,
                unjoined,
            )
        except Exception as exc:
            errors.append(f"{kind}: {type(exc).__name__}: {exc}")
    if "rayleigh" in config["evidence"]:
        try:
            quality_review.attach_rayleigh(
                inputs,
                items,
                records,
                config["evidence"]["rayleigh"],
                networks,
                unjoined,
            )
        except Exception as exc:
            errors.append(f"rayleigh: {type(exc).__name__}: {exc}")
    # Unjoined staging rows keep their own status, so "not run" stays visible.
    for entry in unjoined:
        kind = entry["check"]
        if kind in {"microseism", "local_p"} and kind in config["evidence"]:
            for row in read_registry(config["evidence"][kind]):
                if row["Seed_id"] == entry["seed_id"]:
                    entry["source_status"] = row.get(f"{kind}_status") or row.get(
                        "calibration_check_pressure_vertical_microseisms_status", ""
                    )
                    entry["note"] = row.get(f"{kind}_note", "")
    return items, unjoined, errors, inputs.hashes


# --------------------------------------------------------------------------
# Planning


def task_id(kind, payload):
    return f"{kind}-{json_hash(payload)[:16]}"


def key_dict(row):
    return {k: row[k] for k in KEYS}


def epoch_bounds(row):
    start = timestamp(row["channel_start_time"]) if row["channel_start_time"] else None
    end = timestamp(row["channel_end_time"]) if row["channel_end_time"] else None
    return start, end


def overlaps(a_start, a_end, b_start, b_end):
    low = max(a_start if a_start is not None else -math.inf, b_start or -math.inf)
    high = min(a_end if a_end is not None else math.inf, b_end or math.inf)
    return high > low


def build_plan(
    config,
    campaign_dir,
    *,
    scope="pressure",
    full_headers=False,
    progress=True,
):
    """Inventory everything and derive a deterministic missing-work task list."""
    if scope not in SCOPES:
        raise ValueError(f"Scope must be one of {SCOPES}")
    campaign_dir = Path(campaign_dir).resolve()
    parameters = config["ppsd_parameters"]
    tolerance = config["gap_tolerance_seconds"]
    cache = HeaderCache(campaign_dir / "cache/headers.json")
    files, ignored = scan_archive(
        config, cache, full_headers=full_headers, progress=progress
    )
    selected, duplicates = select_copies(files, config["layout_preference"], tolerance)
    index = {os.path.realpath(f["path"]): f for f in files if not f.get("issue")}
    library = ResponseLibrary()
    epochs, inventory_sources = read_epochs(config, library, set(selected))
    primary = assign_primary(epochs)
    by_seed = defaultdict(list)
    for row in primary:
        by_seed[row["Seed_id"]].append(row)
    products = discover_products(config, library)
    products_by_seed = defaultdict(list)
    for product in products:
        if product.get("seed_id"):
            products_by_seed[product["seed_id"]].append(product)
    evidence, unjoined, evidence_errors, evidence_hashes = load_evidence(
        config, config["networks"]
    )

    channels, tasks, unmatched_streams = [], [], []
    for seed, records in sorted(selected.items()):
        spans = merge_intervals(
            [span for r in records for span in file_intervals(r)], tolerance
        )
        rates = sorted({r["sample_rate_hz"] for r in records})
        matched = [
            row
            for row in by_seed.get(seed, [])
            if any(overlaps(*epoch_bounds(row), s, e) for s, e in spans)
        ]
        covered_by_epochs = []
        for row in matched:
            covered_by_epochs.extend(clip_intervals(spans, *epoch_bounds(row)))
        outside = subtract_intervals(spans, covered_by_epochs)
        if not matched or outside:
            unmatched_streams.append(
                {
                    "seed_id": seed,
                    "reason": "no_inventory_epoch"
                    if not matched
                    else "data_outside_inventory_epochs",
                    "seconds": total_seconds(outside if matched else spans),
                    "intervals": interval_text(outside if matched else spans),
                    "files": len(records),
                }
            )
        for row in matched:
            channels.append(
                plan_epoch(
                    row, records, spans, rates, products_by_seed[seed], config, index
                )
            )

    pressure = [c for c in channels if c["channel_class"] == "pressure"]
    for channel in channels:
        in_scope = channel["channel_class"] == "pressure" or (
            scope == "all-channels"
            and channel["channel_class"] in {"vertical", "horizontal", "other_seismic"}
        )
        channel["ppsd_in_scope"] = in_scope
        new = plan_ppsd(channel, campaign_dir, parameters)
        if in_scope:
            tasks.extend(new)
        elif channel["ppsd_decision"] in {"compute", "supplement"}:
            channel["ppsd_decision"] = "not_planned_out_of_scope"
            channel["ppsd_reason"] = (
                "pressure is the primary deliverable; use --ppsd-scope all-channels"
            )
            channel["ppsd_tasks"] = []
    for channel in channels:
        if channel["ppsd_tasks"]:
            problems = verify_headers(channel, cache)
            if problems:
                dropped = set(channel["ppsd_tasks"])
                tasks = [t for t in tasks if t["task_id"] not in dropped]
                partial = bool(channel.get("ppsd_product"))
                channel["ppsd_decision"] = (
                    "reuse_existing_partial" if partial else "blocked"
                )
                channel["ppsd_reason"] = (
                    "full-header verification: "
                    + "; ".join(problems[:5])
                    + (f" (+{len(problems) - 5} more)" if len(problems) > 5 else "")
                )
                channel["ppsd_tasks"] = []
                channel["ppsd_target"] = channel.get("ppsd_product")
            else:
                channel["full_header_verified"] = True
    cache.save()
    nearby_dependencies = [
        t["task_id"]
        for t in tasks
        if t["kind"] in {"ppsd", "ppsd_merge"} and t["channel_class"] == "pressure"
    ]
    for channel in pressure:
        tasks.extend(plan_spectral_checks(channel, campaign_dir, nearby_dependencies))
        attach_calibration_evidence(channel, evidence, unjoined, library)

    note_colocated_evidence(pressure)
    with_data = {tuple(c[k] for k in KEYS) for c in channels}
    inventory_only = [
        {
            k: row[k]
            for k in (
                "Seed_id",
                "channel_start_time",
                "channel_end_time",
                "inventory_file",
                "channel_class",
            )
        }
        for row in primary
        if tuple(row[k] for k in KEYS) not in with_data
    ]
    plan = {
        "_files": files,
        "plan_version": PLAN_VERSION,
        "created_utc": datetime.now(UTC).isoformat(),
        "scope": scope,
        "full_headers": full_headers,
        "config": config,
        "campaign_dir": str(campaign_dir),
        "header_cache": {"hits": cache.hits, "misses": cache.misses},
        "archive": {
            "data_directory": config["data_directory"],
            "recognized_files": len(files),
            "ignored_files": ignored,
            "files_with_issues": [
                {k: f.get(k) for k in ("path", "filename_seed_id", "seed_id", "issue")}
                for f in files
                if f.get("issue")
            ],
            "excluded_duplicate_copies": [
                {
                    k: f.get(k)
                    for k in ("path", "layout", "seed_id", "start", "end", "exclusion")
                }
                for f in duplicates
            ],
            "layouts": dict(Counter(f["layout"] for f in files)),
            "streams": len(selected),
            "unmatched_streams": unmatched_streams,
        },
        "inventory_sources": inventory_sources,
        "inventory_only_epochs": inventory_only,
        "products": products,
        "evidence_errors": evidence_errors,
        "evidence_unjoined": unjoined,
        "evidence_source_hashes": evidence_hashes,
        "channels": channels,
        "tasks": tasks,
    }
    plan["summary"] = summarize(plan)
    plan["plan_id"] = json_hash(
        [plan["scope"], [{k: t[k] for k in ("task_id", "dependencies")} for t in tasks]]
    )[:16]
    return plan


def public(value):
    """Drop private planning fields (leading underscore) before saving."""
    if isinstance(value, dict):
        return {k: public(v) for k, v in value.items() if not str(k).startswith("_")}
    if isinstance(value, list):
        return [public(v) for v in value]
    return value


def plan_epoch(row, records, spans, rates, products, config, index):
    """Coverage of one primary epoch and its relation to saved products."""
    parameters = config["ppsd_parameters"]
    start, end = epoch_bounds(row)
    available = clip_intervals(spans, start, end)
    epoch_files = [
        r
        for r in records
        if any(overlaps(start, end, a, b) for a, b in file_intervals(r))
    ]
    issues = []
    if row.get("response_status") != "usable":
        issues.append(
            f"response_unusable: {row.get('response_error', 'not evaluated')}"
        )
    expected = float(row["sample_rate_hz"] or 0)
    ppm = max(
        (abs(rate / expected - 1) * 1e6 for rate in rates if expected),
        default=math.inf,
    )
    if not expected or ppm > 1e-3:  # process_epoch requires math.isclose(rel 1e-9)
        issues.append(
            f"header_sample_rate_differs: StationXML {row['sample_rate_hz']} Hz, "
            f"headers {rates} Hz ({ppm:.4g} ppm); existing PPSD code refuses this"
        )
    if expected and parameters["minimum_frequency"] >= expected / 2:
        issues.append("minimum_frequency_not_below_nyquist")
    channel = {
        **{k: v for k, v in row.items() if not k.startswith("_")},
        "network": row["Seed_id"].split(".")[0],
        "station": row["Seed_id"].split(".")[1],
        "channel": row["Seed_id"].split(".")[3],
        "available_intervals": interval_text(available),
        "available_seconds": total_seconds(available),
        "available_start": iso(available[0][0]) if available else "",
        "available_end": iso(available[-1][1]) if available else "",
        "files": len(epoch_files),
        "file_layouts": dict(Counter(r["layout"] for r in epoch_files)),
        "header_sample_rates_hz": rates,
        "prerequisite_issues": issues,
        "_epoch": row,
        "_files": epoch_files,
        "_available": available,
    }
    channel["products"] = [
        relate_product(p, row, available, epoch_files, config, index)
        for p in products
        if overlaps(
            start,
            end,
            *(
                (timestamp(p["epoch_start"]) if p.get("epoch_start") else None),
                (timestamp(p["epoch_end"]) if p.get("epoch_end") else None),
            ),
        )
    ]
    return channel


def relate_product(product, row, available, epoch_files, config, index):
    """Compatibility of one saved product with this epoch, explicitly listed."""
    parameters = config["ppsd_parameters"]
    notes, blockers = [], []
    exact = (product.get("epoch_start"), product.get("epoch_end")) == (
        row["channel_start_time"],
        row["channel_end_time"],
    )
    if not exact:
        blockers.append("epoch_bounds_differ_requires_explicit_window_selection")
    response = same_response(row.get("_response"), product.get("_response"))
    if response is False:
        blockers.append("response_differs_from_primary_inventory")
    elif response is None:
        blockers.append("response_identity_unverified")
    if product.get("inventory_sha256") != row.get("inventory_sha256"):
        notes.append("built_from_other_inventory_file")
    differing = [
        k
        for k, v in product.get("parameters", {}).items()
        if v is None or not math.isclose(float(v), parameters[k])
    ]
    differing += [
        k for k, v in product.get("fixed", {}).items() if v != FIXED_PRODUCT_CONFIG[k]
    ]
    if differing:
        blockers.append("parameters_differ: " + ",".join(differing))
    if product.get("status") not in {"complete", "partial"}:
        blockers.append(f"product_status: {product.get('status')}")
    if product.get("obspy_version") and product["obspy_version"] != _obspy_version():
        notes.append(f"obspy_version: {product['obspy_version']}")
    if product.get("source_counts", {}).get("modified"):
        blockers.append("source_files_modified_since_processing")
    if product.get("source_counts", {}).get("missing"):
        notes.append("recorded_source_files_missing_now")
    for issue in product.get("issues", []):
        (notes if issue.startswith("snapshot") else blockers).append(issue)
    covered = product_coverage(product, parameters)
    used = {s["realpath"] for s in product.get("_sources", [])}
    new_spans = [
        span
        for record in epoch_files
        if os.path.realpath(record["path"]) not in used
        for span in file_intervals(record)
    ]
    # Files the product read, possibly another copy of the selected files.
    used_spans = merge_intervals(
        [span for path in used if path in index for span in file_intervals(index[path])]
    )
    unknown = sum(path not in index for path in used)
    if unknown:
        notes.append(f"recorded_sources_not_in_current_archive: {unknown}")
    uncovered = merge_intervals(
        subtract_intervals(available, covered)
        + subtract_intervals(
            clip_intervals(merge_intervals(new_spans), *epoch_bounds(row)), used_spans
        )
    )
    windowed = [(a, b) for a, b in uncovered if b - a >= parameters["ppsd_length"]]
    return {
        "npz": product["npz"],
        "npz_sha256": product.get("npz_sha256", ""),
        "registered": bool(product["registered_in"]),
        "registered_in": product["registered_in"],
        "status": product.get("status"),
        "segments": product.get("segments"),
        "processed": product.get("processed"),
        "requested": product.get("requested"),
        "exact_epoch": exact,
        "same_response": response,
        "compatible": not blockers,
        "blockers": blockers,
        "notes": notes,
        "uncovered_intervals": interval_text(windowed),
        "uncovered_seconds": total_seconds(windowed),
        "sub_window_remainder_seconds": total_seconds(uncovered)
        - total_seconds(windowed),
        "_product": product,
        "_uncovered": windowed,
    }


def _obspy_version():
    import obspy

    return obspy.__version__


def verify_headers(channel, cache):
    """Every record header of files a new PPSD would read: identity and rate."""
    problems = []
    expected = float(channel["sample_rate_hz"])
    for record in channel["_files"]:
        info, _ = cache.get(Path(record["path"]), full=True)
        name = Path(record["path"]).name
        if info.get("error"):
            problems.append(f"{name}: {info['error']}")
        elif info["seed_id"] != record["seed_id"]:
            problems.append(f"{name}: header {info['seed_id']}")
        elif not math.isclose(info["sample_rate_hz"], expected, rel_tol=1e-9):
            problems.append(f"{name}: {info['sample_rate_hz']} Hz")
    return problems


def ppsd_row(channel):
    return {k: channel[k] for k in KEYS} | {"sample_rate_hz": channel["sample_rate_hz"]}


def file_signature(records):
    return [[os.path.realpath(r["path"]), r["size"], r["mtime_ns"]] for r in records]


def plan_ppsd(channel, campaign_dir, parameters):
    """Reuse a compatible registered product; compute only what is missing."""
    tasks = []
    products = channel["products"]
    compatible = [p for p in products if p["compatible"] and p["registered"]]
    channel["ppsd_product"] = None
    channel["ppsd_tasks"] = []
    if not channel["_available"]:
        channel["ppsd_decision"] = "no_data_in_epoch"
        channel["ppsd_reason"] = "waveform files overlap only outside this epoch"
        return tasks
    if len(compatible) > 1:
        channel["ppsd_decision"] = "review_required"
        channel["ppsd_reason"] = "several compatible registered products: " + ", ".join(
            p["npz"] for p in compatible
        )
        return tasks
    blocked = "; ".join(channel["prerequisite_issues"])
    if not compatible and blocked:
        channel["ppsd_decision"] = "blocked"
        channel["ppsd_reason"] = blocked
        return tasks
    if not compatible and products:
        channel["ppsd_decision"] = "review_required"
        channel["ppsd_reason"] = "; ".join(
            f"{p['npz']}: {', '.join(p['blockers']) or 'unregistered product'}"
            for p in products
        )
        return tasks
    base = {
        "kind": "ppsd",
        "check": "ppsd",
        "channel_class": channel["channel_class"],
        "key": key_dict(channel),
        "seed_id": channel["Seed_id"],
        "parameters": dict(parameters),
        "prerequisites": [
            "usable full response in primary StationXML epoch",
            "header sample rate equal to StationXML",
            "header-named read-only symlink view of selected files",
        ],
        "inputs": {
            "inventory_sha256": channel["inventory_sha256"],
            "files": file_signature(channel["_files"]),
        },
        "dependencies": [],
    }
    if not compatible:
        interval = [channel["available_start"], channel["available_end"]]
        task = base | {"intended_coverage": [interval], "interval": interval}
        task["task_id"] = task_id("ppsd", task)
        task.update(
            output_path=str(campaign_dir / "products" / task["task_id"]),
            reuse_decision="compute",
            reason="no saved PPSD product for this channel epoch",
        )
        channel["ppsd_decision"] = "compute"
        channel["ppsd_reason"] = task["reason"]
        channel["ppsd_tasks"] = [task["task_id"]]
        channel["ppsd_target"] = task["task_id"]
        return [task]
    product = compatible[0]
    channel["ppsd_product"] = product["npz"]
    channel["ppsd_target"] = product["npz"]
    if not product["_uncovered"]:
        channel["ppsd_decision"] = "reuse_existing"
        channel["ppsd_reason"] = "registered product covers all available intervals"
        return tasks
    if blocked:
        channel["ppsd_decision"] = "reuse_existing_partial"
        channel["ppsd_reason"] = (
            f"{product['uncovered_seconds'] / 86400:.2f} uncovered days; "
            f"supplement blocked: {blocked}"
        )
        return tasks
    supplements = []
    for start, end in product["_uncovered"]:
        interval = [iso(start), iso(end)]
        task = base | {
            "kind": "ppsd",
            "intended_coverage": [interval],
            "interval": interval,
            "role": "supplement",
            "base_product": product["npz"],
            "base_product_sha256": product["npz_sha256"],
        }
        task["task_id"] = task_id("ppsd", task)
        task.update(
            output_path=str(campaign_dir / "products" / task["task_id"]),
            reuse_decision="supplement",
            reason="available data outside the saved product's processed interval",
        )
        supplements.append(task)
    merge = {
        "kind": "ppsd_merge",
        "check": "ppsd",
        "channel_class": channel["channel_class"],
        "key": key_dict(channel),
        "seed_id": channel["Seed_id"],
        "parameters": dict(parameters),
        "prerequisites": ["completed supplement PPSDs", "unchanged base product hash"],
        "base_product": product["npz"],
        "base_product_sha256": product["npz_sha256"],
        "dependencies": [t["task_id"] for t in supplements],
        "intended_coverage": channel["available_intervals"],
    }
    merge["task_id"] = task_id("ppsd_merge", merge)
    merge.update(
        output_path=str(campaign_dir / "products" / merge["task_id"]),
        reuse_decision="derive",
        reason="combine preserved base windows with supplement windows, "
        "without recomputing the base",
    )
    channel["ppsd_decision"] = "supplement"
    channel["ppsd_reason"] = (
        f"{len(supplements)} uncovered interval(s), "
        f"{product['uncovered_seconds'] / 86400:.2f} days"
    )
    channel["ppsd_tasks"] = [t["task_id"] for t in supplements] + [merge["task_id"]]
    channel["ppsd_target"] = merge["task_id"]
    return supplements + [merge]


def historical_result(channel, prefix):
    """A saved check result bound to the reused product by registry and hash."""
    product = next(
        (p for p in channel["products"] if p["npz"] == channel.get("ppsd_product")),
        None,
    )
    if product is None:
        return None
    for entry in product["_product"].get("_registry_rows", []):
        row = entry["row"]
        if tuple(row.get(k, "") for k in KEYS) != tuple(channel[k] for k in KEYS):
            continue
        path = row.get(f"{prefix}_result_file", "")
        status = row.get(f"{prefix}_status", "")
        if not path:
            if status:
                return {"status": status, "registry": entry["registry"], "path": ""}
            continue
        result = {"status": status, "registry": entry["registry"], "path": path}
        if not Path(path).is_file():
            return result | {"problem": "missing_saved_result"}
        data = json.loads(Path(path).read_text())
        if data.get("seed_id") not in {None, channel["Seed_id"]}:
            return result | {"problem": "result_seed_mismatch"}
        saved = data.get("source_npz_sha256")
        if saved and saved != product["npz_sha256"]:
            return result | {"problem": "result_bound_to_other_npz_content"}
        return result | {
            "status": data.get("status", status),
            "sha256": digest(path),
            "binding": "registry_row_and_npz_hash" if saved else "registry_row",
            "reference_warnings": quality_review.pressure_reference_warnings(data)
            if prefix == "pressure_expected"
            else [],
        }
    return None


def plan_spectral_checks(channel, campaign_dir, nearby_dependencies):
    tasks = []
    channel["checks"] = {}
    decision = channel.get("ppsd_decision")
    for check, (stage, prefix) in SPECTRAL_CHECKS.items():
        entry = {"check": check, "stage": stage}
        missing_metadata = [
            name
            for name in (
                ("Profondeur",)
                if check == "infragravity"
                else ("Lat", "Lon", "Profondeur")
                if check == "nearby"
                else ()
            )
            if not _positive_or_coordinate(channel.get(name, ""), name)
        ]
        if decision not in {
            "reuse_existing",
            "reuse_existing_partial",
            "compute",
            "supplement",
        }:
            entry.update(decision="blocked", reason=f"no usable PPSD ({decision})")
        elif missing_metadata:
            entry.update(
                decision="blocked",
                reason="missing/invalid metadata: " + ",".join(missing_metadata),
            )
        else:
            saved = (
                historical_result(channel, prefix) if decision != "compute" else None
            )
            usable = (
                saved
                and saved.get("path")
                and not saved.get("problem")
                and saved.get("status") not in {"failed", ""}
            )
            if usable and decision in {"reuse_existing", "reuse_existing_partial"}:
                entry.update(
                    decision="reuse_existing",
                    reason="saved result bound to the reused product",
                    result=saved,
                )
                if saved.get("reference_warnings"):
                    entry["reason"] += "; historical reference protocol: " + ",".join(
                        saved["reference_warnings"]
                    )
            else:
                target = channel["ppsd_target"]
                product_hash = next(
                    (
                        p["npz_sha256"]
                        for p in channel["products"]
                        if p["npz"] == target
                    ),
                    "",
                )
                task = {
                    "kind": check,
                    "check": check,
                    "stage": stage,
                    "channel_class": "pressure",
                    "key": key_dict(channel),
                    "seed_id": channel["Seed_id"],
                    "ppsd_source": target,
                    "ppsd_source_sha256": product_hash,
                    "parameters": {
                        "estimator": f"compy_pressure.network.diagnostic_job:{stage}"
                    },
                    "prerequisites": ["pressure PPSD"]
                    + (
                        [
                            "every pressure PPSD task of this plan attempted "
                            "(reference pool; failures are excluded and recorded)"
                        ]
                        if check == "nearby"
                        else []
                    ),
                    "dependencies": [target] if target in channel["ppsd_tasks"] else [],
                    "intended_coverage": channel["available_intervals"],
                }
                if check == "nearby":
                    task["pool_dependencies"] = [
                        t for t in nearby_dependencies if t not in task["dependencies"]
                    ]
                task["task_id"] = task_id(check, task)
                if saved and decision in {"reuse_existing", "reuse_existing_partial"}:
                    problem = saved.get("problem") or saved["status"]
                    reason = f"saved result unusable ({problem})"
                elif saved:
                    reason = "saved result covers only the base product interval"
                else:
                    reason = "no saved result for this product"
                task.update(
                    output_path=str(
                        campaign_dir
                        / (
                            "views"
                            if target == channel.get("ppsd_product")
                            else "products"
                        )
                    ),
                    reuse_decision="compute",
                    reason=reason,
                    historical=saved,
                )
                tasks.append(task)
                entry.update(
                    decision="compute",
                    reason=reason,
                    task_id=task["task_id"],
                    historical=saved,
                )
        channel["checks"][check] = entry
    return tasks


def _positive_or_coordinate(value, name):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    if not math.isfinite(number):
        return False
    return number > 0 if name == "Profondeur" else True


def attach_calibration_evidence(channel, evidence, unjoined, library):
    """Historical calibration evidence by exact or response-equivalent key."""
    seed = channel["Seed_id"]
    start, end = epoch_bounds(channel)
    own = (
        seed,
        channel["channel_start_time"],
        channel["channel_end_time"],
        channel["inventory_file"],
    )
    matches = []
    for key, item in evidence.items():
        if key[0] != seed:
            continue
        other_start = timestamp(key[1]) if key[1] else None
        other_end = timestamp(key[2]) if key[2] else None
        if not overlaps(start, end, other_start, other_end):
            continue
        matches.append((key, item, evidence_relation(channel, key, own, library)))
    for check in CALIBRATION_CHECKS:
        found = [
            {
                "key": list(key),
                "relation": relation,
                "state": item["checks"][check]["state"],
                "source_status": item["checks"][check]["source_status"],
                "flags": item["checks"][check]["flags"],
                "binding": item["checks"][check]["binding"],
                "sources": [s.get("path") for s in item["checks"][check]["sources"]],
            }
            for key, item, relation in matches
            if item["checks"][check]["binding"] != "none"
        ]
        notes = [u for u in unjoined if u["check"] == check and u["seed_id"] == seed]
        usable = [f for f in found if f["relation"] != "response_differs"]
        if usable:
            states = sorted({f["state"] for f in usable})
            decision = "reuse_historical"
            reason = "historical evidence: " + ",".join(states)
            if check == "local_p" and "legacy_unvalidated" in states:
                reason += "; legacy local-P protocol lacks the v2 noise-window guard"
        elif notes and all(
            str(n.get("source_status", "")).startswith("not_applicable") for n in notes
        ):
            decision = "not_applicable"
            reason = "historical staging row: " + "; ".join(
                f"{n['source_status']} ({n.get('note', '')})" for n in notes
            )
        else:
            decision = "deferred"
            reason = (
                "no joined historical result; a new calibration run needs an "
                "explicitly frozen protocol and is not executed by this workflow"
            )
            if found:
                reason += "; only identities with a different response were found"
            if notes:
                reason += "; unjoined historical rows: " + ",".join(
                    f"{n['reason']}/{n.get('source_status', '')}" for n in notes
                )
        channel["checks"][check] = {
            "check": check,
            "decision": decision,
            "reason": reason,
            "evidence": found,
            "unjoined": notes,
        }


def note_colocated_evidence(pressure):
    """Deferred streams at stations whose other pressure stream has evidence.

    Metadata do not establish that two streams share one physical sensor, so
    this is a pointer for review, never a transfer of the other result.
    """
    by_station = defaultdict(list)
    for channel in pressure:
        by_station[(channel["network"], channel["station"])].append(channel)
    for channels in by_station.values():
        for channel in channels:
            for check in CALIBRATION_CHECKS:
                entry = channel["checks"][check]
                if entry["decision"] != "deferred":
                    continue
                others = sorted(
                    other["Seed_id"]
                    for other in channels
                    if other is not channel
                    and other["checks"][check]["decision"] == "reuse_historical"
                )
                entry["colocated_evidence"] = others
                if others:
                    entry["reason"] += (
                        "; station has historical evidence on "
                        + ",".join(others)
                        + " (not transferred)"
                    )


def evidence_relation(channel, key, own, library):
    """How a saved evidence identity relates to this primary epoch."""
    if key == own:
        return "exact_key"
    span = "same_epoch" if key[1:3] == own[1:3] else "overlapping_epoch"
    payload, _ = library.signature(key[3], *key[:3])
    same = same_response(channel.get("_epoch", {}).get("_response"), payload)
    if same is None:
        return f"{span}_response_unverifiable"
    return f"{span}_same_response" if same else "response_differs"


def summarize(plan):
    channels = plan["channels"]
    tasks = plan["tasks"]
    count = Counter
    return {
        "networks": sorted({c["network"] for c in channels}),
        "stations": len({(c["network"], c["station"]) for c in channels}),
        "channel_epochs_with_data": len(channels),
        "channel_epochs_by_class": dict(count(c["channel_class"] for c in channels)),
        "pressure_epochs": sum(c["channel_class"] == "pressure" for c in channels),
        "ppsd_decisions": dict(
            count(f"{c['channel_class']}:{c['ppsd_decision']}" for c in channels)
        ),
        "check_decisions": dict(
            count(
                f"{name}:{entry['decision']}"
                for c in channels
                for name, entry in c.get("checks", {}).items()
            )
        ),
        "tasks_by_kind": dict(count(t["kind"] for t in tasks)),
        "planned_ppsd_days": sum(
            (timestamp(t["interval"][1]) - timestamp(t["interval"][0])) / 86400
            for t in tasks
            if t["kind"] == "ppsd"
        ),
        "inventory_only_epochs": len(plan["inventory_only_epochs"]),
        "unmatched_streams": len(plan["archive"]["unmatched_streams"]),
        "files_with_issues": len(plan["archive"]["files_with_issues"]),
        "excluded_duplicate_copies": len(plan["archive"]["excluded_duplicate_copies"]),
        "products_seen": len(plan["products"]),
        "evidence_errors": len(plan["evidence_errors"]),
    }


# --------------------------------------------------------------------------
# Output


def write_csv(path, rows, fields):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _cell(row.get(k, "")) for k in fields})


def _cell(value):
    if isinstance(value, (list, dict)):
        return json.dumps(value, sort_keys=True)
    return value


CHANNEL_FIELDS = [
    "network",
    "station",
    "channel",
    "Seed_id",
    "channel_start_time",
    "channel_end_time",
    "inventory_file",
    "channel_class",
    "sample_rate_hz",
    "header_sample_rates_hz",
    "response_status",
    "response_signature",
    "alternative_sources",
    "available_start",
    "available_end",
    "available_seconds",
    "files",
    "file_layouts",
    "prerequisite_issues",
    "ppsd_in_scope",
    "ppsd_decision",
    "ppsd_reason",
    "ppsd_product",
    "ppsd_tasks",
]
TASK_FIELDS = [
    "task_id",
    "kind",
    "check",
    "seed_id",
    "key",
    "interval",
    "intended_coverage",
    "reuse_decision",
    "reason",
    "dependencies",
    "output_path",
    "parameters",
    "base_product",
    "ppsd_source",
]


def write_plan(plan, directory):
    """Write a fresh plan folder; existing plans are never replaced."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    atomic_json(directory / "plan.json", public(plan))
    write_csv(directory / "channels.csv", plan["channels"], CHANNEL_FIELDS)
    write_csv(directory / "tasks.csv", plan["tasks"], TASK_FIELDS)
    rows = []
    for c in plan["channels"]:
        for name, entry in c.get("checks", {}).items():
            rows.append(
                {
                    "Seed_id": c["Seed_id"],
                    "channel_start_time": c["channel_start_time"],
                    "channel_end_time": c["channel_end_time"],
                    "inventory_file": c["inventory_file"],
                    "check": name,
                    "decision": entry["decision"],
                    "reason": entry["reason"],
                    "task_id": entry.get("task_id", ""),
                    "historical_states": sorted(
                        {e["state"] for e in entry.get("evidence", [])}
                    ),
                    "historical_relations": sorted(
                        {e["relation"] for e in entry.get("evidence", [])}
                    ),
                }
            )
    write_csv(directory / "checks.csv", rows, list(rows[0]) if rows else ["check"])
    blocked = [
        {
            "Seed_id": c["Seed_id"],
            "channel_start_time": c["channel_start_time"],
            "item": "ppsd",
            "decision": c["ppsd_decision"],
            "reason": c["ppsd_reason"],
        }
        for c in plan["channels"]
        if c["ppsd_decision"] in {"blocked", "review_required"}
    ] + [
        {
            "Seed_id": r["Seed_id"],
            "channel_start_time": r["channel_start_time"],
            "item": r["check"],
            "decision": r["decision"],
            "reason": r["reason"],
        }
        for r in rows
        if r["decision"] in {"blocked", "deferred"}
    ]
    write_csv(
        directory / "blocked.csv",
        blocked,
        ["Seed_id", "channel_start_time", "item", "decision", "reason"],
    )
    (directory / "summary.txt").write_text(summary_text(plan))
    return directory


def summary_text(plan):
    s = plan["summary"]
    lines = [
        f"Plan {plan['plan_id']} ({plan['scope']} scope), "
        f"created {plan['created_utc']}",
        f"Archive: {plan['archive']['recognized_files']} files, "
        f"{plan['archive']['streams']} header streams; layouts "
        f"{plan['archive']['layouts']}",
        f"Stations {s['stations']}; channel epochs with data "
        f"{s['channel_epochs_with_data']} {s['channel_epochs_by_class']}",
        f"PPSD decisions: {s['ppsd_decisions']}",
        f"Check decisions: {s['check_decisions']}",
        f"Tasks: {s['tasks_by_kind']}; planned new PPSD days "
        f"{s['planned_ppsd_days']:.1f}",
        f"Unmatched streams {s['unmatched_streams']}; header issues "
        f"{s['files_with_issues']}; excluded duplicate copies "
        f"{s['excluded_duplicate_copies']}; inventory-only epochs "
        f"{s['inventory_only_epochs']}; evidence errors {s['evidence_errors']}",
        "No gain, score, threshold or approval is produced by this plan.",
    ]
    return "\n".join(lines) + "\n"


def default_campaign_dir(name):
    return ROOT / "results/COVERAGE" / name


def guard_campaign_dir(campaign_dir, config):
    """Campaign outputs must stay outside raw inputs and saved product roots."""
    campaign_dir = Path(campaign_dir).resolve()
    protected = [Path(config["data_directory"])]
    protected += [Path(p).resolve() for p in config["product_roots"]]
    protected += [Path(p).resolve().parent for p in config["product_registries"]]
    for patterns in config["inventory_sources"].values():
        protected += [Path(p).resolve().parent for p in expand_sources(patterns)]
    for path in protected:
        if campaign_dir == path or path in campaign_dir.parents:
            raise ValueError(f"Campaign folder must not be inside protected {path}")
    return campaign_dir


def add_common_arguments(parser):
    parser.add_argument("--campaign", default="full_data_coverage_v1")
    parser.add_argument("--campaign-dir", type=Path)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument(
        "--source-root",
        type=Path,
        default=ROOT,
        help="Checkout holding the saved results/exports named in --config",
    )
    parser.add_argument("--ppsd-scope", choices=SCOPES, default="pressure")
    parser.add_argument(
        "--full-headers",
        action="store_true",
        help="Read every MiniSEED record header (exact gaps; much slower)",
    )


def prepare(args):
    config = load_config(args.config, args.source_root)
    campaign_dir = guard_campaign_dir(
        args.campaign_dir or default_campaign_dir(args.campaign), config
    )
    return config, campaign_dir


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    plan_parser = commands.add_parser("plan", help="Read-only inventory and dry run")
    add_common_arguments(plan_parser)
    run_parser = commands.add_parser("run", help="Dry run, or --execute missing tasks")
    add_common_arguments(run_parser)
    run_parser.add_argument("--execute", action="store_true")
    run_parser.add_argument("--workers", type=int, default=4)
    run_parser.add_argument("--only-kind", action="append")
    run_parser.add_argument(
        "--only-seed",
        action="append",
        help="Bounded run: tasks of these SEED IDs only (dependencies must exist)",
    )
    report_parser = commands.add_parser("report", help="Coverage report and heat maps")
    add_common_arguments(report_parser)
    report_parser.add_argument("--output-dir", type=Path)
    report_parser.add_argument("--max-percentage", type=float, default=30.0)
    report_parser.add_argument(
        "--curated-dir",
        type=Path,
        help="Also write small shareable tables and detail heat maps here",
    )
    args = parser.parse_args(argv)
    config, campaign_dir = prepare(args)
    if args.command == "run":
        from compy_pressure.coverage_run import run_campaign

        result = run_campaign(
            config,
            campaign_dir,
            scope=args.ppsd_scope,
            full_headers=args.full_headers,
            execute=args.execute,
            workers=args.workers,
            only_kinds=args.only_kind,
            only_seeds=args.only_seed,
        )
        raise SystemExit(0 if result["ok"] else 1)
    plan = build_plan(
        config, campaign_dir, scope=args.ppsd_scope, full_headers=args.full_headers
    )
    if args.command == "report":
        from compy_pressure.coverage_report import write_report

        output = write_report(
            plan, campaign_dir, args.output_dir, max_percentage=args.max_percentage
        )
        print(f"Report: {output}", flush=True)
        if args.curated_dir:
            from compy_pressure.coverage_report import write_curated

            print(f"Curated: {write_curated(plan, output, args.curated_dir)}")
        return
    folder = write_plan(plan, campaign_dir / "plans" / fresh_name(plan["plan_id"]))
    print(summary_text(plan) + f"Saved {folder}", flush=True)


if __name__ == "__main__":
    main()
