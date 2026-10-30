"""Recompute everything from the raw archive into one fresh run folder.

One command, resumable, using only local copies of station metadata and
catalogues (``docs/from_scratch_config.json``). Stages, in order:

1. ``preflight``  committed code, Python/imports, every configured input exists
2. ``workspaces`` per-network calibration masters: YV and YL from their local
                  StationXML; 7D and XO through their preparation scripts fed the
                  saved responses/catalogues (``--offline``: never downloads);
                  then one standalone response snapshot per channel epoch
3. ``identity``   YV Rayleigh with bounded events; must return RR52
                  g = 0.6787278140712512 exactly, otherwise the run stops
4. ``coverage``   fresh coverage campaign with reuse disabled: PPSDs for all
                  channels with usable responses, plus the expected-pressure,
                  nearby and infragravity checks, for every network
5. ``rayleigh``   every network with both event policies (bounded and all)
6. ``microseism`` common-week pressure/vertical check on the fresh masters
7. ``local_p``    reviewed v2 local-P campaign, station table and export
8. ``registry``   master-style CSV of the new PPSDs, for the plot scripts
9. ``summary``    status counts and comparison with the frozen results

Everything is written below ``results/FROM_SCRATCH/run_<stamp>/`` (local P
must live below ``results/LOCAL_P/campaigns/``). Raw data, StationXML, masters,
frozen results and saved products are only read. No gain is applied and no
threshold, score or rating is set.

Examples::

    python scripts/run_from_scratch.py --check
    python scripts/run_from_scratch.py --workers 6
    python scripts/run_from_scratch.py --run-root results/FROM_SCRATCH/run_X  # resume
    python scripts/run_from_scratch.py --networks YV --stations RR52   # smoke test
"""

import argparse
import copy
import csv
import hashlib
import json
import os
import shutil
import statistics
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "docs/from_scratch_config.json"
SCHEMA = "bruitfm-from-scratch-v1"
STAGES = (
    "preflight",
    "workspaces",
    "identity",
    "coverage",
    "rayleigh",
    "microseism",
    "local_p",
    "registry",
    "summary",
)
EPOCH_KEYS = ("Seed_id", "channel_start_time", "channel_end_time", "inventory_file")
SNAPSHOT_COLUMN = "rayleigh_response_snapshot_file"
# Masters are written once (the project refuses to overwrite one); the
# calibration master already carries every standalone snapshot path.
CALIBRATION_MASTER = "_calibration_master"


class IdentityFailure(RuntimeError):
    """The fresh YV reproduction did not return the frozen RR52 multiplier."""


def utc_now():
    return datetime.now(UTC).isoformat(timespec="seconds")


def sha256(path):
    with open(path, "rb") as stream:
        # hashlib.file_digest needs Python 3.11; ComPy supports 3.10
        checksum = hashlib.sha256()
        for block in iter(lambda: stream.read(1 << 20), b""):
            checksum.update(block)
        return checksum.hexdigest()


# --------------------------------------------------------------------------
# Configuration


def resolve(value, root=ROOT):
    path = Path(value).expanduser()
    return path if path.is_absolute() else root / path


def load_config(path=CONFIG, root=ROOT):
    config = json.loads(Path(path).read_text())
    if config.get("schema") != SCHEMA:
        raise ValueError(f"Unsupported from-scratch config schema in {path}")
    policies = config["event_policies"]
    if not policies or set(policies) - {"bounded", "all"}:
        raise ValueError("event_policies must be a subset of bounded, all")
    for name, network in config["networks"].items():
        if network["workspace"] not in {"inventory", "prepare"}:
            raise ValueError(f"{name}: workspace must be inventory or prepare")
        if network.get("frozen_policy") not in {None, "bounded", "all"}:
            raise ValueError(f"{name}: unknown frozen_policy")
    config["_path"] = str(Path(path).resolve())
    config["_sha256"] = sha256(path)
    config["_root"] = str(root)
    return config


def required_inputs(config, networks, root=ROOT):
    """Every local file or folder the selected networks need."""
    paths = [
        config["data_directory"],
        config["coverage_template"],
        config["dispersion_file"],
        *config["local_p_catalogs"],
    ]
    for name in networks:
        network = config["networks"][name]
        paths += [network["catalog"]] if "catalog" in network else []
        paths += network.get("inventories", [])
        paths += network.get("source_inventories", [])
        paths += network.get("local_responses", [])
        paths += [network["local_catalog"]] if "local_catalog" in network else []
        paths += [network["frozen_workflow"]] if network.get("frozen_workflow") else []
        paths += [network["prepare_script"]] if "prepare_script" in network else []
    return [resolve(p, root) for p in paths]


def data_dir(config, name, workspace):
    network = config["networks"][name]
    if network["workspace"] == "prepare":
        return workspace / "_inputs/waveforms"
    return resolve(config["data_directory"])


# --------------------------------------------------------------------------
# Code state and run record


def code_state(root=ROOT):
    """Hash of HEAD plus any uncommitted change in ../compy_pressure/ or scripts/."""
    parts = [
        subprocess.run(
            ["git", "-C", str(root), *args], capture_output=True, text=True
        ).stdout
        for args in (
            ["rev-parse", "HEAD"],
            [
                "status",
                "--porcelain",
                "--untracked-files=all",
                "--",
                "../compy_pressure",
                "scripts",
            ],
            ["diff", "HEAD", "--", "../compy_pressure", "scripts"],
        )
    ]
    return {
        "head": parts[0].strip(),
        "dirty": parts[1].strip(),
        "state": hashlib.sha256("".join(parts).encode()).hexdigest()[:16],
    }


class Run:
    """The run folder and its run.json stage record."""

    def __init__(self, root):
        self.root = Path(root).resolve()
        self.path = self.root / "run.json"
        self.data = json.loads(self.path.read_text()) if self.path.exists() else {}

    def save(self):
        self.root.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(".tmp")
        temporary.write_text(json.dumps(self.data, indent=2, default=str) + "\n")
        temporary.replace(self.path)

    def stage(self, name):
        return self.data.setdefault("stages", {}).get(name, {})

    def done(self, name):
        return self.stage(name).get("status") == "complete"

    def mark(self, name, status, **details):
        record = self.data.setdefault("stages", {}).setdefault(name, {})
        record.update(status=status, updated_utc=utc_now(), **details)
        self.save()


def run_command(command, log, env=None, cwd=ROOT):
    """Run one command, streaming output to the console and a log file."""
    log = Path(log)
    log.parent.mkdir(parents=True, exist_ok=True)
    printable = " ".join(str(c) for c in command)
    print(f"\n$ {printable}", flush=True)
    with log.open("a") as handle:
        handle.write(f"\n# {utc_now()}\n$ {printable}\n")
        process = subprocess.Popen(
            [str(c) for c in command],
            cwd=cwd,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        lines = []
        for line in process.stdout:
            print(line, end="", flush=True)
            handle.write(line)
            lines = (lines + [line])[-200:]
        status = process.wait()
        handle.write(f"# exit {status}\n")
    return status, lines


def child_env(root=ROOT):
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{root.parent}{os.pathsep}{root}"
    env.update(
        OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        MPLBACKEND="Agg",
    )
    return env


# --------------------------------------------------------------------------
# Calibration workspaces


def inventory_master(sources, links):
    """Fresh master rows from local StationXML; rows keep the real inventory path."""
    import pandas as pd

    from compy_pressure.inventory import build_master_dataset

    frames = []
    for index, source in enumerate(sources):
        source = Path(source).expanduser().resolve()
        if source.is_dir():
            frames.append(build_master_dataset(source))
            continue
        folder = links / f"{index:02d}"
        folder.mkdir(parents=True, exist_ok=True)
        link = folder / source.name
        if not link.exists():
            link.symlink_to(source)
        frames.append(build_master_dataset(folder))
    master = pd.concat(frames, ignore_index=True).astype(object)
    master["inventory_file"] = master["inventory_file"].map(
        lambda p: os.path.realpath(p) if p else p
    )
    return master.drop_duplicates(list(EPOCH_KEYS)).reset_index(drop=True)


def epoch_inventory(row, cache):
    """The single exact channel epoch of a master row (as ppsd.select_epoch)."""
    from obspy import read_inventory

    path = row["inventory_file"]
    if path not in cache:
        cache[path] = read_inventory(path, format="STATIONXML")
    network, station, location, code = row["Seed_id"].split(".")
    subset = copy.deepcopy(
        cache[path].select(
            network=network, station=station, location=location, channel=code
        )
    )
    channels = [
        channel
        for net in subset
        for sta in net
        for channel in sta
        if str(channel.start_date or "") == row["channel_start_time"]
        and str(channel.end_date or "") == row["channel_end_time"]
    ]
    if len(channels) != 1:
        raise ValueError(f"Expected one exact epoch for {row['Seed_id']}")
    for net in subset:
        for sta in net:
            sta.channels = [c for c in sta if c is channels[0]]
    return subset, channels[0]


def write_snapshots(master, workspace):
    """One standalone response snapshot per usable channel epoch.

    Calibrations need an exact saved response, not a PPSD. Existing snapshots
    (XO preparation) are kept. Returns the updated master and counts.
    """
    from compy_pressure.ppsd import response_handling

    master = master.astype(object).copy()
    if SNAPSHOT_COLUMN not in master:
        master[SNAPSHOT_COLUMN] = ""
    master[SNAPSHOT_COLUMN] = master[SNAPSHOT_COLUMN].fillna("")
    cache, counts = {}, Counter()
    for index, row in master.iterrows():
        if row[SNAPSHOT_COLUMN] and Path(row[SNAPSHOT_COLUMN]).is_file():
            counts["kept"] += 1
            continue
        try:
            inventory, channel = epoch_inventory(row, cache)
            response_handling(channel)
        except Exception as error:  # no usable response: no snapshot, no calibration
            counts[f"unusable:{type(error).__name__}"] += 1
            continue
        tag = hashlib.sha256(
            json.dumps([str(row[k]) for k in EPOCH_KEYS]).encode()
        ).hexdigest()[:16]
        destination = workspace / "_responses" / f"{row['Seed_id']}__{tag}"
        destination.mkdir(parents=True, exist_ok=True)
        inventory.write(str(destination / "inventory.xml"), format="STATIONXML")
        master.at[index, SNAPSHOT_COLUMN] = str(destination / "inventory.xml")
        counts["written"] += 1
    return master, dict(counts)


def seed_prepare_inputs(network, workspace, root=ROOT):
    """Copy the saved responses and catalogue where the preparation looks."""
    responses = workspace / "_inputs/responses"
    responses.mkdir(parents=True, exist_ok=True)
    for source in network["local_responses"]:
        source = resolve(source, root)
        target = responses / source.name
        if not target.exists():
            shutil.copy2(source, target)
        elif sha256(target) != sha256(source):
            raise ValueError(f"Seeded response differs from its source: {target}")
    catalog = workspace / "_inputs/catalog_usgs_m5p5.json"
    source = resolve(network["local_catalog"], root)
    if not catalog.exists():
        shutil.copy2(source, catalog)
    elif sha256(catalog) != sha256(source):
        raise ValueError(f"Seeded catalogue differs from its source: {catalog}")


# --------------------------------------------------------------------------
# Rayleigh


def rayleigh_command(config, name, policy, workspace, stations=(), python=None):
    network = config["networks"][name]
    common = config["rayleigh_common"]
    catalog = (
        workspace / "_inputs/catalog_usgs_m5p5.json"
        if network["workspace"] == "prepare"
        else resolve(network["catalog"])
    )
    command = [
        python or sys.executable,
        "-m",
        "compy_pressure.network",
        "rayleigh",
        "--network",
        name,
        "--data-dir",
        data_dir(config, name, workspace),
        "--master-dir",
        workspace / CALIBRATION_MASTER,
        "--results-dir",
        workspace / f"rayleigh_{policy}",
        "--rayleigh-catalog",
        catalog,
        "--rayleigh-dispersion-file",
        resolve(config["dispersion_file"]),
        "--rayleigh-frequency-min",
        common["band"][0],
        "--rayleigh-frequency-max",
        common["band"][1],
        "--rayleigh-pressure-sign",
        network["pressure_sign"],
        "--rayleigh-tilt-method",
        common["tilt_method"],
        "--rayleigh-min-events",
        common["min_events"],
        "--rayleigh-max-events",
        common["max_events"],
        "--rayleigh-event-sampling",
        policy,
        "--rayleigh-no-master-update",
    ]
    if network.get("event_resample_rate") is not None:
        command += ["--rayleigh-event-resample-rate", network["event_resample_rate"]]
    if network.get("native_rate_tolerance_ppm") is not None:
        command += [
            "--rayleigh-native-rate-tolerance-ppm",
            network["native_rate_tolerance_ppm"],
        ]
    for station in stations:
        command += ["--station", station]
    return command


def workflow_entries(path):
    entries = json.loads(Path(path).read_text())["entries"].values()
    return {e["seed_id"]: e for e in entries}


def multiplier(entry):
    value = entry["update"].get("rayleigh_pressure_multiplier")
    return None if value in (None, "") else str(value)


def identity_value(workflow, seed_id):
    entry = workflow_entries(workflow).get(seed_id)
    return multiplier(entry) if entry else None


def compare_workflows(new_path, frozen_path):
    """Status/gain agreement of a fresh workflow with the frozen one."""
    new, frozen = workflow_entries(new_path), workflow_entries(frozen_path)
    common = sorted(set(new) & set(frozen))
    identical, changes, differences = 0, [], []
    for seed in common:
        a, b = frozen[seed], new[seed]
        same = a["status"] == b["status"] and multiplier(a) == multiplier(b)
        identical += same
        if not same:
            differences.append(
                [seed, a["status"], multiplier(a), b["status"], multiplier(b)]
            )
        if multiplier(a) and multiplier(b):
            changes.append(float(multiplier(b)) / float(multiplier(a)) - 1)
    return {
        "stations_compared": len(common),
        "identical": identical,
        "only_frozen": sorted(set(frozen) - set(new)),
        "only_new": sorted(set(new) - set(frozen)),
        "differences": differences,
        "estimates_in_both": len(changes),
        "median_gain_change": statistics.median(changes) if changes else None,
        "max_abs_gain_change": max(map(abs, changes)) if changes else None,
    }


# --------------------------------------------------------------------------
# Coverage and registry


def coverage_config(template, networks, destination):
    """The coverage sources with every saved product removed: nothing reused."""
    config = json.loads(Path(template).read_text())
    config["networks"] = [n for n in config["networks"] if n in networks]
    config["product_registries"] = []
    config["product_roots"] = []
    config["evidence"] = {}
    config["from_scratch_note"] = "Reuse disabled: every PPSD and check is computed."
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(config, indent=2) + "\n")
    return destination


def registry_rows(campaign, masters=()):
    """Master-style rows for every completed PPSD of a coverage campaign."""
    metadata = {}
    for master in masters:
        with open(master, newline="") as stream:
            for row in csv.DictReader(stream):
                key = (
                    row["Seed_id"],
                    row["channel_start_time"],
                    row["channel_end_time"],
                )
                metadata.setdefault(key, row)
    products = {}
    for state_file in sorted(Path(campaign).glob("tasks/*/state.json")):
        state = json.loads(state_file.read_text())
        if state.get("kind") not in {"ppsd", "ppsd_merge"}:
            continue
        if state.get("status") != "complete":
            continue
        for output in state.get("outputs", []):
            if not output["path"].endswith("ppsd.npz"):
                continue
            npz = Path(output["path"])
            manifest = npz.with_name("processing.json")
            if not manifest.is_file():
                continue
            cfg = json.loads(manifest.read_text())["config"]
            key = (cfg["seed_id"], cfg["epoch_start"], cfg["epoch_end"])
            # A merged product supersedes its base PPSD for the same epoch.
            if key in products and state["kind"] != "ppsd_merge":
                continue
            row = dict(metadata.get(key, {}))
            row.update(
                Seed_id=cfg["seed_id"],
                channel_start_time=cfg["epoch_start"],
                channel_end_time=cfg["epoch_end"],
                inventory_file=json.loads(manifest.read_text()).get(
                    "source_inventory", row.get("inventory_file", "")
                ),
                ppsd_npz_file=str(npz),
                ppsd_status="complete",
                ppsd_units=cfg["psd_units"],
            )
            products[key] = row
    return [products[k] for k in sorted(products)]


def write_registry(rows, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with open(destination, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return destination


def task_counts(campaign):
    counts = Counter()
    for state_file in Path(campaign).glob("tasks/*/state.json"):
        state = json.loads(state_file.read_text())
        counts[f"{state.get('kind')}:{state.get('status')}"] += 1
    return dict(sorted(counts.items()))


# --------------------------------------------------------------------------
# Stages


class Pipeline:
    def __init__(self, args, config):
        self.args = args
        self.config = config
        self.networks = args.networks or list(config["networks"])
        unknown = set(self.networks) - set(config["networks"])
        if unknown:
            raise SystemExit(f"Unknown networks: {sorted(unknown)}")
        self.stations = args.stations or []
        self.env = child_env()
        self.python = sys.executable

    # ------------------------------------------------------------ helpers
    def workspace(self, name):
        return self.run.root / "calibration" / name

    def log(self, name):
        return self.run.root / "logs" / f"{name}.log"

    def same_code(self, stage):
        current = code_state()["state"]
        if current != self.run.data["code"]["state"]:
            raise SystemExit(
                "Code in ../compy_pressure/ or scripts/ changed during the run; "
                f"stopping before {stage}."
            )

    def selected_seeds(self):
        """SEED IDs of the selected stations in the fresh masters (smoke tests)."""
        seeds = set()
        for name in self.networks:
            master = self.workspace(name) / CALIBRATION_MASTER / "master_inventory.csv"
            if not master.is_file():
                continue
            with open(master, newline="") as stream:
                for row in csv.DictReader(stream):
                    if row["Seed_id"].split(".")[1] in self.stations:
                        seeds.add(row["Seed_id"])
        return sorted(seeds)

    # ------------------------------------------------------------ stages
    def preflight(self):
        problems = []
        if sys.version_info[:2] != (3, 13):
            problems.append(f"Python {sys.version.split()[0]}; supported is 3.13")
        for module in ("numpy", "scipy", "obspy", "pandas", "matplotlib", "disba"):
            try:
                __import__(module)
            except ImportError as error:
                problems.append(f"cannot import {module}: {error}")
        missing = [
            str(p)
            for p in required_inputs(self.config, self.networks)
            if not p.exists()
        ]
        problems += [f"missing input: {p}" for p in missing]
        state = code_state()
        if state["dirty"] and not self.args.allow_dirty:
            problems.append(
                "uncommitted code in ../compy_pressure/ or scripts/ "
                "(commit it or pass --allow-dirty):\n" + state["dirty"]
            )
        free = shutil.disk_usage(ROOT).free / 1e9
        print(
            f"Code {state['head'][:10]} state {state['state']}; free disk {free:.0f} GB"
        )
        print(
            f"Networks: {' '.join(self.networks)}; stations: {self.stations or 'all'}"
        )
        for line in problems:
            print("PROBLEM:", line)
        return problems, state, free

    def workspaces(self):
        from compy_pressure.inventory import load_master_dataset, save_master_dataset

        details = {}
        for name in self.networks:
            network = self.config["networks"][name]
            workspace = self.workspace(name)
            calibration = workspace / CALIBRATION_MASTER
            if (calibration / "master_inventory.csv").exists():
                details[name] = {"calibration_master": "already built"}
                continue
            if network["workspace"] == "inventory":
                master = inventory_master(
                    [resolve(p) for p in network["inventories"]],
                    workspace / "_inventory",
                )
            else:
                # The preparation selects its network's rows from this master and
                # rebuilds its own master from the seeded saved responses.
                source = workspace / "_source"
                if not (source / "_master/master_inventory.csv").exists():
                    save_master_dataset(
                        inventory_master(
                            [resolve(p) for p in network["source_inventories"]],
                            source / "_inventory",
                        ),
                        source / "_master",
                    )
                seed_prepare_inputs(network, workspace)
                status, _ = run_command(
                    [
                        self.python,
                        resolve(network["prepare_script"]),
                        "--results-dir",
                        workspace,
                        "--source-master",
                        source / "_master",
                        "--offline",
                        *network["prepare_args"],
                    ],
                    self.log(f"workspace_{name}"),
                    self.env,
                )
                if status:
                    raise RuntimeError(f"{name} preparation failed (exit {status})")
                master = load_master_dataset(workspace / "_master")
            master, counts = write_snapshots(master, workspace)
            save_master_dataset(master, calibration)
            details[name] = counts | {"channel_epochs": len(master)}
            print(f"{name}: {details[name]}")
        return details

    def rayleigh_runs(self, runs):
        results = {}
        for name, policy in runs:
            self.same_code(f"rayleigh {name} {policy}")
            command = rayleigh_command(
                self.config, name, policy, self.workspace(name), self.stations
            )
            status, _ = run_command(
                command, self.log(f"rayleigh_{name}_{policy}"), self.env
            )
            workflow = (
                self.workspace(name) / f"rayleigh_{policy}/_workflow/rayleigh.json"
            )
            record = {"exit": status, "workflow": str(workflow)}
            if workflow.is_file():
                record["status_counts"] = json.loads(workflow.read_text()).get(
                    "status_counts"
                )
            results[f"{name}:{policy}"] = record
        return results

    def identity(self):
        control = self.config["identity_control"]
        name, seed = control["network"], control["seed_id"]
        station = seed.split(".")[1]
        if name not in self.networks or (
            self.stations and station not in self.stations
        ):
            return {"skipped": "control station not selected"}
        result = self.rayleigh_runs([(name, control["policy"])])
        workflow = (
            self.workspace(name)
            / f"rayleigh_{control['policy']}/_workflow/rayleigh.json"
        )
        value = identity_value(workflow, seed) if workflow.is_file() else None
        result["identity"] = {
            "seed_id": seed,
            "expected": control["g"],
            "actual": value,
        }
        print(f"\nIdentity control {seed}: {value} (frozen {control['g']})")
        if value != control["g"]:
            self.run.mark("identity", "failed", **result)
            raise IdentityFailure(
                f"{seed} returned {value}, not the frozen {control['g']}; "
                "the fresh inputs do not reproduce the reference. Run stopped."
            )
        return result

    def coverage(self):
        campaign = self.run.root / "coverage/campaign"
        sources = coverage_config(
            resolve(self.config["coverage_template"]),
            self.networks,
            self.run.root / "coverage/coverage_sources.json",
        )
        base = [
            self.python,
            "-m",
            "compy_pressure.coverage",
        ]
        common = [
            "--campaign-dir",
            campaign,
            "--config",
            sources,
            "--source-root",
            ROOT,
        ]
        status, _ = run_command(
            [*base, "plan", *common, "--ppsd-scope", "all-channels"],
            self.log("coverage_plan"),
            self.env,
        )
        if status:
            raise RuntimeError(f"Coverage plan failed (exit {status})")
        command = [
            *base,
            "run",
            *common,
            "--ppsd-scope",
            "all-channels",
            "--execute",
            "--workers",
            self.args.workers,
        ]
        for seed in self.selected_seeds() if self.stations else []:
            command += ["--only-seed", seed]
        status, _ = run_command(command, self.log("coverage_execute"), self.env)
        counts = task_counts(campaign)
        if status:
            raise RuntimeError(f"Coverage execution reported failures (exit {status})")
        return {"campaign": str(campaign), "tasks": counts}

    def rayleigh(self):
        control = self.config["identity_control"]
        # The identity stage already ran the control network/policy, if selected.
        ran_control = "identity" in self.run.stage("identity")
        runs = [
            (name, policy)
            for name in self.networks
            for policy in self.config["event_policies"]
            if not (
                ran_control
                and (name, policy) == (control["network"], control["policy"])
            )
        ]
        results = self.rayleigh_runs(runs)
        failed = [k for k, v in results.items() if v["exit"]]
        if failed:
            raise RuntimeError(f"Rayleigh runs failed: {failed}")
        return results

    def microseism(self):
        sources = {
            name: {
                "master": str(self.workspace(name) / CALIBRATION_MASTER),
                "data": str(data_dir(self.config, name, self.workspace(name))),
            }
            for name in self.networks
        }
        folder = self.run.root / "microseism"
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "sources.json").write_text(json.dumps(sources, indent=2) + "\n")
        command = [
            self.python,
            ROOT / "scripts/run_microseism_check.py",
            "--sources",
            folder / "sources.json",
            "--output-root",
            folder,
        ]
        for name in self.networks:
            command += ["--network", name]
        for name, start in self.config.get("microseism_starts", {}).items():
            if name in self.networks:
                command += ["--start", f"{name}={start}"]
        for station in self.stations:
            command += ["--station", station]
        status, lines = run_command(command, self.log("microseism"), self.env)
        if status:
            raise RuntimeError(f"Microseism check failed (exit {status})")
        run_folder = lines[-1].rsplit(": ", 1)[-1].strip() if lines else ""
        record = {"run": run_folder}
        if run_folder and Path(run_folder).is_dir():
            compare_status, _ = run_command(
                [
                    self.python,
                    ROOT / "scripts/compare_microseism_methods.py",
                    run_folder,
                ],
                self.log("microseism_compare"),
                self.env,
            )
            record["comparison_exit"] = compare_status
            summary = Path(run_folder) / "summary.json"
            if summary.is_file():
                record["status_counts"] = json.loads(summary.read_text()).get(
                    "status_counts"
                )
        return record

    def local_p(self):
        campaign = (
            ROOT / "results/LOCAL_P/campaigns" / f"from_scratch_{self.run.root.name}"
        )
        command = [
            self.python,
            ROOT / "scripts/run_local_p_campaign.py",
            "--campaign-dir",
            campaign,
            "--networks",
            *self.networks,
        ]
        if self.stations:
            command += ["--stations", *self.stations]
        status, _ = run_command(command, self.log("local_p"), self.env)
        if status:
            raise RuntimeError(f"Local-P campaign failed (exit {status})")
        table = campaign / "station_table.csv"
        status, _ = run_command(
            [
                self.python,
                ROOT / "scripts/summarize_local_p_campaign.py",
                "--campaign-dir",
                campaign,
                "--csv",
                table,
                "--markdown",
                campaign / "station_table.md",
            ],
            self.log("local_p_table"),
            self.env,
        )
        record = {"campaign": str(campaign), "table_exit": status}
        export_status, _ = run_command(
            [
                self.python,
                ROOT / "scripts/export_local_p_results.py",
                "--campaign-dir",
                campaign,
                "--output-dir",
                f"{campaign}_export",
            ],
            self.log("local_p_export"),
            self.env,
        )
        record["export_exit"] = export_status
        if table.is_file():
            with open(table, newline="") as stream:
                record["status_counts"] = dict(
                    Counter(row["status"] for row in csv.DictReader(stream))
                )
        return record

    def registry(self):
        campaign = self.run.root / "coverage/campaign"
        masters = [
            self.workspace(name) / CALIBRATION_MASTER / "master_inventory.csv"
            for name in self.networks
        ]
        rows = registry_rows(campaign, [m for m in masters if m.is_file()])
        path = write_registry(rows, self.run.root / "registry/master_inventory.csv")
        counts = Counter(
            f"{r['Seed_id'].split('.')[0]}:{r['ppsd_units']}" for r in rows
        )
        return {
            "registry": str(path),
            "ppsds": len(rows),
            "by_network_units": dict(counts),
        }

    def summary(self):
        record = {"run": str(self.run.root), "stages": self.run.data.get("stages", {})}
        comparisons = {}
        for name in self.networks:
            network = self.config["networks"][name]
            policy = network.get("frozen_policy")
            workflow = (
                self.workspace(name) / f"rayleigh_{policy}/_workflow/rayleigh.json"
            )
            frozen = resolve(network["frozen_workflow"]) if policy else None
            if policy and workflow.is_file() and frozen and frozen.is_file():
                comparisons[name] = compare_workflows(workflow, frozen) | {
                    "policy": policy
                }
        record["rayleigh_vs_frozen"] = comparisons
        (self.run.root / "summary.json").write_text(
            json.dumps(record, indent=2, default=str) + "\n"
        )
        lines = [f"# From-scratch run {self.run.root.name}", ""]
        lines.append("| Stage | Status |")
        lines.append("|---|---|")
        for stage in STAGES:
            lines.append(
                f"| {stage} | {self.run.stage(stage).get('status', 'not run')} |"
            )
        lines += ["", "## Rayleigh against the frozen workflows", ""]
        lines.append("| Network | Policy | Stations | Identical | Median gain change |")
        lines.append("|---|---|---:|---:|---:|")
        for name, c in comparisons.items():
            median = (
                f"{100 * c['median_gain_change']:+.2f}%"
                if c["median_gain_change"] is not None
                else "-"
            )
            lines.append(
                f"| {name} | {c['policy']} | {c['stations_compared']} | "
                f"{c['identical']} | {median} |"
            )
        (self.run.root / "SUMMARY.md").write_text("\n".join(lines) + "\n")
        return {"rayleigh_vs_frozen": comparisons}

    # ------------------------------------------------------------ driver
    def selected_stages(self):
        stages = [s for s in STAGES if not self.args.only or s in self.args.only]
        return [s for s in stages if s not in (self.args.skip or [])]

    def execute(self):
        problems, state, free = self.preflight()
        if problems:
            return 2
        if self.args.check:
            print("\nPlanned stages:", " ".join(self.selected_stages()))
            print("Preflight passed. Nothing was run or written.")
            return 0
        if self.args.run_root:
            self.run = Run(self.args.run_root)
            if not self.run.path.exists():
                raise SystemExit(f"No run.json in {self.args.run_root}")
            if (
                self.run.data["networks"] != self.networks
                or self.run.data["stations"] != self.stations
            ):
                raise SystemExit("Resume must keep the original networks and stations")
        else:
            stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
            self.run = Run(
                ROOT / "results/FROM_SCRATCH" / f"run_{stamp}_{uuid4().hex[:6]}"
            )
            self.run.data.update(
                created_utc=utc_now(),
                networks=self.networks,
                stations=self.stations,
                config=self.config["_path"],
                config_sha256=self.config["_sha256"],
                code=state,
                python=sys.version,
            )
            self.run.save()
        if self.run.data["code"]["state"] != state["state"]:
            raise SystemExit("Code differs from the run's recorded code; use a new run")
        print(f"Run folder: {self.run.root}")
        failed = []
        for stage in self.selected_stages():
            if self.run.done(stage):
                print(f"\n== {stage}: already complete")
                continue
            print(f"\n== {stage} ({utc_now()})")
            self.same_code(stage)
            self.run.mark(stage, "running")
            try:
                details = (
                    {"free_disk_gb": round(free)}
                    if stage == "preflight"
                    else getattr(self, stage)()
                )
            except IdentityFailure as error:
                print(f"\nSTOPPED: {error}")
                return 3
            except Exception as error:  # record and continue with independent stages
                self.run.mark(stage, "failed", error=f"{type(error).__name__}: {error}")
                print(f"\n{stage} FAILED: {error}")
                failed.append(stage)
                if stage == "workspaces":
                    print("Calibration stages need the workspaces; stopping.")
                    return 1
                continue
            self.run.mark(stage, "complete", **(details or {}))
        print(f"\nRun folder: {self.run.root}")
        print("Failed stages:", failed or "none")
        return 1 if failed else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--check", action="store_true", help="preflight only")
    parser.add_argument("--run-root", type=Path, help="resume this run folder")
    parser.add_argument("--networks", nargs="+", help="default: all configured")
    parser.add_argument("--stations", nargs="+", help="subset, e.g. a smoke test")
    parser.add_argument("--only", action="append", choices=STAGES)
    parser.add_argument("--skip", action="append", choices=STAGES)
    parser.add_argument("--workers", type=int, default=4, help="coverage workers")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()
    return Pipeline(args, load_config(args.config)).execute()


if __name__ == "__main__":
    sys.exit(main())
