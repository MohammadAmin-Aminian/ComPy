"""Read-only check that a copied BRUIT-FM workspace still supports reuse.

Run it on the old computer first (baseline), then on the new computer after
copying, before planning or computing anything. It checks:

* Python version and the scientific imports;
* the data, inventory, master and results folders;
* every registered saved PPSD: NPZ, ``processing.json`` and each recorded
  source waveform at the same path with identical size and mtime_ns. The
  coverage planner needs all three to reuse a product: a changed mtime blocks
  reuse, and a changed path makes the files look uncovered;
* frozen reference files with their recorded SHA-256, and every absolute path
  the frozen six-check review points to (several live in sibling folders).

Nothing is written unless ``--json-out`` names a new file. Exit status 0 means
every check passed, 1 means something must be fixed before planning.
"""

import argparse
import csv
import hashlib
import json
import os
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = Path.home() / "Desktop/Data_Bruit_FM"
FROZEN = {
    "results/NETWORKS/YV_RAYLEIGH_NO_TILT/_workflow/rayleigh.json": (
        "cbe5af4ae7782903d96ea5f6eeff3c07ca3c57e85b44f4f3d55c68f9f85486d2"
    ),
    "results/EXPERIMENTS/SOCRATES/all_20261013T195356Z_00cf3343/experiment.json": (
        "446b8e36b94c9952736b9a2c9b2353f0dfb0c95e05e8cc3719b9cea72fcc20cf"
    ),
    "results/REVIEWS/SIX_CHECKS/final_20261019/review.json": (
        "241e9ab125a718fe92958b2129a267ebeff43076287142d85c0beb6d8cea1fd2"
    ),
    str(DATA / "INVENTORIES/YV.xml"): (
        "58afc2d6f0056e52cf954bdd5b9fb80e92f9572d44d067b30966a160d2044494"
    ),
}
RR52_G = "0.6787278140712512"


def sha256(path):
    with open(path, "rb") as stream:
        # hashlib.file_digest needs Python 3.11; ComPy supports 3.10
        checksum = hashlib.sha256()
        for block in iter(lambda: stream.read(1 << 20), b""):
            checksum.update(block)
        return checksum.hexdigest()


def check_environment():
    problems = []
    if sys.version_info[:2] != (3, 13):
        problems.append(
            f"Python {sys.version.split()[0]}; the supported version is 3.13"
        )
    for module in ("numpy", "scipy", "obspy", "matplotlib", "pandas", "disba"):
        try:
            __import__(module)
        except ImportError as error:
            problems.append(f"cannot import {module}: {error}")
    return problems


def check_folders(root, data):
    required = [
        data / "Data",
        data / "INVENTORIES",
        root / "exports/master_inventory.csv",
        root / "results",
    ]
    return [f"missing: {p}" for p in required if not p.exists()]


def source_state(item):
    real = os.path.realpath(item["path"])
    if not Path(real).is_file():
        return "missing"
    stat = Path(real).stat()
    same = (stat.st_size, stat.st_mtime_ns) == (item["size"], item["mtime_ns"])
    return "unchanged" if same else "modified"


def registries(root):
    """Product registries the coverage planner reads (docs/coverage_sources.json)."""
    config = json.loads((root / "docs/coverage_sources.json").read_text())
    return [root / r for r in config.get("product_registries", [])]


def check_ppsds(registry_paths, root):
    """Saved PPSD products and their recorded waveform sources, per network.

    Products without processing.json cannot be reused even before a move
    (the planner treats that as a blocker), so they are counted separately.
    """
    products, sources, examples = Counter(), Counter(), []
    seen = set()
    for registry in registry_paths:
        if not registry.is_file():
            products[("-", "missing_registry")] += 1
            examples.append(f"missing registry: {registry}")
            continue
        with open(registry, newline="") as stream:
            for row in csv.DictReader(stream):
                npz = row.get("ppsd_npz_file", "")
                if row.get("ppsd_status") not in {"complete", "partial"} or not npz:
                    continue
                npz = Path(npz).expanduser()
                if not npz.is_absolute():
                    npz = root / npz
                if npz in seen:
                    continue
                seen.add(npz)
                network = row["Seed_id"].split(".")[0]
                manifest = npz.with_name("processing.json")
                if not npz.is_file():
                    products[(network, "missing_npz")] += 1
                    examples.append(f"missing PPSD: {npz}")
                    continue
                if not manifest.is_file():
                    products[(network, "no_processing_json")] += 1
                    continue
                products[(network, "present")] += 1
                for item in json.loads(manifest.read_text()).get("source_files", []):
                    state = source_state(item)
                    sources[(network, state)] += 1
                    if state != "unchanged" and len(examples) < 10:
                        examples.append(f"{state} source: {item['path']}")
    return products, sources, examples


def check_frozen(root):
    rows = []
    for name, expected in FROZEN.items():
        path = Path(name) if Path(name).is_absolute() else root / name
        if not path.is_file():
            rows.append((name, "missing"))
        else:
            rows.append((name, "ok" if sha256(path) == expected else "hash differs"))
    workflow = root / "results/NETWORKS/YV_RAYLEIGH_NO_TILT/_workflow/rayleigh.json"
    if workflow.is_file():
        entries = json.loads(workflow.read_text())["entries"].values()
        g = next(
            (
                e["update"].get("rayleigh_pressure_multiplier")
                for e in entries
                if e["seed_id"] == "YV.RR52.00.BDH"
            ),
            None,
        )
        rows.append(("RR52 frozen g", "ok" if str(g) == RR52_G else f"differs: {g}"))
    return rows


def check_review_pointers(root):
    review = root / "results/REVIEWS/SIX_CHECKS/final_20261019/review.json"
    if not review.is_file():
        return None, []
    paths = set(re.findall(r'"(/[^"]+)"', review.read_text()))
    paths = {p for p in paths if p.startswith(("/Users/", "/home/"))}
    missing = sorted(p for p in paths if not Path(p).exists())
    folders = Counter("/".join(p.split("/")[:5]) for p in missing)
    return len(paths), [f"{n} missing under {folder}" for folder, n in folders.items()]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--root", default=str(ROOT), help="project checkout")
    parser.add_argument("--data", default=str(DATA), help="Data_Bruit_FM folder")
    parser.add_argument("--json-out", help="optional new JSON report file")
    args = parser.parse_args()
    root, data = Path(args.root).resolve(), Path(args.data).expanduser()
    failures = []

    print(f"Project: {root}\nData:    {data}\nUser:    {Path.home()}")
    env = check_environment()
    print("\n[1] Environment:", "OK" if not env else "PROBLEMS")
    for line in env:
        print("   ", line)
    failures += env

    folders = check_folders(root, data)
    print("[2] Folders:", "OK" if not folders else "PROBLEMS")
    for line in folders:
        print("   ", line)
    failures += folders

    master = root / "exports/master_inventory.csv"
    if master.is_file():
        products, sources, examples = check_ppsds(registries(root), root)
        bad_products = sum(
            n
            for (_, s), n in products.items()
            if s in {"missing_npz", "missing_registry"}
        )
        bad_sources = sum(n for (_, s), n in sources.items() if s != "unchanged")
        print(
            "[3] Saved PPSDs reusable:",
            "OK" if not (bad_products or bad_sources) else "PROBLEMS",
        )
        for network in sorted({n for n, _ in products} | {n for n, _ in sources}):
            p = {s: n for (k, s), n in products.items() if k == network}
            s = {st: n for (k, st), n in sources.items() if k == network}
            print(f"    {network}: products {p}; source files {s}")
        for line in examples:
            print("   ", line)
        if bad_products or bad_sources:
            failures.append(
                f"{bad_products} PPSD products and {bad_sources} source files are "
                "missing or changed; the planner would recompute them"
            )
    else:
        products, sources = Counter(), Counter()

    frozen = check_frozen(root)
    print(
        "[4] Frozen references:",
        "OK" if all(s == "ok" for _, s in frozen) else "PROBLEMS",
    )
    for name, state in frozen:
        print(f"    {state:13s} {name}")
    failures += [f"{n}: {s}" for n, s in frozen if s != "ok"]

    total, missing = check_review_pointers(root)
    if total is not None:
        print(
            f"[5] Paths in the frozen review ({total}):",
            "OK" if not missing else "PROBLEMS",
        )
        for line in missing:
            print("   ", line)
        failures += missing

    print(
        "\nRESULT:",
        "READY" if not failures else f"NOT READY ({len(failures)} problems)",
    )
    if args.json_out:
        with open(args.json_out, "x") as handle:
            json.dump(
                {
                    "root": str(root),
                    "data": str(data),
                    "products": {f"{k[0]}:{k[1]}": v for k, v in products.items()},
                    "sources": {f"{k[0]}:{k[1]}": v for k, v in sources.items()},
                    "failures": failures,
                },
                handle,
                indent=2,
            )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
