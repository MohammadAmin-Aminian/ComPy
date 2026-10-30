"""Validate an epoch-resolved local-P export; update only with explicit --apply.

Default is a read-only plan. This never assigns calibration scores or changes
the six scientific check-status columns. Ambiguous/not-run export rows are
reported and skipped, not guessed onto a master epoch.
"""

import argparse
import csv
import json
from pathlib import Path

from compy_pressure.inventory import (
    EPOCH_COLUMNS,
    load_master_dataset,
    update_local_p_results,
)

KEYS = EPOCH_COLUMNS + ["inventory_file"]


def validated_updates(export_path, master_directory):
    master = load_master_dataset(master_directory)
    counts = master.groupby(KEYS, dropna=False).size().to_dict()
    updates, skipped, seen = [], [], set()
    with Path(export_path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        if row.get("master_epoch_match") != "exact":
            skipped.append(row.get("Seed_id", ""))
            continue
        key = tuple(row.get(k, "") for k in KEYS)
        if key in seen or counts.get(key) != 1:
            raise ValueError(f"Duplicate, ambiguous or missing master epoch: {key}")
        seen.add(key)
        updates.append(
            {k: row[k] for k in KEYS}
            | {
                k: v
                for k, v in row.items()
                if k.startswith("local_p_") and not k.endswith("_score")
            }
        )
    return updates, skipped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export", type=Path, required=True)
    parser.add_argument("--master-dir", type=Path, required=True)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Explicitly write local_p_* fields with master backups",
    )
    args = parser.parse_args()
    updates, skipped = validated_updates(args.export, args.master_dir)
    if args.apply:
        update_local_p_results(args.master_dir, updates)
    print(
        json.dumps(
            {
                "applied": args.apply,
                "unique_exact_updates": len(updates),
                "skipped_unresolved_or_not_run": len(skipped),
                "scores_and_check_statuses": "unchanged",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
