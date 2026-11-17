#!/usr/bin/env bash
# Run as `bash scripts/run_three_methods.sh`, never source into the terminal.
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  printf 'Run this file with bash; do not source it.\n' >&2
  return 2
fi
set -Eeuo pipefail

PRESSURE_MAIN="$(cd "$(dirname "$0")/.." && pwd)"
PRESSURE_METHODS="${PRESSURE_MICROSEISM_CHECKOUT:-${PRESSURE_CLAUDE_CHECKOUT:-$PRESSURE_MAIN}}"
PRESSURE_STAGE=preflight
PRESSURE_SESSION=
trap 'PRESSURE_STATUS=$?; printf "\nStopped during %s (exit %s, line %s).\nYour terminal remains open. Results/logs: %s\n" "$PRESSURE_STAGE" "$PRESSURE_STATUS" "$LINENO" "${PRESSURE_SESSION:-not started}" >&2; exit "$PRESSURE_STATUS"' ERR

PRESSURE_CHECK_ONLY=0
PRESSURE_RERUN_RAYLEIGH=0
PRESSURE_RAYLEIGH=
while [[ $# -gt 0 ]]; do
  case "$1" in
    --check) PRESSURE_CHECK_ONLY=1; shift ;;
    --rerun-rayleigh) PRESSURE_RERUN_RAYLEIGH=1; shift ;;
    --rayleigh-run)
      [[ $# -ge 2 ]] || { printf '%s needs a directory\n' "$1" >&2; exit 2; }
      PRESSURE_RAYLEIGH="$2"; shift 2 ;;
    *) printf 'Usage: bash %s [--check] [--rayleigh-run DIR | --rerun-rayleigh]\n' "$0" >&2; exit 2 ;;
  esac
done
if [[ -n "$PRESSURE_RAYLEIGH" && "$PRESSURE_RERUN_RAYLEIGH" -eq 1 ]]; then
  printf 'Choose --rayleigh-run or --rerun-rayleigh.\n' >&2
  exit 2
fi

for PRESSURE_FILE in \
  "$PRESSURE_MAIN/scripts/run_calibration_experiment.py" \
  "$PRESSURE_METHODS/scripts/run_microseism_check_all.sh" \
  "$PRESSURE_METHODS/scripts/compare_microseism_methods.py"; do
  if [[ ! -f "$PRESSURE_FILE" ]]; then
    printf 'Required file missing: %s\n' "$PRESSURE_FILE" >&2
    exit 2
  fi
done

if command -v conda >/dev/null; then
  PRESSURE_PY=(conda run --no-capture-output -n compy python)
else
  PRESSURE_PY=(/opt/homebrew/Caskroom/miniforge/base/envs/compy/bin/python)
fi
PYTHONPATH="$PRESSURE_METHODS/..:$PRESSURE_METHODS" "${PRESSURE_PY[@]}" -c \
  'import compy_pressure.microseism_check, compy_pressure.microseism_diagnostic, scripts.run_microseism_check; import disba; print("Method imports: OK")'

if [[ -z "$PRESSURE_RAYLEIGH" && "$PRESSURE_RERUN_RAYLEIGH" -eq 0 ]]; then
  PRESSURE_RAYLEIGH="$("${PRESSURE_PY[@]}" - "$PRESSURE_MAIN" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

metadata = json.loads((Path(sys.argv[1]) / "docs/socrates_all_station_comparison.json").read_text())
source = Path(metadata["source_experiment"])
if hashlib.sha256(source.read_bytes()).hexdigest() != metadata["source_experiment_sha256"]:
    raise ValueError("Saved Rayleigh experiment hash changed")
print(source.parent)
PY
)"
fi
if [[ "$PRESSURE_RERUN_RAYLEIGH" -eq 0 ]]; then
  test -f "$PRESSURE_RAYLEIGH/experiment.json"
  test -f "$PRESSURE_RAYLEIGH/comparison.csv"
  printf 'Reusing saved Rayleigh results: %s\n' "$PRESSURE_RAYLEIGH"
fi

if [[ "$PRESSURE_CHECK_ONLY" -eq 1 ]]; then
  printf 'Preflight passed. No analysis was run.\n'
  exit 0
fi

PRESSURE_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
PRESSURE_SESSION="$PRESSURE_MAIN/results/EXPERIMENTS/THREE_METHODS/run_${PRESSURE_STAMP}_$$"
mkdir -p "$(dirname "$PRESSURE_SESSION")"
mkdir "$PRESSURE_SESSION"
printf 'Session: %s\nMain commit: %s\nMicroseism code commit: %s\n' \
  "$PRESSURE_SESSION" \
  "$(git -C "$PRESSURE_MAIN" rev-parse --short HEAD)" \
  "$(git -C "$PRESSURE_METHODS" rev-parse --short HEAD)" \
  | tee "$PRESSURE_SESSION/session.log"

PRESSURE_STAGE=Rayleigh
if [[ "$PRESSURE_RERUN_RAYLEIGH" -eq 1 ]]; then
  PRESSURE_RAYLEIGH="$PRESSURE_SESSION/rayleigh"
  PYTHONPATH="$PRESSURE_MAIN/..:$PRESSURE_MAIN" "${PRESSURE_PY[@]}" \
    "$PRESSURE_MAIN/scripts/run_calibration_experiment.py" \
    --network all --all-stations --output-dir "$PRESSURE_RAYLEIGH" \
    2>&1 | tee "$PRESSURE_SESSION/rayleigh.log"
else
  printf 'Reused: %s\n' "$PRESSURE_RAYLEIGH" | tee "$PRESSURE_SESSION/rayleigh.log"
fi

PRESSURE_STAGE=Rayleigh-comparison-table
"${PRESSURE_PY[@]}" - "$PRESSURE_RAYLEIGH" "$PRESSURE_SESSION/gain_table.csv" <<'PY' \
  2>&1 | tee "$PRESSURE_SESSION/gain_table.log"
import csv
import sys
from pathlib import Path

run = Path(sys.argv[1])
by_station = {}
with (run / "comparison.csv").open(newline="") as stream:
    for row in csv.DictReader(stream):
        key = row["network"], row["station"]
        variants = by_station.setdefault(key, {})
        if row["variant"] in variants:
            raise ValueError(f"Duplicate Rayleigh variant: {key} {row['variant']}")
        variants[row["variant"]] = row

output = []
for (network, station), variants in sorted(by_station.items()):
    def gain(*names):
        for name in names:
            value = variants.get(name, {}).get("pressure_multiplier_g", "")
            if value:
                return value
        return ""

    saved = variants["saved_before"]
    output.append({
        "network": network,
        "station": station,
        "previous_g": gain("saved_before"),
        "pole_trial_g": gain("pole_absolute_normal_state", "pole_absolute"),
        "weighted_trial_g": gain("pole_weighted_log_normal_state", "pole_weighted_log"),
        "comparison_direction": saved["coefficient_direction"],
        "published_coefficient": saved["published_coefficient"],
    })

if len(output) != 89:
    raise ValueError(f"Expected 89 station epochs, found {len(output)}")
destination = Path(sys.argv[2])
with destination.open("x", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(output[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(output)
print(f"Rayleigh comparison table: {destination}")
PY

PRESSURE_STAGE=microseism-methods
bash "$PRESSURE_METHODS/scripts/run_microseism_check_all.sh" \
  --gain-table "$PRESSURE_SESSION/gain_table.csv" \
  --output-root "$PRESSURE_SESSION/microseism" \
  2>&1 | tee "$PRESSURE_SESSION/microseism.log"

PRESSURE_MICRO="$(sed -n 's/^Run folder: //p' "$PRESSURE_SESSION/microseism.log" | tail -n 1)"
test -n "$PRESSURE_MICRO"
test -f "$PRESSURE_MICRO/method_comparison.csv"
printf '\nFinished.\nRayleigh: %s\nComparison: %s\nLogs: %s\n' \
  "$PRESSURE_RAYLEIGH" "$PRESSURE_MICRO/method_comparison.csv" "$PRESSURE_SESSION"
