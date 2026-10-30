#!/usr/bin/env bash
# Run as `bash scripts/run_three_methods.sh`, never source into the terminal.
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  printf 'Run this file with bash; do not source it.\n' >&2
  return 2
fi
set -Eeuo pipefail

BRUIT_MAIN="$(cd "$(dirname "$0")/.." && pwd)"
BRUIT_METHODS="${BRUIT_MICROSEISM_CHECKOUT:-${BRUIT_CLAUDE_CHECKOUT:-$BRUIT_MAIN}}"
BRUIT_STAGE=preflight
BRUIT_SESSION=
trap 'BRUIT_STATUS=$?; printf "\nStopped during %s (exit %s, line %s).\nYour terminal remains open. Results/logs: %s\n" "$BRUIT_STAGE" "$BRUIT_STATUS" "$LINENO" "${BRUIT_SESSION:-not started}" >&2; exit "$BRUIT_STATUS"' ERR

BRUIT_CHECK_ONLY=0
BRUIT_RERUN_RAYLEIGH=0
BRUIT_RAYLEIGH=
while [[ $# -gt 0 ]]; do
  case "$1" in
    --check) BRUIT_CHECK_ONLY=1; shift ;;
    --rerun-rayleigh) BRUIT_RERUN_RAYLEIGH=1; shift ;;
    --rayleigh-run)
      [[ $# -ge 2 ]] || { printf '%s needs a directory\n' "$1" >&2; exit 2; }
      BRUIT_RAYLEIGH="$2"; shift 2 ;;
    *) printf 'Usage: bash %s [--check] [--rayleigh-run DIR | --rerun-rayleigh]\n' "$0" >&2; exit 2 ;;
  esac
done
if [[ -n "$BRUIT_RAYLEIGH" && "$BRUIT_RERUN_RAYLEIGH" -eq 1 ]]; then
  printf 'Choose --rayleigh-run or --rerun-rayleigh.\n' >&2
  exit 2
fi

for BRUIT_FILE in \
  "$BRUIT_MAIN/scripts/run_calibration_experiment.py" \
  "$BRUIT_METHODS/scripts/run_microseism_check_all.sh" \
  "$BRUIT_METHODS/scripts/compare_microseism_methods.py"; do
  if [[ ! -f "$BRUIT_FILE" ]]; then
    printf 'Required file missing: %s\n' "$BRUIT_FILE" >&2
    exit 2
  fi
done

if command -v conda >/dev/null; then
  BRUIT_PY=(conda run --no-capture-output -n bruitfm python)
else
  BRUIT_PY=(/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python)
fi
PYTHONPATH="$BRUIT_METHODS/..:$BRUIT_METHODS" "${BRUIT_PY[@]}" -c \
  'import compy_pressure.microseism_check, compy_pressure.microseism_diagnostic, scripts.run_microseism_check; import disba; print("Method imports: OK")'

if [[ -z "$BRUIT_RAYLEIGH" && "$BRUIT_RERUN_RAYLEIGH" -eq 0 ]]; then
  BRUIT_RAYLEIGH="$("${BRUIT_PY[@]}" - "$BRUIT_MAIN" <<'PY'
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
if [[ "$BRUIT_RERUN_RAYLEIGH" -eq 0 ]]; then
  test -f "$BRUIT_RAYLEIGH/experiment.json"
  test -f "$BRUIT_RAYLEIGH/comparison.csv"
  printf 'Reusing saved Rayleigh results: %s\n' "$BRUIT_RAYLEIGH"
fi

if [[ "$BRUIT_CHECK_ONLY" -eq 1 ]]; then
  printf 'Preflight passed. No analysis was run.\n'
  exit 0
fi

BRUIT_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
BRUIT_SESSION="$BRUIT_MAIN/results/EXPERIMENTS/THREE_METHODS/run_${BRUIT_STAMP}_$$"
mkdir -p "$(dirname "$BRUIT_SESSION")"
mkdir "$BRUIT_SESSION"
printf 'Session: %s\nMain commit: %s\nMicroseism code commit: %s\n' \
  "$BRUIT_SESSION" \
  "$(git -C "$BRUIT_MAIN" rev-parse --short HEAD)" \
  "$(git -C "$BRUIT_METHODS" rev-parse --short HEAD)" \
  | tee "$BRUIT_SESSION/session.log"

BRUIT_STAGE=Rayleigh
if [[ "$BRUIT_RERUN_RAYLEIGH" -eq 1 ]]; then
  BRUIT_RAYLEIGH="$BRUIT_SESSION/rayleigh"
  PYTHONPATH="$BRUIT_MAIN/..:$BRUIT_MAIN" "${BRUIT_PY[@]}" \
    "$BRUIT_MAIN/scripts/run_calibration_experiment.py" \
    --network all --all-stations --output-dir "$BRUIT_RAYLEIGH" \
    2>&1 | tee "$BRUIT_SESSION/rayleigh.log"
else
  printf 'Reused: %s\n' "$BRUIT_RAYLEIGH" | tee "$BRUIT_SESSION/rayleigh.log"
fi

BRUIT_STAGE=Rayleigh-comparison-table
"${BRUIT_PY[@]}" - "$BRUIT_RAYLEIGH" "$BRUIT_SESSION/gain_table.csv" <<'PY' \
  2>&1 | tee "$BRUIT_SESSION/gain_table.log"
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

BRUIT_STAGE=microseism-methods
bash "$BRUIT_METHODS/scripts/run_microseism_check_all.sh" \
  --gain-table "$BRUIT_SESSION/gain_table.csv" \
  --output-root "$BRUIT_SESSION/microseism" \
  2>&1 | tee "$BRUIT_SESSION/microseism.log"

BRUIT_MICRO="$(sed -n 's/^Run folder: //p' "$BRUIT_SESSION/microseism.log" | tail -n 1)"
test -n "$BRUIT_MICRO"
test -f "$BRUIT_MICRO/method_comparison.csv"
printf '\nFinished.\nRayleigh: %s\nComparison: %s\nLogs: %s\n' \
  "$BRUIT_RAYLEIGH" "$BRUIT_MICRO/method_comparison.csv" "$BRUIT_SESSION"
