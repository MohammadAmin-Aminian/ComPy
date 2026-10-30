#!/usr/bin/env bash
# Run all 46 frozen-baseline microseism diagnostics sequentially.
set -uo pipefail

cd "$(dirname "$0")/.."

BRUIT_PY=(conda run --no-capture-output -n bruitfm python)
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_ROOT="results/DIAGNOSTICS/MICROSEISM_BASELINE_LOGS/run_${RUN_ID}"
SUMMARY="$LOG_ROOT/summary.tsv"
mkdir -p "$LOG_ROOT"
printf 'station\tstatus\tlog\n' >"$SUMMARY"

STATIONS=()
while IFS= read -r station; do
  STATIONS+=("$station")
done < <(
  PYTHONPATH=..:. "${BRUIT_PY[@]}" -m compy_pressure.microseism_campaign --list |
    awk '{print $1}'
)

if [[ ${#STATIONS[@]} -ne 46 ]]; then
  printf 'Expected 46 baseline stations, found %s\n' "${#STATIONS[@]}" >&2
  exit 2
fi

failures=0
for station in "${STATIONS[@]}"; do
  log="$LOG_ROOT/${station}.log"
  printf '\n===== %s =====\n' "$station"
  if PYTHONPATH=..:. "${BRUIT_PY[@]}" -m compy_pressure.microseism_campaign \
    --station "$station" 2>&1 | tee "$log"; then
    status=complete
  else
    status=failed
    failures=$((failures + 1))
  fi
  printf '%s\t%s\t%s\n' "$station" "$status" "$log" >>"$SUMMARY"
done

printf '\nSummary: %s\n' "$SUMMARY"
printf 'Completed: %s; failed: %s\n' "$((46 - failures))" "$failures"
[[ $failures -eq 0 ]]
