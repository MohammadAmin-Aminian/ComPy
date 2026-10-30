#!/usr/bin/env bash
# Repeat the microseism check on further weeks: each station's week in a
# completed base run, moved by whole weeks (no fallback; stations outside
# coverage stay unresolved), then the method comparison for each new run.
# Diagnostic only; one new run folder per shift. Default shifts: -12 -6 6 12.
# Usage: bash scripts/run_microseism_weeks.sh BASE_RUN [SHIFT ...]
set -eo pipefail
if [ $# -lt 1 ]; then
  echo "usage: bash scripts/run_microseism_weeks.sh BASE_RUN [SHIFT ...]" >&2
  exit 2
fi
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BASE="$1"
shift
SHIFTS=("$@")
[ ${#SHIFTS[@]} -gt 0 ] || SHIFTS=(-12 -6 6 12)
if command -v conda >/dev/null; then
  BRUIT_PY=(conda run --no-capture-output -n bruitfm python)
else
  BRUIT_PY=(/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python)
fi
export PYTHONPATH="$ROOT/..:$ROOT"  # this branch's code, not the main install
LOG="$(mktemp -t microseism_weeks)"
trap 'rm -f "$LOG"' EXIT
for K in "${SHIFTS[@]}"; do
  if [[ ! "$K" =~ ^-?[0-9]+$ || "$K" -eq 0 ]]; then
    printf 'Week shifts must be nonzero integers: %s\n' "$K" >&2
    exit 2
  fi
  # Resume by exact base/shift binding; completed waves/spectra are not rerun.
  EXISTING="$("${BRUIT_PY[@]}" - "$BASE" "$K" <<'PY'
import sys
from scripts.combine_microseism_weeks import week_runs

for shift, run in week_runs(sys.argv[1]):
    if shift == int(sys.argv[2]):
        print(run)
PY
)"
  if [[ -n "$EXISTING" ]]; then
    printf 'Week shift %s: reusing completed %s\n' "$K" "$EXISTING"
    continue
  fi
  "${BRUIT_PY[@]}" "$ROOT/scripts/run_microseism_check.py" \
    --base-run "$BASE" --shift-weeks "$K" | tee "$LOG"
  RUN="$(tail -1 "$LOG" | sed 's/.*: //')"
  "${BRUIT_PY[@]}" "$ROOT/scripts/compare_microseism_methods.py" "$RUN"
  echo "Week shift $K: $RUN"
done
