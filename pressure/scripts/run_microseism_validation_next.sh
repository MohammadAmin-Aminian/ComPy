#!/usr/bin/env bash
# Continue only the missing validation weeks; completed shifts are reused.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BASE_ROOT="${PRESSURE_MICROSEISM_BASE_ROOT:-$(dirname "$ROOT")/pressure-qc-claude/results/EXPERIMENTS/MICROSEISM_ALL}"

YV="$BASE_ROOT/run_20261014T183508Z_1f096dca"
YL="$BASE_ROOT/run_20261014T183501Z_0a1ff170"
SEVEN_D="$BASE_ROOT/run_20261014T183456Z_61c4966a"
XO="$BASE_ROOT/run_20261014T183418Z_736f654a"
for BASE in "$YV" "$YL" "$SEVEN_D" "$XO"; do
  test -f "$BASE/settings.json"
done

# 7D has only 33–41 days on disk. YV has the longer deployment coverage.
bash "$ROOT/scripts/run_microseism_weeks.sh" "$SEVEN_D" -2 -1 1 2
bash "$ROOT/scripts/run_microseism_weeks.sh" "$YV" -12 -6 6 12

if command -v conda >/dev/null; then
  PRESSURE_PY=(conda run --no-capture-output -n compy python)
else
  PRESSURE_PY=(/opt/homebrew/Caskroom/miniforge/base/envs/compy/bin/python)
fi
PYTHONPATH="$ROOT/..:$ROOT" "${PRESSURE_PY[@]}" \
  "$ROOT/scripts/combine_microseism_weeks.py" "$YV" "$YL" "$SEVEN_D" "$XO" \
  --docs "$ROOT/docs"
