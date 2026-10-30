#!/usr/bin/env bash
# Microseism check at every YV/YL/7D/XO station, then the four-column method
# comparison (Claude microseism | Rayleigh | Codex diagnostic | published).
# Diagnostic only: reads masters/waveforms in place, writes one new run folder
# under this checkout's results/EXPERIMENTS/MICROSEISM_ALL/. ~25-35 min.
# Usage: bash scripts/run_microseism_check_all.sh [--network YL ...] [--resume RUN_DIR]
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
if command -v conda >/dev/null; then
  BRUIT_PY=(conda run --no-capture-output -n bruitfm python)
else
  BRUIT_PY=(/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python)
fi
export PYTHONPATH="$ROOT/..:$ROOT"  # this branch's code, not the main install
LOG="$(mktemp -t microseism_all)"
trap 'rm -f "$LOG"' EXIT
"${BRUIT_PY[@]}" "$ROOT/scripts/run_microseism_check.py" "$@" | tee "$LOG"
RUN="$(tail -1 "$LOG" | sed 's/.*: //')"
"${BRUIT_PY[@]}" "$ROOT/scripts/compare_microseism_methods.py" "$RUN"
echo "Run folder: $RUN"
echo "Table: $RUN/method_comparison.md  Figure: $RUN/method_comparison.png"
