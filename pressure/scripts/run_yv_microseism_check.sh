#!/usr/bin/env bash
# YV-only microseism check on the common 2013-06-10 week (the completed
# run_20261014T175553Z_30002f23 used this setting). Diagnostic only.
# Usage: bash scripts/run_yv_microseism_check.sh [--open-plots] [--resume RUN_DIR]
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
if command -v conda >/dev/null; then
  BRUIT_PY=(conda run --no-capture-output -n bruitfm python)
else
  BRUIT_PY=(/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python)
fi
PYTHONPATH="$ROOT/..:$ROOT" "${BRUIT_PY[@]}" "$ROOT/scripts/run_microseism_check.py" \
  --network YV --start YV=2013-06-10 \
  --output-root "$ROOT/results/EXPERIMENTS/MICROSEISM_YV" "$@"
