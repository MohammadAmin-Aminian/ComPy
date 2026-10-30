#!/usr/bin/env bash
# Independent RR52 microseism pressure–vertical check (diagnostic only).
# Reads the saved paired spectra and the fixed Rayleigh reference; writes a new
# run folder. No gain, master, score or saved input is changed.
# Usage: bash scripts/run_rr52_microseism_check.sh [--show]
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# Saved inputs live in the main checkout; override with BRUITFM_RESULTS if moved.
SAVED="${BRUITFM_RESULTS:-$HOME/Projects/bruit-fm/results}"
OUTPUT="${MICROSEISM_OUTPUT_ROOT:-$ROOT/results/EXPERIMENTS/MICROSEISM}"
R="$SAVED/NETWORKS/TEST/YV/RR52/YV.RR52.00.BDH__20121019T090000__58afc2d6f005"
if command -v conda >/dev/null; then
  BRUIT_PY=(conda run --no-capture-output -n bruitfm python)
else
  BRUIT_PY=(/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python)
fi
# PYTHONPATH selects this branch's code over the editable install in main.
PYTHONPATH="$ROOT/.." "${BRUIT_PY[@]}" -m compy_pressure.microseism_check \
  --paired-run "$R/paired_spectra/run_bfbef17f2c13" \
  --reference-result "$R/rayleigh_calibration/run_392920872d9c/comparison.json" \
  --reference-sha256 0d6f7552d45fe1ad59194a6c45dd00501b3bf6ca6c54c8f8f0a324d93a5314e9 \
  --output-root "$OUTPUT" "$@"
