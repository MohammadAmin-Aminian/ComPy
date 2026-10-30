#!/usr/bin/env bash
# All downloaded XO stations; separate masters/results, no tilt or applied gains.
set -eo pipefail
cd "$(dirname "$0")/.."
BRUIT_PY=(conda run --no-capture-output -n bruitfm python)
RUN_ROOT=results/NETWORKS/XO_ALL_NO_TILT
"${BRUIT_PY[@]}" scripts/prepare_xo_run.py \
  --fetch-responses --exclude-incompatible-files --results-dir "$RUN_ROOT"
ARGS=(--network XO --data-dir "$RUN_ROOT/_inputs/waveforms" \
  --master-dir "$RUN_ROOT/_master" --results-dir "$RUN_ROOT")
"${BRUIT_PY[@]}" -m compy_pressure.network ppsd "${ARGS[@]}" --workers 4
"${BRUIT_PY[@]}" -m compy_pressure.network plots "${ARGS[@]}"
"${BRUIT_PY[@]}" -m compy_pressure.network rayleigh "${ARGS[@]}" \
  --rayleigh-catalog "$RUN_ROOT/_inputs/catalog_usgs_m5p5.json" \
  --rayleigh-dispersion-file ../compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --rayleigh-frequency-min 0.03 --rayleigh-frequency-max 0.07 \
  --rayleigh-pressure-sign 1 --rayleigh-tilt-method none \
  --rayleigh-event-sampling all --rayleigh-no-master-update
