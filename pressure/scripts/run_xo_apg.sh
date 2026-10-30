#!/usr/bin/env bash
# Separate HDH APG / HHZ run; preserve the completed EDH/DPG results.
set -eo pipefail
cd "$(dirname "$0")/.."
BRUIT_PY=(conda run --no-capture-output -n bruitfm python)
RUN_ROOT=results/NETWORKS/XO_APG_NO_TILT
STATIONS=(--station LA33 --station LA34 --station LT02)
"${BRUIT_PY[@]}" scripts/prepare_xo_run.py \
  "${STATIONS[@]}" --pressure-kind apg --event-resample-rate 2 \
  --native-rate-tolerance-ppm 0.1 \
  --fetch-responses --exclude-incompatible-files --results-dir "$RUN_ROOT"
ARGS=(--network XO "${STATIONS[@]}" --data-dir "$RUN_ROOT/_inputs/waveforms" \
  --master-dir "$RUN_ROOT/_master" --results-dir "$RUN_ROOT")
# Exact response snapshots are prepared above; no PPSD computation is needed.
"${BRUIT_PY[@]}" -m compy_pressure.network rayleigh "${ARGS[@]}" \
  --rayleigh-catalog "$RUN_ROOT/_inputs/catalog_usgs_m5p5.json" \
  --rayleigh-dispersion-file ../compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --rayleigh-frequency-min 0.03 --rayleigh-frequency-max 0.07 \
  --rayleigh-pressure-sign 1 --rayleigh-tilt-method none \
  --rayleigh-event-resample-rate 2 --rayleigh-event-sampling all \
  --rayleigh-native-rate-tolerance-ppm 0.1 \
  --rayleigh-no-master-update
"${BRUIT_PY[@]}" -m compy_pressure.network view "${ARGS[@]}" --view-kind rayleigh \
  --view-workflow "$RUN_ROOT/_workflow/rayleigh.json"
