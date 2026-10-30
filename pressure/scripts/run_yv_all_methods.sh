#!/usr/bin/env bash
# Calibrate every YV station with each reviewed method, using this checkout's code.
#
#   1. Rayleigh P/Z   `python -m compy_pressure.network rayleigh` with the frozen YV
#                     no-tilt settings (0.03-0.07 Hz, sign -1, tilt none, 3-8
#                     events, saved M5.5 catalogue, RR52 dispersion) into a NEW
#                     workflow folder; the master is never updated.
#                     --rayleigh-events bounded (default) reproduces the frozen
#                     2026-10-09 event selection and must return RR52
#                     g = 0.6787278140712512 exactly, or the run stops.
#                     --rayleigh-events all evaluates every eligible event (the
#                     package default since a19bc3b); a new result, not a reproduction.
#   2. Microseism P/Z the common 2013-06-10 base week at every YV station;
#                     --weeks adds shifted weeks and the frozen frequency/week rule.
#   3. Local P        reviewed v2 campaign for every eligible YV pair, then the
#                     export and station table.
#   (opt) Pole trial  --pole-trial replays the opt-in pressure-pole-reflection
#                     experiment on the preserved YV workflow (saved spectra).
#
# Every stage writes NEW folders only. No raw waveform, response, master,
# saved result or score is modified and no gain is applied. Outputs are
# provisional method results, not accepted calibrations (docs/deployment.md).
#
# Usage: bash scripts/run_yv_all_methods.sh [--check] [--only STAGE]
#          [--weeks] [--shifts "-12 -6 -2 -1 1 2 6 12"] [--pole-trial]
#   --check        verify inputs, print the plan and the local-P dry run; run nothing
#   --only STAGE   rayleigh | microseism | local-p (repeatable)
#   --rayleigh-events bounded|all   see above (default bounded)
#   --allow-dirty  run with uncommitted code in ../compy_pressure/ or scripts/ (recorded);
#                  by default the run refuses, and it always stops if the code
#                  changes while it runs.
# Set BRUITFM_PYTHON to use another interpreter and BRUITFM_SOURCE_ROOT to run
# another checkout (default: this one). Run with bash; never source it.
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  printf 'Run this file with bash; do not source it.\n' >&2
  return 2
fi
set -Eeuo pipefail

# BRUITFM_SOURCE_ROOT selects the checkout whose code, exports and results are used.
ROOT="${BRUITFM_SOURCE_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
PY="${BRUITFM_PYTHON:-/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python}"
CHECK_ONLY=0
WEEKS=0
POLE_TRIAL=0
ALLOW_DIRTY=0
RAYLEIGH_EVENTS=bounded
RR52_FROZEN_G=0.6787278140712512
SHIFTS="-12 -6 -2 -1 1 2 6 12"
ONLY=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --check) CHECK_ONLY=1; shift ;;
    --weeks) WEEKS=1; shift ;;
    --pole-trial) POLE_TRIAL=1; shift ;;
    --allow-dirty) ALLOW_DIRTY=1; shift ;;
    --rayleigh-events)
      [[ $# -ge 2 && "$2" =~ ^(bounded|all)$ ]] || {
        printf -- '--rayleigh-events needs bounded or all\n' >&2; exit 2; }
      RAYLEIGH_EVENTS="$2"; shift 2 ;;
    --shifts) [[ $# -ge 2 ]] || exit 2; SHIFTS="$2"; shift 2 ;;
    --only)
      [[ $# -ge 2 && "$2" =~ ^(rayleigh|microseism|local-p)$ ]] || {
        printf -- '--only needs rayleigh, microseism or local-p\n' >&2; exit 2; }
      ONLY+=("$2"); shift 2 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; exit 2 ;;
  esac
done
wanted() {
  [[ ${#ONLY[@]} -eq 0 ]] && return 0
  local s
  for s in "${ONLY[@]}"; do [[ "$s" == "$1" ]] && return 0; done
  return 1
}
STAGE=preflight
SESSION=
trap 'STATUS=$?; printf "\nStopped during %s (exit %s, line %s). Logs: %s\n" "$STAGE" "$STATUS" "$LINENO" "${SESSION:-not started}" >&2; exit "$STATUS"' ERR
cd "$ROOT"
export PYTHONPATH="$ROOT/..:$ROOT"

CATALOG=results/RAYLEIGH_REVIEW/catalog_usgs_m5p5_2012_2013.json
DISPERSION=../compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv
LOCALP_CATALOG=results/LOCAL_P/catalogs/usgs_local_YV_20121231_20131127_m2p5.csv
command -v "$PY" >/dev/null || { printf 'Python not found: %s\n' "$PY" >&2; exit 2; }
for FILE in exports/master_inventory.csv "$CATALOG" "$DISPERSION" \
  "$LOCALP_CATALOG" docs/local_p_timing.json \
  scripts/run_microseism_check.py scripts/run_microseism_weeks.sh \
  scripts/combine_microseism_weeks.py scripts/run_local_p_campaign.py \
  scripts/export_local_p_results.py scripts/summarize_local_p_campaign.py \
  scripts/run_calibration_experiment.py; do
  [[ -e "$FILE" ]] || { printf 'Required input missing: %s/%s\n' "$ROOT" "$FILE" >&2; exit 2; }
done
for SHIFT in $SHIFTS; do
  [[ "$SHIFT" =~ ^-?[0-9]+$ && "$SHIFT" -ne 0 ]] || {
    printf 'Week shifts must be nonzero integers: %s\n' "$SHIFT" >&2; exit 2; }
done
"$PY" -c 'import compy_pressure.network, compy_pressure.rayleigh_network, compy_pressure.microseism_check, compy_pressure.local_p_check, scripts.run_microseism_check; print("Imports: OK")'
# The executed code is HEAD plus any uncommitted src/scripts changes. Python
# loads modules at process start, so an edit during a run mixes code versions.
code_state() {
  {
    git rev-parse HEAD
    git status --porcelain --untracked-files=all -- ../compy_pressure scripts
    git diff HEAD -- ../compy_pressure scripts
  } | shasum -a 256 | cut -c1-16
}
CODE_STATE="$(code_state)"
DIRTY="$(git status --porcelain --untracked-files=all -- ../compy_pressure scripts)"
printf 'Code: %s  (%s, %s, state %s)\n' "$ROOT" "$(git rev-parse --short HEAD)" \
  "$(git branch --show-current)" "$CODE_STATE"
if [[ -n "$DIRTY" ]]; then
  printf 'Uncommitted code in ../compy_pressure/ or scripts/:\n%s\n' "$DIRTY"
  if [[ "$ALLOW_DIRTY" -eq 0 ]]; then
    printf 'Refusing to run: commit, stash or review it first, or pass --allow-dirty.\n' >&2
    exit 2
  fi
fi
same_code() {
  [[ "$(code_state)" == "$CODE_STATE" ]] || {
    printf 'Code in ../compy_pressure/ or scripts/ changed during the run; stopping before %s.\n' \
      "$1" >&2
    exit 3
  }
}

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
RAYLEIGH_DIR="results/NETWORKS/YV_RAYLEIGH_NO_TILT_${RAYLEIGH_EVENTS}_$STAMP"
LOCALP_CAMPAIGN="results/LOCAL_P/campaigns/yv_$STAMP"
RAYLEIGH_CMD=("$PY" -m compy_pressure.network rayleigh --network YV
  --results-dir "$RAYLEIGH_DIR"
  --rayleigh-catalog "$CATALOG" --rayleigh-dispersion-file "$DISPERSION"
  --rayleigh-frequency-min 0.03 --rayleigh-frequency-max 0.07
  --rayleigh-pressure-sign -1 --rayleigh-tilt-method none
  --rayleigh-min-events 3 --rayleigh-max-events 8 --rayleigh-event-sampling "$RAYLEIGH_EVENTS"
  --rayleigh-no-master-update)
MICRO_CMD=("$PY" scripts/run_microseism_check.py --network YV --start YV=2013-06-10
  --output-root results/EXPERIMENTS/MICROSEISM_YV)

STAGE=preflight-local-p
printf '\nLocal-P plan (dry run, nothing written):\n'
"$PY" scripts/run_local_p_campaign.py --campaign-dir "$LOCALP_CAMPAIGN" \
  --networks YV --dry-run | tail -3
if [[ "$CHECK_ONLY" -eq 1 ]]; then
  printf '\nPlanned stages:\n'
  wanted rayleigh && printf '  Rayleigh:   %s\n' "${RAYLEIGH_CMD[*]}"
  wanted microseism && printf '  Microseism: %s\n' "${MICRO_CMD[*]}"
  wanted microseism && [[ "$WEEKS" -eq 1 ]] && printf '  Weeks:      %s\n' "$SHIFTS"
  wanted local-p && printf '  Local P:    %s\n' "$LOCALP_CAMPAIGN"
  [[ "$POLE_TRIAL" -eq 1 ]] && printf '  Pole trial: run_calibration_experiment.py --network YV --all-stations\n'
  printf '\nPreflight passed. Nothing was run.\n'
  exit 0
fi

SESSION="results/EXPERIMENTS/YV_ALL_METHODS/run_$STAMP"
mkdir -p "$(dirname "$SESSION")"
mkdir "$SESSION"
{
  printf 'Session: %s/%s\n' "$ROOT" "$SESSION"
  printf 'Code: %s state %s; uncommitted: %s\n' "$(git rev-parse HEAD)" "$CODE_STATE" \
    "${DIRTY:-none}"
  printf 'Rayleigh events: %s\n' "$RAYLEIGH_EVENTS"
} | tee "$SESSION/session.log"

RAYLEIGH_OUT="not run"
if wanted rayleigh; then
  STAGE=rayleigh
  same_code rayleigh
  [[ ! -e "$RAYLEIGH_DIR" ]] || { printf 'Exists: %s\n' "$RAYLEIGH_DIR" >&2; exit 2; }
  "${RAYLEIGH_CMD[@]}" 2>&1 | tee "$SESSION/1_rayleigh.log"
  RAYLEIGH_OUT="$RAYLEIGH_DIR/_workflow/rayleigh.json"
  same_code rayleigh-check
  RR52_G="$("$PY" -c 'import json, sys
entries = json.load(open(sys.argv[1]))["entries"].values()
print(next(e["update"]["rayleigh_pressure_multiplier"] for e in entries
           if e["seed_id"] == "YV.RR52.00.BDH"))' "$RAYLEIGH_OUT")"
  printf 'RR52 g = %s (frozen %s)\n' "$RR52_G" "$RR52_FROZEN_G" | tee -a "$SESSION/session.log"
  if [[ "$RAYLEIGH_EVENTS" == bounded && "$RR52_G" != "$RR52_FROZEN_G" ]]; then
    printf 'RR52 does not reproduce the frozen reference; stopping.\n' \
      | tee -a "$SESSION/session.log" >&2
    exit 3
  fi
fi

POLE_OUT="not run"
if [[ "$POLE_TRIAL" -eq 1 ]]; then
  STAGE=pole-trial
  same_code pole-trial
  POLE_OUT="results/EXPERIMENTS/SOCRATES/yv_pole_trial_$STAMP"
  "$PY" scripts/run_calibration_experiment.py --network YV --all-stations \
    --output-dir "$POLE_OUT" 2>&1 | tee "$SESSION/1b_pole_trial.log"
fi

MICRO_BASE="not run"
FROZEN="not run"
if wanted microseism; then
  STAGE=microseism
  same_code microseism
  "${MICRO_CMD[@]}" 2>&1 | tee "$SESSION/2_microseism.log"
  MICRO_BASE="$(sed -n 's/^{.*}: //p' "$SESSION/2_microseism.log" | tail -n 1)"
  [[ -d "$MICRO_BASE" ]] || { printf 'Microseism run folder not found in log\n' >&2; exit 2; }
  if [[ "$WEEKS" -eq 1 ]]; then
    STAGE=microseism-weeks
    # shellcheck disable=SC2086  # intentional word splitting of the shift list
    bash scripts/run_microseism_weeks.sh "$MICRO_BASE" $SHIFTS 2>&1 \
      | tee "$SESSION/3_microseism_weeks.log"
    STAGE=microseism-frozen-rule
    FROZEN="$SESSION/microseism_frozen_rule"
    "$PY" scripts/combine_microseism_weeks.py "$MICRO_BASE" --output-dir "$FROZEN" 2>&1 \
      | tee "$SESSION/4_microseism_frozen_rule.log"
  fi
fi

LOCALP_TABLE="not run"
if wanted local-p; then
  STAGE=local-p-campaign
  same_code local-p
  "$PY" scripts/run_local_p_campaign.py --campaign-dir "$LOCALP_CAMPAIGN" \
    --networks YV 2>&1 | tee "$SESSION/5_local_p.log"
  STAGE=local-p-export
  "$PY" scripts/export_local_p_results.py --campaign-dir "$LOCALP_CAMPAIGN" \
    --master exports/master_inventory.csv --output-dir "${LOCALP_CAMPAIGN}_export" 2>&1 \
    | tee "$SESSION/6_local_p_export.log"
  STAGE=local-p-table
  LOCALP_TABLE="$LOCALP_CAMPAIGN/station_table.md"
  "$PY" scripts/summarize_local_p_campaign.py --campaign-dir "$LOCALP_CAMPAIGN" \
    --csv "$LOCALP_CAMPAIGN/station_table.csv" --markdown "$LOCALP_TABLE" 2>&1 \
    | tee "$SESSION/7_local_p_table.log"
fi

STAGE=finished
same_code finish
{
  printf '\nFinished %s\n' "$(date -u +%FT%TZ)"
  printf 'Rayleigh workflow:      %s\n' "$RAYLEIGH_OUT"
  printf 'Pole-reflection trial:  %s\n' "$POLE_OUT"
  printf 'Microseism base week:   %s\n' "$MICRO_BASE"
  printf 'Microseism frozen rule: %s\n' "$FROZEN"
  printf 'Local-P station table:  %s\n' "$LOCALP_TABLE"
  printf 'No master, response, saved result or score was modified.\n'
} | tee -a "$SESSION/session.log"
