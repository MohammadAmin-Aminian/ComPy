#!/usr/bin/env bash
# Run all three pressure-vertical calibration methods on all available data.
#
#   1. Rayleigh   (main checkout) saved results for all 89 station epochs;
#                 --rerun-rayleigh explicitly replays the saved event windows.
#   2. Microseism (main checkout) base week at all stations, repeated weeks
#                 (--new-weeks / --shifts) and the frozen frequency/week rule.
#   3. Local P    (this checkout, reviewed v2) all eligible
#                 pairs, export and station table.
#
# Every stage writes NEW folders only. No master, response, raw waveform,
# saved result or score is modified. Stages reuse each project's tested
# entry points; nothing here changes a scientific setting.
#
# Usage:  bash scripts/run_all_calibrations.sh [--check] [--rerun-rayleigh]
#           [--new-weeks] [--shifts "-12 -6 -2 -1 1 2 6 12"] [--skip-weeks]
#           [--main DIR] [--localp DIR]
# Default: reuse completed Rayleigh, no new week shifts. --new-weeks can take
# several hours with eight shifts; --rerun-rayleigh explicitly replays trials.
# This covers configured stations/selected windows, not every archive day.
# Run with bash, never source it.
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  printf 'Run this file with bash; do not source it.\n' >&2
  return 2
fi
set -Eeuo pipefail

MAIN="$HOME/Projects/bruit-fm"
LOCALP="$(cd "$(dirname "$0")/.." && pwd)"
CHECK_ONLY=0
RERUN_RAYLEIGH=0
SKIP_WEEKS=1
# 7D has only ~5-10 weeks on disk, so short shifts are needed there; long
# shifts give independent weeks for YV/YL/XO. Each shift is a full new run.
SHIFTS="-12 -6 -2 -1 1 2 6 12"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --check) CHECK_ONLY=1; shift ;;
    --reuse-rayleigh) RERUN_RAYLEIGH=0; shift ;;
    --rerun-rayleigh) RERUN_RAYLEIGH=1; shift ;;
    --new-weeks) SKIP_WEEKS=0; shift ;;
    --skip-weeks) SKIP_WEEKS=1; shift ;;
    --shifts) [[ $# -ge 2 ]] || exit 2; SHIFTS="$2"; shift 2 ;;
    --main) [[ $# -ge 2 ]] || exit 2; MAIN="$2"; shift 2 ;;
    --localp) [[ $# -ge 2 ]] || exit 2; LOCALP="$2"; shift 2 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; exit 2 ;;
  esac
done

STAGE=preflight
SESSION=
trap 'STATUS=$?; printf "\nStopped during %s (exit %s, line %s). Logs: %s\n" "$STAGE" "$STATUS" "$LINENO" "${SESSION:-not started}" >&2; exit "$STATUS"' ERR

if command -v conda >/dev/null; then
  PY=(conda run --no-capture-output -n bruitfm python)
else
  PY=(/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python)
fi

for FILE in \
  "$MAIN/scripts/run_three_methods.sh" \
  "$MAIN/scripts/run_microseism_weeks.sh" \
  "$MAIN/scripts/combine_microseism_weeks.py" \
  "$LOCALP/scripts/run_local_p_campaign.py" \
  "$LOCALP/scripts/export_local_p_results.py" \
  "$LOCALP/scripts/summarize_local_p_campaign.py" \
  "$LOCALP/docs/local_p_timing.json"; do
  [[ -f "$FILE" ]] || { printf 'Required file missing: %s\n' "$FILE" >&2; exit 2; }
done
for CATALOG in \
  usgs_local_XO_20180520_20190101_m2p5.csv \
  usgs_local_7D_20111019_20120101_m2p0.csv \
  usgs_local_YV_20121231_20131127_m2p5.csv; do
  [[ -f "$LOCALP/results/LOCAL_P/catalogs/$CATALOG" ]] || {
    printf 'Local-P catalogue missing: %s (see docs/local_p_check.md section 3.4)\n' \
      "$CATALOG" >&2
    exit 2
  }
done
for SHIFT in $SHIFTS; do
  [[ "$SHIFT" =~ ^-?[0-9]+$ && "$SHIFT" -ne 0 ]] || {
    printf 'Week shifts must be nonzero integers: %s\n' "$SHIFT" >&2; exit 2; }
done

printf 'Main checkout:    %s  (%s, %s)\n' "$MAIN" \
  "$(git -C "$MAIN" rev-parse --short HEAD)" "$(git -C "$MAIN" branch --show-current)"
printf 'Local-P checkout: %s  (%s, %s)\n' "$LOCALP" \
  "$(git -C "$LOCALP" rev-parse --short HEAD)" "$(git -C "$LOCALP" branch --show-current)"
for REPO in "$MAIN" "$LOCALP"; do
  if [[ -n "$(git -C "$REPO" status --porcelain -- ../compy_pressure scripts)" ]]; then
    printf 'Warning: uncommitted code in %s; provenance will record a dirty tree.\n' "$REPO"
  fi
done

STAGE=preflight-rayleigh-microseism
bash "$MAIN/scripts/run_three_methods.sh" --check
STAGE=preflight-local-p
(cd "$LOCALP" && PYTHONPATH="$LOCALP/..:$LOCALP" "${PY[@]}" scripts/run_local_p_campaign.py \
  --campaign-dir results/LOCAL_P/campaigns/preflight --dry-run | tail -5)
if [[ "$CHECK_ONLY" -eq 1 ]]; then
  printf '\nPreflight passed. Nothing was run.\n'
  exit 0
fi

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
SESSION="$MAIN/results/EXPERIMENTS/ALL_CALIBRATIONS/run_${STAMP}"
mkdir -p "$(dirname "$SESSION")"
mkdir "$SESSION"
LOCALP_CAMPAIGN="$LOCALP/results/LOCAL_P/campaigns/all_${STAMP}"
printf 'Session: %s\n' "$SESSION" | tee "$SESSION/session.log"

# 1-2. Rayleigh and microseism base week (fresh session under the main checkout).
STAGE=rayleigh-and-microseism
THREE_ARGS=()
[[ "$RERUN_RAYLEIGH" -eq 1 ]] && THREE_ARGS+=(--rerun-rayleigh)
bash "$MAIN/scripts/run_three_methods.sh" ${THREE_ARGS[@]+"${THREE_ARGS[@]}"} 2>&1 \
  | tee "$SESSION/1_rayleigh_microseism.log"
RAYLEIGH="$(sed -n 's/^Rayleigh: //p' "$SESSION/1_rayleigh_microseism.log" | tail -n 1)"
COMPARISON="$(sed -n 's/^Comparison: //p' "$SESSION/1_rayleigh_microseism.log" | tail -n 1)"
[[ -n "$RAYLEIGH" && -f "$COMPARISON" ]]
MICRO_BASE="$(dirname "$COMPARISON")"

# 2b. Repeated microseism weeks and the frozen frequency/week rule.
FROZEN="$SESSION/microseism_frozen_rule"
if [[ "$SKIP_WEEKS" -eq 0 ]]; then
  STAGE=microseism-weeks
  # shellcheck disable=SC2086  # intentional word splitting of the shift list
  bash "$MAIN/scripts/run_microseism_weeks.sh" "$MICRO_BASE" $SHIFTS 2>&1 \
    | tee "$SESSION/2_microseism_weeks.log"
  STAGE=microseism-frozen-rule
  (cd "$MAIN" && PYTHONPATH="$MAIN/..:$MAIN" "${PY[@]}" \
    scripts/combine_microseism_weeks.py "$MICRO_BASE" --output-dir "$FROZEN") 2>&1 \
    | tee "$SESSION/3_microseism_frozen_rule.log"
else
  FROZEN="skipped (--skip-weeks)"
fi

# 3. Local P (writes only under the local-P checkout's results/LOCAL_P/).
STAGE=local-p-campaign
(cd "$LOCALP" && PYTHONPATH="$LOCALP/..:$LOCALP" "${PY[@]}" scripts/run_local_p_campaign.py \
  --campaign-dir "$LOCALP_CAMPAIGN") 2>&1 | tee "$SESSION/4_local_p.log"
STAGE=local-p-export
(cd "$LOCALP" && PYTHONPATH="$LOCALP/..:$LOCALP" "${PY[@]}" scripts/export_local_p_results.py \
  --campaign-dir "$LOCALP_CAMPAIGN" --master "$MAIN/exports/master_inventory.csv" \
  --output-dir "${LOCALP_CAMPAIGN}_export") 2>&1 | tee "$SESSION/5_local_p_export.log"
STAGE=local-p-table
(cd "$LOCALP" && PYTHONPATH="$LOCALP/..:$LOCALP" "${PY[@]}" scripts/summarize_local_p_campaign.py \
  --campaign-dir "$LOCALP_CAMPAIGN" \
  --csv "$LOCALP_CAMPAIGN/station_table.csv" \
  --markdown "$LOCALP_CAMPAIGN/station_table.md") 2>&1 | tee "$SESSION/6_local_p_table.log"

STAGE=finished
{
  printf '\nFinished %s\n' "$(date -u +%FT%TZ)"
  printf 'Rayleigh experiment:        %s\n' "$RAYLEIGH"
  printf 'Microseism base week:       %s\n' "$MICRO_BASE"
  printf 'Microseism frozen rule:     %s\n' "$FROZEN"
  printf 'Local-P campaign:           %s\n' "$LOCALP_CAMPAIGN"
  printf 'Local-P table:              %s\n' "$LOCALP_CAMPAIGN/station_table.md"
  printf 'Logs:                       %s\n' "$SESSION"
  printf 'No master, response, saved result or score was modified.\n'
} | tee -a "$SESSION/session.log"
