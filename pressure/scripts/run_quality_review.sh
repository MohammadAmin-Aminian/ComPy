#!/usr/bin/env bash
# Read saved six-check evidence and write a new review + noise selection plan.
# No waveform calculation, approval, gain application or production master write.
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  printf 'Run with bash; do not source this script.\n' >&2
  return 2
fi
set -Eeuo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SOURCES="$ROOT/docs/quality_review_sources.json"
DECISIONS=
OUT=
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sources) SOURCES="$2"; shift 2 ;;
    --decisions) DECISIONS="$2"; shift 2 ;;
    --output-dir) OUT="$2"; shift 2 ;;
    --help|-h)
      printf 'Usage: bash scripts/run_quality_review.sh [--sources JSON] [--decisions JSON] [--output-dir NEW_DIR]\n'
      exit 0 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; exit 2 ;;
  esac
done
cd "$ROOT"
OUT="${OUT:-$ROOT/results/REVIEWS/SIX_CHECKS/run_$(date -u +%Y%m%dT%H%M%SZ)_$RANDOM}"
if command -v conda >/dev/null; then
  PY=(conda run --no-capture-output -n bruitfm python)
else
  PY=(/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python)
fi
export PYTHONPATH="$ROOT/..:$ROOT"
ARGS=(--sources "$SOURCES" --output-dir "$OUT")
[[ -z "$DECISIONS" ]] || ARGS+=(--decisions "$DECISIONS")
"${PY[@]}" -m compy_pressure.quality_review "${ARGS[@]}"
"${PY[@]}" -m compy_pressure.noise_model --review "$OUT/review.json" \
  --output-dir "$OUT/noise_plan" --plan
printf '\nSaved review: %s/station_summary.csv\n' "$OUT"
printf 'Exact channel review: %s/station_channel_review.csv\n' "$OUT"
printf 'No production master, score, response or waveform was changed.\n'
