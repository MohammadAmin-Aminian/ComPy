#!/bin/bash
cd "$(dirname "$0")/.."
PY=/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
for s in RR28 RR29 RR34 RR36 RR38 RR40 RR50 RR52; do
  $PY scripts/a_invert.py $s vpvs > logs/a_vpvs_$s.log 2>&1; echo "$s exit $?"; tail -1 logs/a_vpvs_$s.log
done
echo ALL_DONE
