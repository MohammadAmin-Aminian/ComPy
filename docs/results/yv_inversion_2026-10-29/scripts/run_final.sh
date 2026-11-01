#!/bin/bash
# Gain-sampled, joint-epoch, free-Vp/Vs runs (and the H2 bracket).
cd "$(dirname "$0")/.."
PY=/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
for s in RR50 RR38 RR36 RR40 RR28 RR29 RR34 RR52; do
  $PY scripts/a_invert.py $s gain final final_stable final_h2 > logs/a_final_$s.log 2>&1
  echo "$s exit $?"; grep -E 'chi2' logs/a_final_$s.log
done
echo ALL_DONE
