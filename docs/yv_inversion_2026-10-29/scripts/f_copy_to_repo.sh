#!/bin/bash
# Copy scripts, figures, tables and JSON summaries (not the npz draws) into the ComPy repository.
cd "$(dirname "$0")/.."
DEST=../ComPy/docs/yv_inversion_2026-10-29
mkdir -p $DEST/scripts $DEST/figures $DEST/results/tables
cp scripts/*.py scripts/*.sh $DEST/scripts/
cp figures/*.png $DEST/figures/
cp results/tables/* $DEST/results/tables/
for s in RR28 RR29 RR34 RR36 RR38 RR40 RR50 RR52; do
  mkdir -p $DEST/results/$s && cp results/$s/*.json $DEST/results/$s/
done
du -sh $DEST
