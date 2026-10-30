# Baseline microseism pressure–vertical comparison

Explicit campaign: `/Users/mohammadamin/Projects/bruit-fm/results/DIAGNOSTICS/MICROSEISM_BASELINE_LOGS/run_20261014T170400Z/summary.tsv`.

## Result

All 46 frozen-baseline stations completed. The saved P/Z diagnostics provide coherence, phase and transfer-amplitude checks, but **no independent microseism pressure gain**: `H_nominal=(g_pressure/g_vertical)*eta(f; water, structure)`. The Rayleigh multiplier is fixed and only changes the displayed transfer through `H_fixed=H_nominal/g`.

Across the campaign, 37/46 stations have weekly median coherence² at least 0.8, and 30/46 contain at least one 0.1–0.2 Hz bin above the Ruan et al. 0.95 reference. The later Saikia daily reference selects 149/275 station-days across 36/46 stations. These are descriptive comparisons, not BRUIT-FM pass scores.

Highest weekly medians: YL.C05W (0.982), YL.A02W (0.979), YL.B11W (0.975), 7D.J45A (0.975), YL.C17W (0.972). Lowest: YV.RR11 (0.583), YV.RR50 (0.630), YV.RR29 (0.644), YV.RR07 (0.687), YL.C11W (0.693).

RR52, the pilot reference, is close to the YV network median: coherence² 0.802, 29/361 Ruan-reference bins, 6/7 selected days and phase +4.1°. Its nominal/fixed display amplitudes remain 4.986e-08/7.346e-08 m/Pa.

## Network comparison

| Network | Stations | Segments min/median/max | Median coherence² (network median) | Weekly median ≥0.8 | Any bins >0.95 | Selected days | Phase center | Nominal → fixed amplitude spread |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| YV | 16 | 143/167/168 | 0.819 | 11 | 4 | 31/111 | -3.5° | ×1.579 → ×1.342 |
| YL | 22 | 24/132/167 | 0.953 | 20 | 20 | 81/109 | -175.8° | ×1.629 → ×1.587 |
| 7D | 4 | 98/155/167 | 0.939 | 4 | 4 | 27/27 | +177.3° | ×1.292 → ×1.362 |
| XO | 4 | 154/161/168 | 0.799 | 2 | 2 | 10/28 | +175.0° | ×1.484 → ×1.459 |

The phase separates into 16 near-0° and 30 near-180° stations. The split is by network: positive scalar gains cannot cause it. Treat it as an unresolved pressure/vertical sign or response-convention difference until metadata and processing conventions are audited; do not pool signed admittance curves yet.

The nominal-to-fixed spread factors are descriptive only. A lower cross-station spread after displaying the frozen gains is not an independent validation because water depth, seafloor structure, vertical response and wavefield differ among stations.

## First method and papers

The already-frozen first-method values retain their prior paper comparisons: 8 YV values differ from [Aminian et al. (2025)](https://www.ipgp.fr/~stutz/2025-aminian-crawford-stutzmann-al-GJI.pdf) by +0.35% to +6.08% (median +3.20%) in the same `g` direction. The 21 YL peers excluding C07W differ from [Zha & Webb (2016)](https://doi.org/10.1002/2015JB012688) by -7.04% to -3.37% (median -5.43%) in reciprocal `1/g`; C07W remains the known mixed-state outlier at +86.76%. These numbers come from method 1, not from a microseism fit.

[Ruan et al. (2014)](https://doi.org/10.1002/2014JB011162) is represented only by its >0.95 coherence reference: reproducing its 7D scale factors requires the coupled modal structure/water model and joint inversion. [Deng et al. (2022)](https://doi.org/10.1029/2021GL096952) uses frequency-dependent body-wave bands for XO, so no scalar microseism/paper difference is reported. The [Saikia et al. (2020)](https://doi.org/10.1007/s11001-020-09402-0) daily >0.8 rule is shown only as a selection reference.

Primary-source interpretation and links are retained in [the RR52 method note](rr52_microseism_diagnostic.md); the frozen publication-direction rules are in [the campaign note](microseism_baseline_run.md).

## Limits and preservation

- Campaign intervals were selected independently and are not synchronized; network/source-season differences remain.
- Coverage ranges from 24 to 168 hourly segments, so short records should not be ranked as equivalent week-long estimates.
- The diagnostics contain no horizontal/tilt correction, local-event/coda screen, modal prediction, gain uncertainty or classification.
- Every explicit log, diagnostic, saved artifact and stored aggregate was verified. No waveform, response, master, prior result or score was changed; no fixed gain was applied to original data.

Machine-readable station values: [CSV](microseism_baseline_comparison.csv) and [JSON](microseism_baseline_comparison.json). [Comparison figure](figures/microseism_baseline_comparison.png).
