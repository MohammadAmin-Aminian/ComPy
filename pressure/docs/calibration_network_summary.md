# Calibration checkpoint across four networks

Checkpoint: 12 October 2026. This report reads completed saved workflows only.
No earthquake batch, spectrum or gain was recomputed. Repeated runs of the same
station are counted once. All estimates are provisional; no gains are applied.

## Network coverage

| Network | Distinct stations | Provisional gains | Unresolved | Reference-band exclusions |
|---|---:|---:|---:|---:|
| YV | 43 | 16 | 25 | 2 |
| YL | 25 | 22 | 3 | 0 |
| 7D | 15 | 4 | 11 | 0 |
| XO | 6 | 4 | 2 | 0 |
| **Total** | **89** | **46** | **41** | **2** |

YV uses the completed separate no-tilt workflow, YL the all-station no-tilt
workflow, 7D its all-station workflow, and XO its latest APG/DPG gauge workflow.
The original YV and earlier XO runs remain preserved. Counts describe processing
outcomes, not instrument pass/fail classifications or whole-project completion.

## Completed XO gauge results

The original six-station XO run used EDH hydrophones at LA33/LA34/LT02 and yielded
two DPG estimates. The separate HDH APG run at 2 Hz completed with two estimates
and one unresolved station. The latest all-station gauge run uses APGs where
hardware metadata identifies APG, DPGs otherwise, and 1 Hz working copies:

| Station | Sensor / pair | Pressure multiplier g | Sensitivity multiplier 1/g | Accepted events |
|---|---|---:|---:|---:|
| LA33 | APG HDH / HHZ | 0.934275 | 1.070348 | 7 |
| LA34 | APG HDH / HHZ | 0.962680 | 1.038767 | 8 |
| LT02 | APG HDH / HHZ | — | — | 0 |
| WD48 | DPG LDH / LHZ | — | — | 0 |
| WD66 | DPG LDH / LHZ | 0.909498 | 1.099508 | 16 |
| WD70 | DPG LDH / LHZ | 0.929889 | 1.075398 | 21 |

`g` multiplies nominal pressure; `1/g` multiplies counts/Pa sensitivity. Pressure
polarity is separately assumed. Scalars do not validate the response outside the
fit band. Exact saved values, rejection counts and source hashes are in the
[numerical checkpoint](calibration_network_summary.json).

| Station | Comparison baseline | Change in g | Accepted events, before → after |
|---|---|---:|---:|
| LA33 | Separate 2 Hz APG run | +0.000025% | 7 → 7 |
| LA34 | Separate 2 Hz APG run | +0.000013% | 8 → 8 |
| WD66 | Original native-grid DPG run | +0.199244% | 10 → 16 |
| WD70 | Original native-grid DPG run | −0.007922% | 16 → 21 |

All previously accepted events in these comparisons are retained. WD66 adds six
and WD70 five events. Numerical stability between processing variants is not an
external accuracy measurement. LT02 rejects 40 attempts for band coherence;
separate logger timing remains unverified. WD48 has 29 missing windows and seven
coherence rejections. Clustering and window guards also exclude events. Do not
label either sensor faulty or relax QC merely to obtain a gain.

## Comparison with published methods

Our current Rayleigh method removes nominal complex responses, converts vertical
velocity to upward acceleration, measures P/a in event windows, fits a scalar
pressure multiplier to a prescribed water-column/dispersion prediction, and takes
the median accepted-event gain. The no-tilt runs use P/Z only. Event sampling uses
all eligible catalogue events under existing magnitude/depth/clustering/QC rules.

- **YL / Zha & Webb:** 21 of 22 reciprocal comparisons are 3.4–7.0% below the
  publication, median −5.4%; C07W is +86.8%. This outlier remains part of any
  assessment of stations with estimates. Our 0.03–0.07 Hz band, borrowed RR52
  dispersion, accepted M≥6 events and median estimator differ from their regional
  model, 0.04–0.1 Hz procedure and M≥7 selection. [Full station comparison](yl_all_zha_comparison.md)
  and [primary paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012688).
- **7D / Ruan et al.:** the paper fits a DPG amplitude scale factor jointly with
  sediment structure using D/P ratios. The original numerical station factors
  have not been located in the reviewed paper, author page, related studies or
  indexed public repository searches. This does not imply they were not computed
  or that no public copy exists. Table 3 contains sediment delays, thicknesses and
  depths. The 2015 follow-up corrects the 2014 *starting* profile parameter c to
  0.029, not 0.29; it is not a gain correction. [Source-search notes](7d_xo_runs.md)
  and [primary paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014JB011162).
- **XO / Deng et al.:** calibration uses body waves and separate frequency bands,
  rather than our scalar Rayleigh fit. The text reports WD48 about 0.95 and WD70
  0.91–1.71 across bands. Exact matched coefficients and their correction direction
  remain unverified; supplementary access was blocked. Do not infer agreement from
  a scalar falling inside a frequency-dependent range. [Primary paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021GL096952).
- **YV:** the subsequent Socrates review verified eight station-matched pressure
  multipliers in Aminian et al. (2025) Table 1. Compared with our **original tilt
  YV baseline**, all eight differ by +0.3–6.1% (median +3.0%). This is a shared
  method-lineage reproducibility comparison; event, glitch and tilt processing
  differ. It does not establish absolute accuracy or comparison of the latest
  no-tilt gains. [Review and qualifications](reviews/codex_integration_review.md)
  and [primary Table 1](https://www.ipgp.fr/~stutz/2025-aminian-crawford-stutzmann-al-GJI.pdf).

There is no measured overall accuracy percentage. The conversation's informal
7/10 assessment was a judgment of research-prototype maturity, not a scientific
score or acceptance threshold.

## Verification and next work

The reporting checkpoint had **143 tests passed**, Ruff check/format and runner syntax passed.
The reporting check rederived accepted-event medians and reciprocals from saved
JSONs, verified counts, and rechecked workflow/result hashes after reading them.
The original master and both new XO masters have no nonblank scientific scores;
their hashes stayed unchanged during reporting. This is a bounded saved-result
check, not a repeat of the comprehensive raw/archive audit.

Next: review C07W's matched event responses/windows; inspect LT02 timing and WD48
coverage; obtain Ruan's fitted scale factors and Deng's coefficient convention and
band values; validate local dispersion/depth/polarity and held-out uncertainty.
Preserve the complete YV batch and all raw data, inventories, masters and scores.


## Subsequent independent-branch integration review

The reviewed integration passes **155 tests**. All 416 saved accepted-event fits
reproduce exactly, 46 station medians remain unchanged, and bounded raw probes
confirm candidate gains and rejection diagnoses. C07W's relative vertical states
now have supporting evidence; its mixed-state median remains unsuitable as an
accepted DPG calibration. B01W is phase-limited rather than solvable by widening
scalar bounds alone; RR55 requires a separate response-identity review.
[Evidence, preservation scope and next task](reviews/codex_integration_review.md).
This qualification leaves the original numerical checkpoint JSON unchanged.
