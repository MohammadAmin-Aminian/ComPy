# Completed independent-week microseism evaluation

15 October 2026. Analysis code: main `cc323b9`. Frozen rule:
`claude-frequency-v1-20261015`, specified in Claude's `b603948` handoff.
No threshold, model, response or Rayleigh reference was changed after seeing
these new weeks. YL remains the development set; YV, 7D and XO are subsequent
temporal tests under the same model assumptions.

## Execution and preservation

All eight requested new shifts completed with exit code zero:

| Network | Shifts (weeks) | Station-weeks attempted | Paired / checked | Unresolved pairing | No base interval |
|---|---|---:|---:|---:|---:|
| 7D | -2, -1, +1, +2 | 60 | 44 | 16 | 0 |
| YV | -12, -6, +6, +12 | 172 | 144 | 24 | 4 |
| Total | | 232 | 188 | 40 | 4 |

All 188 new paired cases have the model-free diagnostic, with zero diagnostic
errors. Unresolved cases retain their original pairing/coverage reasons. YL
and XO were read from their completed runs; no calibration batch was rerun.
The final two YV shifts were scheduled concurrently in separate output folders;
the sequential launcher reused the completed +12 run instead of duplicating it.

Every new run records unchanged masters and no applied gain. The final read-only
preservation check verified 5,274 distinct protected files across the prior
three-method experiment and this review: zero missing or changed hashes.
This is the recorded-input scope, not an archive-wide byte audit.

The combined review now contains 24 run folders, including the four old 7D
shifts that were out of coverage. Those contribute no passing weeks. Four
identical YL repeats are deduplicated. The current review hashes 2,422 selected
input files; exact paths and hashes are in [microseism_weeks.json](microseism_weeks.json).
The full local review is
`results/REVIEWS/MICROSEISM_INTEGRATION/run_20261015T141224Z_2a1d41d2/`.

## Frozen-rule outcome

| Network | Stable (<=5%) | Moderate (5-10%) | Fewer than 3 passing weeks | No passing weeks |
|---|---:|---:|---:|---:|
| YV | 2 | 0 | 12 | 29 |
| YL (development) | 2 | 0 | 16 | 7 |
| 7D | 2 | 1 | 4 | 8 |
| XO | 0 | 0 | 3 | 3 |
| Total | 6 | 1 | 35 | 47 |

Seven stations meet the conditional repeatability rule. All seven have stable
relative polarity in their passing weeks and non-overlapping intervals.
None is an accepted absolute pressure calibration.

| Station | Passing / checked weeks | Conditional median g | Weekly spread | Status | Rayleigh baseline g | Difference from baseline |
|---|---:|---:|---:|---|---:|---:|
| YV.RR40 | 4 / 5 | 0.74853 | 3.15% | stable | 0.74259 | +0.80% |
| YV.RR52 | 3 / 5 | 0.68531 | 1.08% | stable | 0.67873 | +0.97% |
| YL.C05W | 3 / 5 | 0.57455 | 4.35% | stable | 0.58466 | -1.73% |
| YL.C09W | 3 / 5 | 0.58433 | 1.95% | stable | 0.57530 | +1.57% |
| 7D.J45A | 3 / 4 | 0.79432 | 4.84% | stable | 0.82866 | -4.14% |
| 7D.J54A | 3 / 4 | 1.01984 | 2.73% | stable | 1.07962 | -5.54% |
| 7D.J67A | 3 / 4 | 0.97301 | 5.08% | moderate | unavailable | unavailable |

Five qualifying stations are outside the YL design set. This supports temporal
and cross-band consistency for that subset, conditional on the existing
vertical response, water assumptions and historical RR52 Earth model. It does
not validate those shared assumptions or establish instrument accuracy.

The two YV medians are within 1% of their fixed Rayleigh references. Their
published Aminian values are 0.740 and 0.670, respectively; those publications
share the Rayleigh method family. The 7D offsets remain 4-6% and should be
investigated, not tuned away. H1-H2 envelopes and weekly spread are descriptive,
not confidence intervals or evidence of instrument failure.

## Decision for continuation

The subsequent [gain-band phase/model audit](microseism_gain_band_audit.md) is
complete. It retains all frozen estimates, adds RR40 and daily phase flags,
and supplies a read-only 89-row table staging export. Tilt and absolute response
are still unresolved; the bounded P/H1/H2/Z pilot is the next task.

Keep the method experimental and the frozen rule unchanged. Preserve all
earlier best-band results as legacy diagnostics. Do not rerun these completed
weeks, refit Rayleigh, apply gains, or populate calibration scores.

Next, use the exact passing-week paths in the JSON to check phase in the actual
gain bands and assess coherent tilt/vertical-response effects on RR52/RR40 and
the 7D controls. A bounded P/H/Z product would be needed to assess tilt; the
current saved P/Z pairs alone cannot exclude it. Investigate the 7D model/band
offset and C07W's time states separately. Report a second rule explicitly if
one is later proposed; do not retrospectively tune v1.

[All-station frozen-rule table](microseism_frequency_table.md) ·
[Legacy week table](microseism_weeks_table.md) ·
[Integration decision](reviews/microseism_integration.md).
