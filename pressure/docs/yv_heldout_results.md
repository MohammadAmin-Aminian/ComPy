# Completed frozen-gain held-out validation

User run: `run_20261013T213126Z_87e24fe6`, code `026ebf1`.
All 50 selected unused earthquakes were attempted at RR01/RR06/RR48/RR52.
The saved pole-only gains were frozen: no validation refit or source correction.
[CSV](yv_heldout_results.csv) · [Provenance and counts](yv_heldout_results.json).

| Station | Frozen g | Isolated usable events | Median prediction error | Conclusion |
|---|---:|---:|---:|---|
| RR01 | 1.041380 | 0 | — | Insufficient independent isolated evidence |
| RR06 | 1.139679 | 1 | −1.041% | Encouraging single event; insufficient |
| RR48 | 1.386710 | 0 | — | Epoch coverage and coherence limit the test |
| RR52 | 0.678728 | 14 | +0.664% | Good conditional held-out consistency |

Error means `g_frozen * observed_P/a / predicted_P/a - 1`, not a difference from
a published coefficient. RR52's event-error 16th–84th percentiles are −2.35%
and +2.07%; its full isolated-event range is −5.96%…+4.91%. These are event
spread under the shared response/model/estimator, not absolute calibration
confidence. Fundamental-mode isolation and hardware calibration remain pending.

## Why the other stations remain unvalidated

- **RR01:** 39 coherence failures, ten input/window refusals and one body-arrival
  review case. That review case has error −2.058%, but is not an isolated event.
  Input refusals include five epoch limits, two within-channel native-grid offsets
  and three windows beginning before origin. No safeguard was relaxed.
- **RR06:** 37 coherence failures, ten input/window refusals, two body-arrival
  review cases and one isolated event. The three QC-usable events have errors
  −0.526%, −1.041% and +1.317%, but only the −1.041% event is in the primary
  isolated summary. This is insufficient for independent calibration acceptance.
- **RR48:** its selected pressure and vertical epochs end on
  `2013-06-10T23:28:43Z`. Thirty-eight selected earthquakes are outside that epoch;
  three predicted windows fail the origin guard, and nine fail coherence.
  Absence of usable validation does not prove the gain is wrong.
- **RR52 control:** 14 isolated events, eight body-arrival review cases,
  24 coherence failures and four input/window refusals. The control's held-out
  consistency supports this calculation's measurement procedure, while sharing
  the reference and nominal-response assumptions of the training method.

Nominal response QC yields zero usable events at RR01/RR06; pole reflection yields
one/three including body-arrival review. This supports the phase-convention
experiment, but is not enough independent isolated data to certify those gains.

## Preservation and next step

The completion report and post-run check confirm unchanged **71 source/metadata/
code hashes** and **2,208 indexed raw-file size/mtime records**. Raw byte hashes
are not claimed. Waveforms, XMLs, snapshots, masters, prior results and all scores
remain intact. The current suite has 173 passing tests; no source code changed
while interpreting this run.

Keep the new YV gains provisional. Review RR01/RR06's flagged arrivals in the
waveforms and select a longer-epoch station to replace RR48 in further unused-
event validation. Do not widen epochs or lower QC to manufacture usable events.
The original fitter remains the default. This study does not complete absolute
calibration validation for the 15 newly estimated YV stations.
