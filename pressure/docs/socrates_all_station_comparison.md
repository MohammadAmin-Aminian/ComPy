# All-station Socrates trial comparison

Completed user run: `all_20261013T195356Z_00cf3343`, code `e79f773`.
Scope: the saved no-tilt event windows of 89 distinct YV/YL/7D/XO stations.
Nominal replay reproduces the previous results. Alternative gains are provisional;
no correction was applied to original waveforms, responses, masters or scores.
All 2,432 protected source hashes were unchanged at completion and at comparison.

[Full 89-station table](socrates_all_station_comparison.csv) ·
[Provenance](socrates_all_station_comparison.json).

## Coverage

| Network | Previous station estimates | Pole trial / state-separated estimates | Ordinary unresolved | Reference exclusions |
|---|---:|---:|---:|---:|
| YV | 16 | 31 | 10 | 2 |
| YL | 22 | 22* | 3 | 0 |
| 7D | 4 | 4 | 11 | 0 |
| XO | 4 | 4 | 2 | 0 |
| Total | 46 | 61* | 26 | 2 |

*C07W uses three normal-state candidates instead of its mixed-state station
median. This is a diagnostic selection, not an accepted hardware correction.
Counts describe numerical QC outcomes, not final instrument acceptance.

## YV: pole reflection and fitting objective are distinct

Pole reflection gives 15 previously unresolved stations at least three numerical
QC candidates. Candidate events increase from 49 to 133 in the saved windows.
No previously estimated YV station loses its estimate. Gain changes at the
original 16 stations range from −2.82% to +1.60% under pole-only refitting.

| Newly estimated station | Pole-only pressure multiplier g | Candidates |
|---|---:|---:|
| RR01 | 1.041380 | 4 |
| RR05 | 1.120144 | 4 |
| RR06 | 1.139679 | 6 |
| RR10 | 1.083792 | 6 |
| RR16 | 1.078309 | 8 |
| RR19 | 1.129611 | 4 |
| RR20 | 1.341531 | 7 |
| RR21 | 1.114118 | 6 |
| RR26 | 1.088703 | 4 |
| RR30 | 1.185682 | 8 |
| RR43 | 1.140108 | 4 |
| RR44 | 1.307477 | 4 |
| RR46 | 1.015093 | 7 |
| RR48 | 1.386710 | 5 |
| RR55 | 1.592099 | 3 |

Remaining YV unresolved: RR12/RR17/RR18/RR22/RR31/RR42/RR45/RR47/RR53/RR56.
RR45 has two candidates, RR53 one; the others have zero. RR41/RR57 retain their
reference-band exclusions. More numerical QC passes do not verify the nominal
sensor model, pressure sign, vertical absolute response or mode isolation.

RR55's three newly accepted trial events were previously phase-review cases,
with nominal median phases about 29–31 degrees and reflected medians about
14–15 degrees. These are different from the four original boundary events.
The pressure response may represent the wrong hardware type; resolving that
identity remains necessary even though the trial now passes the existing mask.
RR22 stays unresolved; its coherent-event phase worsens under reflection.

The eight matched INSU pressure multipliers in
[Aminian et al. (2025), Table 1](https://www.ipgp.fr/~stutz/2025-aminian-crawford-stutzmann-al-GJI.pdf)
have median signed discrepancy **+3.20% before and +3.50% with weighted fitting**.
Pole reflection does not affect these INSU responses. Weighting improves four
individual comparisons and worsens four. This run uses the no-tilt baseline,
not the original tilt baseline used by the earlier Socrates V11 calculation.
Glitch processing, event sets and model choices still differ from the publication.

## YL: state separation helps; weighting increases the systematic discrepancy

Use **1/g** to compare with
[Zha & Webb (2016), Table 1](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012688).
All 21 estimated peer stations excluding C07W move farther from the published
coefficients under weighted fitting. Median discrepancy changes from **−5.43%
to −6.72%**; the range changes from −7.04…−3.37% to −7.86…−4.99%.
Pole reflection has no effect on YL.

C07W's mixed-state reciprocal was 3.71649, +86.76% from the published 1.99.
The three normal-state candidates give 1.90135 (−4.45%) with the original fit,
and 1.85971 (−6.55%) with weighted fitting. The smaller discrepancy supports
state-separated review; it does not establish the hardware mechanism or an
accepted correction. Do not multiply the recorded vertical samples by two.

A14W and B01W retain zero candidates; C10W retains one. B01W's coherent legacy
boundary events still have no phase-valid bins. Its high-pass fits remain
diagnostic and do not supply an accepted scalar calibration.

## 7D and XO: numerical sensitivity without a matched external benchmark

The same four station estimates remain in each network. Pole reflection has no
effect. Weighted fitting changes g by −0.43…+3.57% in 7D and −1.37…+0.98% in XO.
No verified original 7D factor table or directly comparable scalar XO coefficients
are available, so these changes cannot be ranked by published agreement.
LT02/WD48 remain unresolved. Their timing/coverage diagnoses are preserved.

## Integration recommendation

Keep the original fitting objective as the production default. The pressure-pole
alternative is promising as an explicit experiment; retain hardware/response and
held-out phase/arrival review before treating its 15 new station gains as valid
calibrations. Keep C07W state and B01W response-shape diagnostics separate from
automatic corrections. The engineering fixes are already merged; these scientific
options were subsequently integrated into main by the user's explicit request.
This integration retains their opt-in trial scope and the original production
fitter; it does not accept or apply the experimental gains.
