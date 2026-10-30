# Completed YL run compared with Zha & Webb

All 25 available pressure epochs completed computational processing. **22** have
provisional station estimates; A14W/B01W/C10W remain unresolved. No gain was applied
and the calibration master was not updated. Source workflow:
`results/NETWORKS/YL_ALL_NO_TILT/_workflow/rayleigh.json`.

The comparison uses **1/g**, where BRUIT-FM fits `g * measured_P/a = predicted_P/a`.
Zha's Figure A2 describes scaling the theoretical ratio by the reported coefficient
to compare with the observations; reducing pressure would use its reciprocal.
This matches correction direction, not necessarily fit objectives or nominal
response versions. The paper does not supply its software implementation here.

Published values and event counts below are from
[Zha & Webb (2016), Table 1](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012688).

| Station | Our reciprocal 1/g | Published coefficient | Relative difference | Our / published events |
|---|---:|---:|---:|---:|
| A02W | 1.436 | 1.50 | −4.2% | 8 / 19 |
| A06W | 1.446 | 1.53 | −5.5% | 3 / 19 |
| A08W | 1.369 | 1.43 | −4.2% | 9 / 17 |
| A10W | 1.709 | 1.82 | −6.1% | 9 / 17 |
| A12W | 1.564 | 1.67 | −6.3% | 19 / 17 |
| A14W | Unresolved | 1.97 | — | 0 / 9 |
| B01W | Unresolved | 0.19 | — | 0 / 18 |
| B04W | 1.354 | 1.41 | −4.0% | 11 / 19 |
| B09W | 1.261 | 1.32 | −4.5% | 12 / 19 |
| B11W | 1.788 | 1.90 | −5.9% | 6 / 17 |
| B14W | 1.539 | 1.62 | −5.0% | 15 / 19 |
| C01W | 1.736 | 1.85 | −6.2% | 5 / 17 |
| C03W | 1.454 | 1.55 | −6.2% | 7 / 19 |
| C05W | 1.710 | 1.77 | −3.4% | 20 / 18 |
| C07W | **3.716** | **1.99** | **+86.8%** | 21 / 11 |
| C08W | 1.687 | 1.78 | −5.2% | 10 / 19 |
| C09W | 1.738 | 1.84 | −5.5% | 31 / 19 |
| C10W | Unresolved | 1.51 | — | 1 / 7 |
| C11W | 1.731 | 1.83 | −5.4% | 14 / 16 |
| C16W | 1.692 | 1.79 | −5.5% | 14 / 19 |
| C17W | 1.368 | 1.43 | −4.4% | 24 / 19 |
| F02W | 1.367 | 1.47 | −7.0% | 17 / 17 |
| N02W | 1.639 | 1.73 | −5.3% | 16 / 19 |
| N03W | 1.721 | 1.80 | −4.4% | 20 / 19 |
| S01W | 1.683 | 1.81 | −7.0% | 5 / 16 |

The 21 numerical comparisons excluding C07W are consistently 3.37–7.04% below
the published coefficients (median −5.43%). These are descriptive discrepancies,
not uncertainty intervals, significance tests or instrument acceptance rules.
The current 0.03–0.07 Hz borrowed model, event/window policy, phase/timing masks
and median estimator differ from the paper's 0.04–0.1 Hz regional-model procedure.
Our counts include accepted M>=6 events; their Table 1 reports M>=7 sampling.
Thus counts exceeding the paper's counts do not establish equivalent sampling.

## C07W needs targeted investigation

Its median pressure multiplier is 0.269071 (reciprocal 3.71649). The 21 saved
candidate event multipliers range from 0.259217 to 0.542620. The four already
accepted M>=7 events alone give a diagnostic reciprocal of 3.70847, so merely
restricting those accepted events to larger earthquakes does not remove this
discrepancy. This subset is not a newly accepted station fit or the paper's set.
Some individual events are near the published reciprocal direction; many differ
by roughly a factor of two. The dates are interleaved, not evidence of a simple
chronological gain step.

The supplied pressure sensitivity/first-stage gain match C08W/C09W
(−753.235 / −0.00179585), and vertical sensitivities are similar. This rules out
an obvious factor-of-two difference in those listed metadata values; it does not
validate actual hardware calibration, response shape, waveform identity or
Rayleigh isolation. Inspect matched events and response conventions before
attributing the discrepancy to hardware or modifying any sensitivity.

## Three unresolved stations

- **A14W:** no accepted candidates. 52 attempts fail coherence; timing, coverage
  and default windows beginning before origin also contribute.
- **B01W:** no accepted candidates. 35 coherence rejections, many timing/data
  rejections and four gain-search boundary failures. The published coefficient
  0.19 would imply `g = 1/0.19 = 5.26` under the reciprocal convention, outside
  the current pressure-multiplier search range 0.1–3.0. That is a software-policy
  limitation to investigate, not proof that this instrument cannot be calibrated.
- **C10W:** one accepted event, below the existing three-event minimum. 51
  coherence rejections, timing/data problems and overlap/phase review contribute.

Exact saved reasons and per-event values are retained in the numerical audit.
No threshold, timing policy, gain bound or model was changed for this comparison.

## Preservation and reproduction

`scripts/compare_yl_zha.py` reads saved comparison JSONs only, verifies source/
master hashes, and writes a new `run_*` folder under
`results/VALIDATION/YL_Zha_comparison/`. It does not read raw waveforms, recompute
spectra, refit gains or write scientific scores. The initial snapshot is
`run_ef52555843c8`; arithmetic/count checks and Ruff passed.

[CSV](yl_all_zha_comparison.csv) · [Numerical audit](yl_all_zha_comparison.json).
Next prioritize C07W's matched-event response/window checks, then timing and
published-model/band/event/estimator reconciliation. Preserve this snapshot.

> **Codex integration qualification:** the daily relative vertical states are
> supported, but a C07W hardware fault is not proven. Peer-event ratios require
> independent peer QC: C08W's 2010-05-31 event has median coherence 0.285,
> so near-unity pressure ratios cannot be assumed for every shared event ID.
> See [the integration review](reviews/codex_integration_review.md) for bounded raw checks and limits.

## Sentinel follow-up: independent evidence on C07W and B01W

These notes add evidence; the table and saved results above are unchanged.
Details and reproduction are in [the Sentinel audit](SENTINEL_AUDIT.md)
(SENT-001/004) and [its evidence files](sentinel/).

- **C07W:** the 21 event gains form two groups (17 near g≈0.27, 4 near g≈0.53),
  not a continuum. Against neighbor C08W on the same events, pressure amplitude
  agrees (ratio ≈1) in both groups, while C07W's **vertical (LHZ)** amplitude is
  ≈0.5× in the low-g group and ≈1× in the high-g group. Daily microseism-band
  ratios against C08W and C09W are bimodal over the deployment (251 days ≈0.44,
  68 days ≈0.85; C08W/C09W steady ≈0.86), with every high-g event inside a
  "normal" interval. Normal-state events give 1/g ≈1.90, consistent with the
  other stations' offset from 1.99. The current C07W station median reflects an
  intermittent vertical-channel amplitude state, not the pressure gauge.
- **B01W:** the four boundary events are highly coherent (median γ² 0.97–0.99)
  but have **no phase-valid bins**: P/a phase is +45° to +70° at every frequency,
  while B04W/A02W are 0±10° on the same events. B01W's R/T rises from 0.11 at
  0.03 Hz to 0.23 at 0.10 Hz; its 0.04–0.1 Hz average reproduces the published
  0.19. The binding limitation is a frequency-dependent pressure response, not
  the 0.1–3 search bound; widening the bound would not yield accepted events.
