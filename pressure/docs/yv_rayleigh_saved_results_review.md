# YV Rayleigh: initial review of saved results

The user completed the 43-station computation. This review examines those saved
products without repeating the network run or applying gains. The result is an
arithmetic audit and an initial scientific review, not final calibration acceptance.

## What was checked

All 43 station reports and recorded event reasons were read. For all **481 saved
event spectral products**, pressure/acceleration ratios, H1 transfers, squared
coherence, phase and accepted-bin masks were rederived from saved auto/cross power
arrays, without calling the production functions. The water-column expression
was separately evaluated from the supplied phase table and declared water/depth
parameters. The **48 candidate-event gains**, reciprocal sensitivities, signed
gains, all **16 station medians** and all **16 median-curve figure gains** were
recomputed with NumPy least squares. No discrepancies were found beyond numerical
roundoff (maximum event gain difference approximately 1.1e-15).

These are downstream arithmetic checks. They reuse saved spectra, the adopted
water equation and historical velocity table. They do not independently test raw
FFT normalization, complex response removal, the Earth/ocean reference, or absolute
sensor sensitivity. Those separate pipeline and field checks remain pending.

All 16 combined Figure 5 plots were visually inspected. Representative broad
arrival and waveform plots were inspected for RR03/usb000ief9 and
RR07/RR40/RR52/usc000f40j, selecting the candidate with the largest recorded
relative RMS residual at each station. This is not a full inspection of every
candidate event or quiet-training interval. No review flag or score was cleared.

## Candidate estimates and sensitivity

The table gives the **median event pressure multiplier**, distinct from the
median-curve gain shown in the plot title. Event CV is sample scatter divided by
the median, not a confidence interval. Relative RMS is the largest event fit
residual norm divided by the theory norm in its accepted bins; it is not an
adopted pass/fail threshold. Depth perturbations hold spectra and all other
parameters fixed and are illustrative, not measured depth uncertainty.

| Station | Median gain | Event CV (%) | Largest event RMS (%) | Gain at depth −5% | Gain at depth +5% |
|---|---:|---:|---:|---:|---:|
| RR03 | 1.1224 | 2.49 | 8.0 | 0.9851 | 1.3038 |
| RR07 | 1.0822 | 3.92 | 14.2 | 0.9612 | 1.2326 |
| RR08 | 1.1579 | 2.17 | 6.9 | 1.0392 | 1.3059 |
| RR09 | 1.2059 | 3.16 | 3.0 | 1.1244 | 1.2923 |
| RR11 | 1.0637 | 2.14 | 4.3 | 0.9725 | 1.1655 |
| RR13 | 1.2087 | 1.28 | 9.7 | 1.0836 | 1.3572 |
| RR14 | 1.1035 | 4.24 | 5.8 | 1.0172 | 1.1979 |
| RR25 | 1.1648 | 1.23 | 8.4 | 1.0241 | 1.3749 |
| RR28 | 0.8195 | 3.45 | 7.9 | 0.7273 | 0.9358 |
| RR29 | 0.7299 | 0.84 | 11.6 | 0.6512 | 0.8292 |
| RR34 | 0.5848 | 2.36 | 10.3 | 0.5274 | 0.6530 |
| RR36 | 0.6513 | 1.20 | 8.5 | 0.6027 | 0.7042 |
| RR38 | 0.6364 | 2.16 | 8.5 | 0.5632 | 0.7134 |
| RR40 | 0.7419 | 4.95 | 18.2 | 0.6602 | 0.8854 |
| RR50 | 0.7165 | 0.52 | 11.8 | 0.6517 | 0.7912 |
| RR52 | 0.6790 | 0.69 | 6.5 | 0.6357 | 0.7244 |

Across these stations, a −5% depth perturbation moves gains about −6% to −12%;
a +5% perturbation moves them about +7% to +19%. Event scatter is only about
0.5–5%. The small scatter therefore cannot establish absolute calibration accuracy.
RR52 remains approximately 0.679; its illustrative depth range is 0.636–0.724.

RR07 and RR40 deserve priority for scalar-fit shape/window review: maximum event
RMS residuals are 14.2% and 18.2%. Other candidate curves also diverge from the
reference outside the red fit-band limits. No correction is justified there by
this scalar fit. The dashed “Corrected” curve represents a plotted multiplication,
not a gain applied to waveforms or responses.

The selected RR03/usb000ief9 window begins before the predicted S arrival, which
lies inside the selected interval. The strongest displayed packet is later.
The automatic guard for other earthquakes does not establish isolation from the
same earthquake's body-wave energy or modal mixtures. This event remains a
candidate requiring window/mode review; it has not been silently rejected or
replaced. The other sampled waveforms show similar pressure/vertical packet
structure, which supports comparison but does not identify a pure mode by itself.

## Why the other stations remain unresolved

All 25 unresolved stations have zero accepted candidate events in this run.
They have combinations of the following **event-level** reasons; these counts
include rejected attempts at candidate stations and do not count instruments:

| Recorded cause | Event attempts |
|---|---:|
| Insufficient band coherence | 303 |
| Phase review | 118 |
| Predicted other-earthquake packet overlap | 12 |
| TiSKitPy rotation-angle limit | 72 |
| Missing, gapped, nonfinite or flat window | 44 |
| No catalogue-screened quiet training | 25 |
| Invalid rotation projection | 9 |
| Sample-grid offset | 6 |
| Gain-grid boundary | 4 |

RR12/RR42 primarily lack usable windows, consistent with the earlier flat-data
and short-coverage investigations. RR45 has 12 rotation-angle rejections and one
quiet-training rejection; RR53 has 19 angle rejections plus phase/coherence and
quiet-training issues. These need geometry, metadata and training review, not
an automatic increase of the allowed rotation.

RR55 reaches the gain-search boundary in four attempts. The failed fits were
not saved as valid event gains, so this review cannot certify their extrapolated
solution or justify increasing the search range.

The machine label `polarity_review` also covers general phase disagreement.
It does not prove reversed wiring. For example, RR22/ak0138esnzr phases range
approximately −57° to −92°; RR55/usb000g7x7 ranges +30° to +53°, despite high
coherence in many bins. A 180° sign change does not restore either to zero phase.
RR56 has a different, frequency-dependent phase pattern. Actual pressure and
vertical responses, timing, instrument type and wavefield must be investigated;
do not flip signs per event or relax phase QC to manufacture candidates.

## RR41/RR57 reference-band flags

The preflight guard uses the **maximum phase velocity of the entire table** at
the requested upper frequency, giving a conservative resonance bound. At the
inventory depths, that guard's maximum safe frequencies are approximately
0.06690 Hz for RR41 (5430 m) and 0.06986 Hz for RR57 (5200 m).

Using the frequency-dependent table instead, the 0.07 Hz water-model argument
is 0.938 of the first pole for RR41 and 0.898 for RR57. The configured cutoff is
0.900. RR57 is consequently a borderline conservative preflight exclusion, not
an established unusable exact model or sensor failure. Depth/model uncertainty
is unresolved. Preserve both saved flags; consider a justified lower band, such
as 0.03–0.06 Hz, only as a separate comparison after reviewing the reference.
No band, guard or saved status was changed by this audit.

## Follow-up work

1. Inspect all candidate-event arrival/phase/quiet-training products and document
   packet choices, starting with RR03, RR07 and RR40. Test held-out events and
   justified alternate windows instead of stopping at three successes.
2. Independently estimate representative raw spectra/complex responses with a
   separate implementation. Verify actual hardware, vertical sensitivity, depth,
   FIR warnings and station/path velocity model.
3. Investigate coherent phase-disagreement stations (for example RR22/RR55/RR56)
   and rotation-limited RR45/RR53. Review a lower-band comparison for RR41/RR57.
4. Quantify systematic and sampling uncertainty and agree acceptance rules before
   applying any gain or assigning a score. Then complete microseism/local-P checks.

The masters were only read/hashed. Existing spectra, responses, scores and
workflow statuses were preserved. Source computation and scientific implementation
were not changed; no long batch was started.

[Detailed numerical audit](yv_rayleigh_saved_results_review.json) ·
[Completed run](yv_rayleigh_completed_summary.json) ·
[Assumptions/checklist](calibration_assumptions.md) ·
[View existing plots](yv_rayleigh_calibration.md#view-without-recomputing).
