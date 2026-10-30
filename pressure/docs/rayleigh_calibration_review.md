# Teleseismic pressure calibration: source review and RR52 comparison

The next calibration pilot follows the earthquake method in Aminian's thesis
and manuscript. It is separate from the continuous-week pressure/displacement
measurement. A measured pressure/acceleration amplitude ratio is compared with
a water-column prediction to estimate a **candidate pressure multiplier**.
No response, waveform, master field or scientific score has been changed.

## Current update: provisional sign and actual TiSKitPy default

The user recalls that pressure and vertical motion had opposite signs and
authorized proceeding with the most probable convention. **RR52 therefore
defaults to pressure sign −1**, recorded as `assumed_needs_confirmation` in each
new JSON result. This is a working assumption supported by the old plotting
code and the observed approximately 180° phase difference, not verified wiring
or response metadata. A single sign is used for all RR52 events; it is not
chosen independently to make each event pass. Original responses and samples
are preserved. `--pressure-sign 1` makes an explicit comparison; other stations
require an explicit sign.

The default tilt method now invokes **TiSKitPy CleanRotator 2.3.1 itself**.
The earlier local eigenvector method remains available only by explicit choice,
with `needs_to_be_tested` status. [Method and short references](tilt_methods.md).

The actual-library eight-mainshock run gives:

| Event or group | Pressure amplitude multiplier | Current outcome |
|---|---:|---|
| Alaska, `ak0138esnzr` | 0.67648 | Candidate under assumed sign; visual window review pending |
| Santa Cruz, `usc000f1s0` | 0.66353 | Amplitude review only: predicted packet overlap |
| Papua, `usb000g7x7` | 0.66440 | Candidate under assumed sign; visual window review pending |
| Three other mainshocks | — | Insufficient band-mean coherence |
| Two September mainshocks | — | Quiet-training grids differ by 0.0001 s |

There are **two accepted pilot candidates**, below the chosen minimum of three
independent events, so station status remains `unresolved`. Their two-event
median is 0.67044; this is not a validated station calibration. The three
amplitude-review cases, including the overlapping event, have median 0.66440.
Signed pressure multipliers would be −g under the assumed convention; the
positive amplitude factor and reciprocal sensitivity factor are saved separately.
No factor has been applied. Tilt angles are about 0.33–0.50°.

Wide arrival-review plots now show predicted P/S, the group-arrival range,
the chosen window and other predicted packets. The phase plot shows both the
metadata sign and provisional reversal. [How to review windows and timing,
rerun one event, and open plots](rayleigh_window_review.md).

[Current validation](rr52_tiskitpy_validation_summary.json) ·
[Current comparison plot](figures/rr52/rayleigh_tiskitpy_comparison.png) ·
[Arrival overview](figures/rr52/rayleigh_tiskitpy_arrival_review.png).

The initial results below are retained as a historical comparison; their
local-method and unresolved-sign findings are superseded by this update.

### Subsequent Figure 5 comparison and magnitude fallback

The runner now adds the two-panel paper Figure 5 / thesis Figure 3.13 layout,
including its axes, labels, colors and aggregate gain objective. It also searches
Mw>7 → Mw≥6.5 → Mw≥6 while fewer than three QC candidates have been found,
with a bounded per-level attempt budget and unchanged scientific checks.
At the existing 0.02–0.06 Hz default, 24 attempts still give two candidates.
A separate run using the original figure's **0.03–0.07 Hz** band gives three
candidates, event-gain median **0.678998** and gain of the median ratio **0.674374**.
This remains a candidate comparison, without a response/master correction.
[Exact results, plot matching, package explanation and run commands](rayleigh_paper_figure.md).

## Material actually reviewed

The relevant calibration, tilt, model, figure and uncertainty sections were
read; this is not a claim to have audited every chapter or publication.

| Source | Relevant evidence |
|---|---|
| Aminian thesis, *Quantitative Analysis of Oceanic Crust Shear Velocity Structure Using Seafloor Compliance Measurement Data* | Section 3.1.6, printed pp. 53–60; Figs. 3.11–3.13; Table 3.1; equations 3.7–3.12; preceding tilt method and station/sensor description. Calibration equation and event/table pages were visually inspected. |
| Local *Complete manuscript file.pdf* | Section 4.1.5, Fig. 5 and Table 1; teleseismic Mw>7 selection, pressure/acceleration fitting, grid search and RR52 multiplier 0.67. It is a manuscript, not a verified byte-identical copy of the published article. |
| Aminian tilt milestone note | Variance minimization, transfer-function correction and coherence minimization. |
| BRUIT-FM WP2 proposal supplied locally | Separate PSD/cross-spectra measurements, calibration checks and later acceptance/classification. |
| [Zha & Webb (2016)](https://doi.org/10.1002/2015JB012688) | Section 2.3, Table 1 and Appendix A: teleseismic pressure/acceleration reference and frequency-independent amplitude correction assumption. |
| [Doran et al. (2019)](https://doi.org/10.1029/2019EA000783) | Independent in situ calibration context; event-centered Rayleigh packets and limits of scalar gain versus capillary time constant. |
| [Ruan et al. (2014)](https://doi.org/10.1002/2014JB011162) | Distinguishes microseism and earthquake transfer measurements; tilt and wavefield/model effects prevent treating all amplitudes as calibration offsets. |
| [Wayne Crawford's BRUIT-FM material](https://www.bruit-fm.org/crawford_noise_reduction.html), [proposal](https://www.bruit-fm.org/Documents/aapg-anr-2021-full_proposal_BRUITFM.pdf), and [TiSKitPy](https://github.com/WayneCrawford/tiskitpy) | Response validation, signal/noise separation, `CleanRotator`, event avoidance and `DataCleaner`. Current rotation implementation was inspected directly. |
| [disba](https://github.com/keurfonluu/disba) | Distinct fundamental-mode `PhaseDispersion` and `GroupDispersion`, period ordering and model units. |

Local originals remain outside this repository. Published article identifier:
Aminian et al. (2025), DOI `10.1093/gji/ggaf253`. Publisher access attempts failed;
the available manuscript and current public ComPy documentation were used, with
that version distinction retained. No private Wayne correspondence was inspected.

## Earlier code: what can be reused

| Code / inspected revision | Finding and use |
|---|---|
| `DPG_Calibration`, `7719622` | Original disba model, amplitude-ratio method and direct pressure-multiplier grid-search objective. Old code has reversed low/high-pass cutoffs, picks the largest signed vertical sample, assumes a particular trace length, omits a final event, and uses fixed frequency slices/reversals. The `velocity_model` argument is not passed through in its theoretical-ratio function. These findings belong to this older version. |
| Current public `ComPy`, `8b1431816f7720ac15855d38b17b3e91f2610fa2` | Several older filtering, trace-length and frequency-grid problems have already been repaired. It picks the largest absolute sample within an event span; that remains a candidate locator, not proof of a fundamental Rayleigh arrival. Its full cleaning chain should not be copied into earthquake windows without checking preservation of genuine horizontal/vertical Rayleigh motion. |
| `PressureCorrection`, `08d553c` | SAC-based phase-window and pressure-gain workflow, with manual P/S/R ranges and time-domain correlations. Inspected as a methodological precursor; not imported or executed because its file/SAC operations are unsuitable for this archive. |
| `SeafloorCal` (`dpgPy`), `520e3d9` | Strong explicit units, polarity/phase checks, physical model, group-velocity candidate windows and event-level uncertainty. Its historical fixes/validation notes were read before testing. Focused model/guard tests: **47 passed**. Numerical water-model comparison agrees with this pilot. This is not a new full audit or field certification of SeafloorCal. |
| `pressure-gauge-calibration-project`, `29a2db2` | Small acoustic-reference and spectral-ratio prototype. Its **10 scientific tests pass** in the BRUIT-FM environment. Useful normalization, grid, polarity and reciprocal-gain safeguards; no claimed previous field validation. |
| TiSKitPy source `3bb51b3eb82fa51f2e0888a063e484aa1c9ef5c3` | `CleanRotator` estimates a minimum-variance vertical axis in a low-frequency band while avoiding earthquakes. The initial pilot used a local covariance/eigenvector method. The current default runs the installed `CleanRotator` 2.3.1 class; this is not the full `DataCleaner` pipeline. |

The two recent packages were useful references rather than untested replacements
for the current archive reader. Their repositories/environments were not edited
or published. Existing dependency installations were preserved.

## Important differences in the historical descriptions

The thesis describes a 20-minute window after a PhaseNet P pick; the manuscript
describes amplitude analysis after P/S arrivals; the older code selects a maximum
vertical sample. These are different selection procedures. A P pick identifies
a body-wave arrival, not the dispersive Rayleigh packet. The pilot therefore
uses group velocity to propose a Rayleigh window and records a separate TauP P
prediction. It still requires visual/event-isolation review.

Table 3.1's frequency column is explicitly the **compliance** band, not the
calibration fit band. It must not be copied into the Rayleigh gain fit. The old
DPG default uses 0.02–0.06 Hz; the local paper plotting script supplies 0.03–0.07 Hz.
The current pilot uses the former for a documented first comparison. Published
studies' bands and coherence conventions are not universal acceptance criteria.

The old plotting code negates pressure in a normalized waveform display, while
the gain estimate uses unsigned spectral amplitudes. Thus visual agreement and
an amplitude fit do not establish pressure polarity. In the current raw-count
pilot, pressure and upward acceleration are nearly opposite in phase. This may
be a DPG/coordinate/response convention issue; it does not identify a defective
instrument or justify silently reversing a channel.

The table uncertainties are approximately five percent of their gain factors.
Without the original event-level uncertainty record, this review cannot certify
those values as empirical standard deviations or confidence intervals. Agreement
with the historical multiplier does not validate every previous preprocessing step.

## Physical reference and gain convention

For positive compressional pressure and upward acceleration, the adopted
homogeneous-water reference is

```text
r = 2πf * sqrt(1/α² - 1/c_phase(f)²)
T(f) = P/a = ρH * tan(rH)/(rH)
```

Its zero-frequency limit is ρH. The removable limit is evaluated stably; the
imaginary-r continuation uses tanh. A conservative band bound avoids the first
acoustic resonance. This pressure/acceleration reference has SI units kg/m².
It is not the inverse of the previously saved displacement/pressure estimator:
transfer direction, acceleration and noise bias must be handled explicitly.

The main measurement is `R = sqrt(Spp/Saa)`. Complex cross-power and squared
coherence are retained, with `H_PA = mean(P*conj(A))/mean(A*conj(A))` as an
additional estimator. The least-squares objective is `sum((g*R-T)²)`; the positive
grid 0.1–3.0 with step 0.001 is checked against `g = dot(R,T)/dot(R,R)`. A boundary
solution remains unresolved. This preserves the original direct-fit objective;
fitting `G*T` to R and reporting 1/G is a different loss under scatter.

`g` multiplies nominal pressure; `1/g` multiplies nominal sensitivity in counts/Pa.
For example, an amplitude-only g≈0.664 would imply pressure amplitudes about one
third lower and a counts/Pa sensitivity about 1.5 times larger. **Neither is
applied here.** A gain-only fit does not calibrate the low-frequency leak time
constant or justify extending a correction across the entire spectrum.

## What the new runner does

`rayleigh_calibration.py` is a **standalone station pilot**, not a network/master
calibration stage. It:

1. (Current default `--event-sampling all` evaluates every eligible event at all
   three levels; the target-and-budget behavior below is `bounded` mode.)
   Starts with Mw>7, source depth<150 km events from an explicitly supplied cached
   USGS catalogue, with a 24-hour mainshock separation proposal. If there are too
   few QC candidates, it expands to Mw≥6.5 and then Mw≥6, without repeating
   attempted events or relaxing QC. Counts and budgets are saved.
2. Centers a 30-minute candidate window on the predicted 40 mHz packet using
   **group** velocity. **Phase** velocity is used only in the water model.
3. Finds a preceding two-hour catalogue-screened quiet interval; estimates a
   minimum-variance near-vertical axis with TiSKitPy (constructor band
   0.001–0.01 Hz; internal optimizer 0.001–0.005 Hz). This does not guarantee
   noise-only data or a perfectly stable tilt axis. See [details](tilt_methods.md).
4. Reads synchronized raw pressure/BH1/BH2/BHZ samples without interpolation or
   implicit event resampling; compares exact response epochs/snapshots. All seismic
   coefficients are converted into the same physical velocity units before
   projection; acceleration is obtained with iω. Direct evalresp ACC conversion
   was checked independently to exclude a conversion-induced phase reversal.
5. Uses nine nonoverlapping 200-second periodic-Hann segments per event, removes
   segment means, corrects full complex responses once and averages auto/complex
   cross-power before deriving ratios/coherence. Frequency spacing is 0.005 Hz.
6. Requires mean/median squared coherence ≥0.95 and at least three coherent bins.
   Positive-convention candidates also require cross-phase within ±20° and no
   other predicted Mw≥6 packet overlapping the window. These numerical safeguards
   are explicit pilot choices, not an instrument pass/fail specification.
7. Applies the explicitly recorded provisional pressure sign to diagnostic
   cross-power and waveform copies, retaining nominal phase for comparison.
   Saves amplitude-only fits with a review flag when phase/arrival QC remains
   unresolved. No source response, score or master update is performed.

The RR52 hardware axes are treated as orthogonal for the quiet projection;
both horizontal inventory azimuths are zero, so geographic geometry is not
verified. The small projection is conditional on that hardware assumption.
Pressure is not used to predict/subtract the target Rayleigh signal from Z.

The dispersion table was generated with **disba 0.7.0** in the already-existing
SeafloorCal Python 3.12 environment. Its historical CRUST2/PREM-like layer block
was extracted without executing the old script; the last layer is declared a
half-space. It is a reproduction reference, not a verified local Earth model or
a verbatim PREM profile. Model/dispersion CSVs and provenance are shipped under
`src/bruitfm/reference_data/rr52_historical_rayleigh_*`. The BRUIT-FM environment
does not need disba to use the prepared table. The later TiSKitPy installation
added only the optional TiSKitPy/PyYAML/xarray packages; core scientific versions
were preserved.

## Initial RR52 field comparison (local method, metadata sign)

Eight selected mainshocks were examined. Three had adequate mean/median
coherence for amplitude review; three failed the band-mean coherence criterion;
two had unequal sample grids and remained unresolved. Aftershock candidates
within the 24-hour selection separation were separately recorded/excluded.

| Event | Median squared coherence | Amplitude-only pressure multiplier | Status |
|---|---:|---:|---|
| 2013-01-05 Alaska (`ak0138esnzr`) | 0.9860 | 0.67649 | Polarity review |
| 2013-02-06 Santa Cruz (`usc000f1s0`) | 0.9929 | 0.66351 | Polarity review; another predicted arrival overlaps |
| 2013-04-16 Papua (`usb000g7x7`) | 0.9966 | 0.66449 | Polarity review |

Median amplitude-only multiplier: **0.66449**, compared with the historical
**0.67**. The across-three-event sample standard deviation is **0.00722**;
this is scatter, not a confidence interval or total calibration uncertainty.
For Santa Cruz, the predicted first P is approximately 01:25:36 UTC, while the
candidate Rayleigh window is 01:41:11–02:11:11 UTC. A fixed P+20-minute endpoint
would capture only the beginning of this proposed window. The selected waveforms
and phase plot were visually reviewed; full isolated-mode validation remains pending.

The quiet projections are about 0.44–0.52°. Removing that small tilt changes
these events' amplitude multipliers by less than 0.001. Quiet low-frequency power
is reduced, but this is not evidence that every target-band contaminant is removed.
Using H1 instead of the PSD ratio changes the reported amplitude fit slightly;
both are saved so noise-estimator sensitivity remains inspectable.

All three amplitude-review cases fail the adopted near-zero phase convention.
Consequently **no signed calibration candidate passes and station status remains
`unresolved`**. The amplitude agreement supports a partial reproduction of the
historical result, not confirmation of a fully calibrated station.

A sensitivity calculation illustrates reference dependence: perturbing assumed
depth by ±5% moves the median amplitude multiplier from about **0.626 to 0.705**;
using only 0.03–0.06 Hz gives **0.672**. These perturbations are illustrative,
not estimates of actual bathymetric uncertainty. This dependence is larger than
the three-event scatter and must be represented before quoting total uncertainty.
[Full sensitivity table](rr52_rayleigh_model_sensitivity.csv).

[Validation summary](rr52_rayleigh_validation_summary.json) ·
[Amplitude/coherence/phase plot](figures/rr52/rayleigh_calibration_comparison.png) ·
[Waveform window](figures/rr52/rayleigh_calibration_waveform.png).

## Run, outputs and preservation

On this machine the explicit catalogue snapshot is already saved locally:

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python -m pip install -e '.[rayleigh]'
conda run --no-capture-output -n bruitfm python -m bruitfm.rayleigh_calibration \
  --catalog results/RAYLEIGH_REVIEW/catalog_usgs_m5p5_2012_2013.json \
  --dispersion-file src/bruitfm/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --network YV --station RR52 --max-events 8
```

Use `--event-id usc000f1s0 --max-events 1` to inspect only the thesis's example.
`--show` opens the comparison plot. Each run creates a new pressure-epoch folder
`rayleigh_calibration/run_<unique-id>/`, saving event spectra NPZ/CSV, grid-loss
CSV, waveform, wide arrival-review and comparison PNGs, a window-review template
and a JSON audit. The final `paper_figure5_comparison.png` combines accepted
events in the original figure's layout. Existing runs are preserved.
No standalone mode writes masters or responses. Other stations require their own
validated reference/parameters; this pilot does not generalize horizontal pairing.

To fetch a new catalogue into a **new filename**, use the USGS FDSN event API:
`https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson&starttime=2012-12-29&endtime=2013-11-25&minmagnitude=5.5&orderby=time-asc`.
The ≥5.5 snapshot has 426 events for quiet-interval screening; the calibration
selection starts at Mw>7 and can fall back to Mw≥6.5/6. Catalogue/input hashes are saved. Full raw products
and that snapshot remain local, excluded from Git.

The full BRUIT-FM suite passes **87 tests**. Checks independently integrate
the acoustic ODE, verify grid/analytical and reciprocal gain conventions, group
versus phase velocity, exact synthetic rotation, raw-response gain recovery,
phase reversal, actual-library projection, provisional-sign gain recovery,
manual event exclusion, timing diagnostics and master/raw preservation. Source waveform size/mtime
checks and inventory/snapshot hashes pass; production CSV/Excel hashes are unchanged.
Magnitude fallback, stopping at the QC-candidate target, cross-level clustering,
deduplication and explicit-event/no-fallback controls are also checked.

## Decision before applying any correction

Proceed with the documented −1 pressure assumption for the RR52 pilot. Seek
hardware/response confirmation before a final instrument correction, inspect
isolated Rayleigh arrivals, confirm the reference/depth and
quantify event/window/model uncertainty. Keep all scores unchanged until there
is a defensible acceptance rule. Do not automatically apply 0.664 or the old 0.67,
refit a sensor time constant, or launch a network calibration batch from this pilot.
