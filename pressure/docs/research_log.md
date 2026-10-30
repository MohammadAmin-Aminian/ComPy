# BRUIT-FM: completed steps

The [README completed steps](../BRUITFM_README.md#completed-steps) lists forty completed
steps in order. Its [YV commands](../BRUITFM_README.md#run-every-completed-stage-for-the-yv-network)
show the current execution stages. This log keeps the detailed calculations
and validation evidence.

## Inventory and data structure

### 1. Read the StationXML inventories

Read inventories organized by network or station, including nested folders,
and collected the available station, channel and instrument information.

### 2. Keep channel epochs and identify duplicates

Kept a separate record for each channel's metadata validity period. Flagged
repeated SEED IDs and time periods so different instrument configurations are
not mixed or silently discarded.

### 3. Build the master metadata table

Created a master table with **16,784 channel-epoch records**, including SEED ID,
coordinates, estimated depth, sensor, sampling rate and response information.
Depth is estimated from inventory elevation and still needs bathymetric validation.

### 4. Export the master table to CSV and Excel

Saved the same metadata in both formats, preserving identifiers and complete
text fields. The master files provide the link between the inventories,
waveform processing and saved results.

### 5. Prepare the six quality-test fields

Added score and status columns for expected pressure levels, nearby stations,
predicted infragravity levels, and pressure–vertical comparisons using Rayleigh
waves, microseisms and local P waves. **The fields are ready; the tests have
not received acceptance scores.**

## RR52 processing and results

### 6. Read the four RR52 channels

Processed station **YV.RR52**, location **00**: `BDH` for pressure, `BH1` and `BH2`
for horizontal motion, and `BHZ` for vertical motion. Identified **329 daily
MiniSEED files per channel**, or **1,316 files in total**, covering the available
recording period.

Read the data in bounded time intervals using actual file-header times. Handled
daily-file overlaps and missing samples without interpolation or zero filling.
The original waveform files remain unchanged.

### 7. Apply instrument responses and calculate PPSDs

Used each channel's matching StationXML epoch and full instrument response.
ObsPy applies response correction once to the raw counts. Pressure spectra
use **Pa²/Hz**; seismic spectra use acceleration units **m²/s⁴/Hz**. Values
are stored in dB relative to these units.

Used **3-hour windows with 50% overlap**, giving a new spectrum approximately
every **90 minutes**. Applied one-octave smoothing with 1/8-octave steps. At the
62.5 Hz sampling rate, the saved frequency centers extend from approximately
**0.000954 Hz to 31.25 Hz**, meeting the 0.001 Hz computational target.
The instrument's usable low-frequency limit still needs verification.

### 8. Complete the full calculation and save the results

Ran the four channels in parallel using **four worker processes**. Obtained
**5,249 PSD windows per channel**, or **20,996 in total**, with **121 frequency
centers per spectrum**.

Saved reusable PPSD NPZ files, response snapshots and processing settings under
`results/TEST/YV/RR52/`. Added checkpoints so interrupted calculations can resume
without duplicating spectra. Updated the four matching master CSV/Excel rows
with result paths, units, counts, dates, parameters and processing status,
while preserving existing metadata and quality-test fields.

### 9. Produce three plots for every channel

Created **12 figures**, three for each channel:

- **PPSD distribution:** how often each spectral level occurs at each frequency.
- **Time spectrogram:** how spectral levels change through the recording period.
- **Individual spectra:** all 5,249 PSD curves drawn together for comparison.

Added the **New Low Noise Model (NLNM)** and **New High Noise Model (NHNM)** to
the seismic distribution and curve plots. Pressure plots retain pressure
units. [View the twelve figures](figures/rr52/README.md).

### 10. Check the outputs and assess temporal resolution

Reloaded all four NPZ files and verified finite spectra, unique window times and
agreement between saved data and plots. Widened the histogram range to
**−250 to 100 dB** to include the quiet seismic tail. Final processing records
contain no chunk errors or warnings. Checked preservation of the master data
and source file attributes. The original processing checks passed; the expanded
suite now contains **33 automated tests**.

The current time plots represent **3-hour averages**, even though columns are
spaced 90 minutes apart. Checked shorter-window alternatives: 1 hour at 0.003 Hz,
30 minutes at 0.01 Hz, and 15 minutes at 0.02 Hz satisfy the current FFT-support
check. These alternatives have **not been calculated on the full dataset**;
increasing overlap alone does not shorten the averaging window.

### 11. Compare saved pressure spectra with expected levels

Implemented the first pressure-level comparison using the existing BDH NPZ.
Calculated 5th/50th/95th percentile spectra, compared selected frequency bins
with the supplied Webb bounds and approximately digitized Brown hydroacoustic
curves, and saved a comparison plot, frequency CSV and JSON settings. Reused
all **5,249 saved pressure spectra**, without waveform reprocessing or another
instrument-response correction. Converted Brown's µPa²/Hz reference levels
to the Pa²/Hz units of the saved spectra.

Selected a nominal **0.00857–7.395 Hz** band from the saved instrument response
and excluded smoothing bins outside it. Nineteen bins contributed; sixteen median levels
were within the reference bounds and three above. Recorded the diagnostics
in the master CSV/Excel, with status `comparison_only` and no calibration score.
The reference values and usable band still need review. Webb's bounds were
transcribed from the supplied instructions; Brown's curves are approximate
digitizations. These comparisons have not established calibration acceptance.
Tests cover reference units, bin selection, response-band estimation, score
opt-in and preservation of existing master data. The expanded suite passes
the pressure-comparison checks.
[Comparison details and commands](pressure_comparison.md).

### 12. Make the pressure comparison visible on the PPSD

Added a pressure PPSD distribution with Webb limits in magenta and Brown
limits in orange. Added separate reference panels showing the median's dB
clearance from the bounds and the percentage of PSD windows within range.
The same command displays all three comparison figures and saves their paths
in the master CSV/Excel. The comparison values and calibration status are unchanged.

The three views are:

- [PPSD with limits](figures/rr52/BDH_ppsd_with_limits.png): occurrence probability
  with lower/upper bounds from both references.
- [Test details](figures/rr52/BDH_test_details.png): reference-specific median
  clearance and window percentages, including digitization review markers.
- [Percentile comparison](figures/rr52/BDH_pressure_comparison.png): observed
  5th/50th/95th percentile spectra, reference bounds and nominal response band.

From the project checkout:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.pressure_test
```

Verified the plots against saved comparison values and checked that existing
master values, scores and source NPZ were preserved. All **21 automated tests**
and GitHub checks passed. The source implementations are recorded in
[Step 11](https://github.com/MohammadAmin-Aminian/bruit-fm/commit/ffbd8233f0bc923088e842100aac1394fb283958)
and [Step 12](https://github.com/MohammadAmin-Aminian/bruit-fm/commit/cbcded709edfea17c5140d33c9aa2143f58582a2).

### 13. Compare pressure spectra with the predicted infragravity model

Implemented Godin's depth-dependent pressure model from Equation 6, solved
the gravity-wave dispersion relation, and matched the model's frequency-bin
averaging to the saved pressure spectra. Reused all 5,249 BDH spectra and
the 2,918 m inventory-derived depth estimate.

With the nominal sensor band and conservative 0.001–0.02 Hz model range,
only two full smoothing bins can be compared. Their median observed-minus-
predicted difference is +16.13 dB. Saved a percentile/model view, a PPSD model
overlay and a residual plot, with numerical results and links in master CSV/Excel.
Status is `comparison_only`; no calibration score or correction was inferred.

Verified model normalization and dispersion, matching bin smoothing, frequency
exclusions, empty-overlap handling, plot values and preservation of master data
and NPZs. The expanded suite passes **33 tests**. Depth, usable instrument band
and local applicability still need review. [Model details and command](infragravity_comparison.md).

### 14. Independently recheck the Godin calculation and clarify Figure 3

Verified the model with an independent scalar root solver and exact ObsPy
FFT-period masks. Recomputed the first pressure PSD from raw counts and
matched the stored spectrum exactly. Rechecked units, percentiles, residual
signs, frequency selection, source hashes and master-file preservation.

The two retained residuals are +15.83 and +16.43 dB; their summary remains
+16.13 dB. Revised Figure 3 to show the full diagnostic residual curve and
the two eligible bins separately, with frequency widths and window spectral
spread. Excluded values do not contribute to the summary, and the strongly
overlapping bins are not independent measurements. No calibration acceptance
or correction follows from this diagnostic. The expanded suite passes 33 tests.

### 15. Test full-length Hann spectra on a week near mid-recording

Calculated pressure spectra directly from raw BDH counts using three
overlapping 12.42-hour Hann windows per day. Used a middle-of-recording week
to avoid deployment startup data. All 21 windows were continuous; screening
the highest notch powers at 0.05–0.06 Hz retained 20 and rejected one.
Applied the full pressure response once and averaged retained linear powers.

Saved reusable spectra, screening decisions, model residuals and three plots.
Compared with 110 existing PPSD windows from the same week in matching
frequency bins. The two fully eligible matched bins are 1.79 and 2.32 dB
lower than the same-week PPSD mean. The primary raw FFT mean has a median
model difference of +17.79 dB across the eligible frequency band. This is a
different estimator and interval from the earlier diagnostic.

Added separate pilot diagnostics and file paths to master CSV/Excel,
preserving every existing field and leaving calibration scores blank.
Verified FFT normalization, screening, sample alignment, response correction,
linear averaging and file preservation; the suite passes 40 tests.
[Method, outputs and commands](godin_pilot.md).

### 16. Run the completed workflow across automatically discovered networks

Added one batch entry point for safe inventory synchronization, availability
checks, full-period PPSDs, all three base plot types, expected pressure limits,
the Godin model comparison and the screened full-window pilot. Network,
station, location and channel selection come from the waveform archive;
response epochs and sample rates are handled independently.

Added channel-level parallel PPSD jobs, bounded reads, saved-stage reuse,
cooperative stopping and isolated failure/retry handling. Only the parent
updates CSV/Excel, with backups and preservation of manual results and extra
sheets. Adding another network's waveforms and StationXML uses the same
commands. Original RR52 products remain available beside the new network
output directory. Calibration scores remain blank.

Verified the stages on real pressure and vertical-motion data at 50, 62.5
and 100 Hz in an isolated export: four channels produced 110 PPSD windows
each; three pressure pilots retained 20 of 21 windows each. Confirmed empty
Godin-bin handling without assigning a failed calibration. The full archive
is ready for processing; this verification did not compute all its recordings.
The expanded suite passes 47 tests.
[Input layout, commands and outputs](network_workflow.md).

### 17. Compare nearby pressure stations on common recording days

Added the second spectral diagnostic using saved response-corrected pressure
PPSDs. Screened geodesic distance, relative water depth, declared sensor labels
and overlapping qualifying UTC days. Proposed 50 km/5% preferred and
100 km/10% exploratory criteria remain configurable and independent of network
spacing. Require at least two distinct reference stations and 30 jointly
supported days per frequency. Single-pair and unsupported comparisons do not
receive a failed calibration or acceptance score.

Calculated daily medians with coverage checks, aligned only supported frequency
bins and summarized paired differences against an equally weighted station
median. Saved candidate audits, per-frequency and daily results, source hashes
and three plots: spectra/residuals, changes over time and neighbor selection.
Added the stage to automatic network discovery, restart and saved-plot viewing.
CSV/Excel receive diagnostic values and paths while preserving existing scores,
manual fields and extra sheets. Changed reference products invalidate saved
comparisons. The workflow checks for unfinished PPSD processing before writing.

Verified a known 6.0206 dB synthetic offset, geographic/depth screening,
frequency alignment, duplicate-station exclusion and joint coverage. Isolated
copies of five completed YV pressure products gave four exploratory comparisons
and one without suitable neighbors in the selected reference pool. RR44 used
three peers over 153 days, with a -1.99 dB median difference; this is a spectral
diagnostic, not a gain correction. Source products and scores were preserved.
The expanded suite passes 62 tests.
[Criteria, outputs and commands](nearby_comparison.md).

## Step 18 — paired pressure and vertical displacement

Implemented bounded synchronized raw reads, full complex pressure/displacement
response correction, upward polarity, one-sided hourly Hann auto/cross-spectra,
linear/complex averaging, coherence and pressure-to-displacement transfer.
Separate run folders preserve existing products. Added restricted `paired_*`
master updates, station discovery, dependency-aware restart and saved-plot viewing.
Unresolved coverage/pairing differs from software failure; scores remain unchanged.

RR52's interior week retained 167/168 pairs. Median squared coherence in 0.1–0.2 Hz
is 0.80214, median |D/P| is 4.98579e−8 m/Pa and median wrapped phase is +4.119°.
Independent SciPy spectra and cross-spectra reproduced the first pair; independent
averages reproduced transfer/coherence. The figure was visually reviewed.
Production master hashes are unchanged. Temporary-export network tests preserve
existing cells/sheets and invalidate caches when the partner source changes.
The expanded suite passes 72 tests. Tilt correction and acceptance remain pending.
[Method and evidence](paired_spectra.md).

Investigated the existing Godin pilot errors without rerunning them: RR12's
21 selected windows and four sampled alternative interior windows are constant;
RR42's single pressure file supplies only about one day. No instrument failure
or correction follows from these observations.
[Bounded-read evidence](pilot_error_investigation.json).

## Step 19 — teleseismic calibration source review and RR52 comparison

Read the calibration/tilt sections of the thesis and manuscript, tilt milestone,
BRUIT-FM proposals, Zha/Webb and related calibration literature, original and
current ComPy/DPG code, both newer calibration packages and Wayne's TiSKitPy
rotation implementation. Distinguished repaired current code from historical
filtering/indexing/picking weaknesses and retained exact source revisions.

Implemented a standalone candidate runner with explicit catalogue/dispersion
inputs, group-velocity packet windows, quiet-only minimum-variance projection,
full complex response correction, pressure/upward-acceleration spectra and
analytically checked pressure-multiplier grid search. Generated the historical
model's phase/group table with existing disba 0.7.0; no environment change.

Eight mainshocks: three amplitude-review cases, three insufficient band-mean
coherence cases and two unequal sample grids. The amplitude-only median is
0.66449 versus historical 0.67, with three-event scatter 0.00722. All amplitude
cases have phase near the opposite convention, so no signed candidate passes.
The Santa Cruz case also has a predicted overlapping aftershock packet. No
response/master/score was changed. Both plots were visually inspected; model
sensitivity is quantified separately and is not a confidence interval.

Independent acoustic ODE, gain/reciprocal/grid, synthetic raw-response/rotation,
arrival and sign-reversal checks pass. The full suite passes 79 tests. The small
prototype's 10 checks and SeafloorCal's 47 focused model/guard checks also pass.
[Full source review and limitations](rayleigh_calibration_review.md).

## Step 19 follow-up — TiSKitPy default, provisional sign and window review

The user recalls opposite pressure/vertical signs and authorized the probable
convention while retaining the uncertainty. RR52 now defaults to pressure −1,
saved as `assumed_needs_confirmation`; one station-wide sign is used. Nominal
cross-phase is retained alongside provisional-sign phase. Source responses,
raw files, production masters and scores are unchanged.

The actual TiSKitPy CleanRotator 2.3.1 class is now the default. Only the optional
TiSKitPy/PyYAML/xarray dependencies were added; existing core packages were
preserved. Quiet training uses physical velocity units and anti-aliased ≤1 Hz
copies. The fitted vertical weights come from public `apply()`, not unvalidated
`tfs()`. The earlier local eigenvector method remains explicit opt-in with
`needs_to_be_tested` metadata and a [short reference](tilt_methods.md).

The new eight-event run has two candidates under the assumed sign, one coherent
overlapping-packet review case, three low band-mean coherence rejections and two
quiet-training timing exclusions. The three amplitude-review factors are
0.67648/0.66353/0.66440 (median 0.66440); the two candidates' median is 0.67044.
Neither is a finalized station calibration. Minimum three accepted independent
events is a pilot rule; station status remains unresolved. September BH1 offsets
are 0.0001 s, implying only 0.00216° phase at 0.06 Hz; no silent clock adjustment
was performed.

Added wide P/S/group-arrival/overlap plots, auditable manual-window decisions,
UTC bounds and source/override provenance. All three curated Santa Cruz plots
were visually inspected. Actual-library projection, preservation, provisional
sign/gain recovery, experimental-method flag, manual exclusion and timing
diagnostics checks pass. Full BRUIT-FM suite: **84 passed**, Ruff checks pass.
Original run/plots remain available as historical comparisons.
[Current validation](rr52_tiskitpy_validation_summary.json) ·
[How to inspect windows and open plots](rayleigh_window_review.md).

## Step 19 follow-up — original figure layout and magnitude fallback

Rendered the manuscript Figure 5 (PDF p. 16) and thesis Figure 3.13 (PDF p. 79),
and inspected the old plotting function and ComPy font settings. Added their
two-panel semilog 0.005–0.1 Hz layout, labels, green event curves, blue medians,
purple theory, black dashed gain-adjusted median and red fit limits. Broad display
spectra are separate from the fit/QC band and saved as `paper_*` NPZ fields.
Theory is not extrapolated beyond the prepared phase-velocity reference.
Both new curated figures were visually inspected against the originals.

Added Mw>7 → Mw≥6.5 → Mw≥6 fallback while fewer than three QC-accepted candidates
exist. Eight new attempts per level by default; no repeated events, cross-level
24-hour clusters or relaxed scientific checks. Search counts/reasons are saved.
An explicit event disables fallback; the original default fit band is preserved.

Real default-band verification tried 24 events, still yielding two candidates.
The separate original-figure 0.03–0.07 Hz run reached three candidates after four
Mw>7 attempts (one overlapping packet remains excluded). Individual gains are
0.679787/0.671255/0.678998, median 0.678998 and sample scatter 0.004715.
The aggregate fit to the median ratio is 0.674374, a different estimator.
No response/master/source correction was applied; production master hashes are
unchanged. Corrected an earlier prose date typo: `usb000g7x7` is April 16.
Full suite: **87 passed**. New tests verify fallback by QC count, inclusive 6.5/6
floors, stopping, deduplication/clustering, controls and aggregate synthetic gain.
[Method, plots and commands](rayleigh_paper_figure.md) ·
[Validation](rr52_paper_validation_summary.json).

## Step 20 — available-network Rayleigh comparisons

Added the explicit `network rayleigh` stage and a serial, checkpointed station
driver. It binds pressure/horizontal/vertical source epochs without first-match
selection, requires saved responses, and uses actual TiSKitPy plus the current
water-reference/fit method and bounded magnitude fallback. Network defaults use
the original-figure 0.03–0.07 Hz band and explicitly assumed pressure −1.
The shared historical RR52 model is a flagged network assumption, not a verified
local model for every station. Source/stat/response/model/code changes invalidate
reuse. Missing metadata and unsuitable resonance bands are saved unresolved;
only `rayleigh_*` diagnostics may be appended to the master, with backups.

A live YV run selected 43 available pressure epochs. It was interrupted at the
user's request to run it themselves after seven completed station comparisons.
RR03/RR07/RR08/RR09 each have three automatic candidates; RR01/RR05/RR06 remain
unresolved. Candidate event-gain medians are 1.12244/1.08223/1.15786/1.20592;
they remain provisional. The master hash checks are unchanged because no final
batch write ran. Existing PPSD/diagnostic/source products were preserved.
RR03's combined plot was visually inspected and copied as a curated example.
[Commands and partial run](yv_rayleigh_calibration.md) ·
[Validation](yv_rayleigh_validation_summary.json).

The user subsequently completed the full 43-epoch YV Rayleigh run: 16 candidate
comparisons, 25 unresolved and two unusable reference bands (RR41/RR57).
The final master update appended `rayleigh_*` diagnostics. A read-only comparison
verified that **all pre-existing CSV columns are exactly unchanged**. All six
scores remain blank and no response/waveform/PPSD correction was applied.
[Completed network audit](yv_rayleigh_completed_summary.json) supersedes the
earlier partial execution snapshot, which is retained as historical evidence.

## Step 21 — named assumption flags and validation checklist

Recorded implementation status as working with provisional assumptions, and
added named polarity, absolute-Z response, depth/model, water, geometry, quiet
screen/coda, tilt stationarity, packet/mode, QC, magnitude, scalar-gain and
uncertainty flags. Generic hardware, timing discrepancies and experimental local
rotation carry specific additional review flags where applicable.
Saved previously fixed pilot choices in station settings and documented what
needs study, testing or improvement in [the checklist](calibration_assumptions.md).
No numerical target was introduced into gain fitting and no score was assigned.

Full suite: **92 passed**; Ruff checks pass. New network checks exercise real
synthetic raw processing, preservation of original master columns and extra
Excel sheets, restricted scores, restart/source invalidation, pressure-epoch
ambiguity, missing channels, unusable reference bands and active-workflow guards.
Existing synthetic/ODE tests support software correctness under tested
assumptions; they do not certify station absolute accuracy.

## Step 22 — independent spectra and event/training review

Supplied the [teleseismic DPG equations and source report](teleseismic_dpg_calibration_report.md),
with a local frozen source archive for the completed-run version. Inspected all
48 candidate arrival plots and screened 48 quiet intervals; 46 passed continuity,
while RR11/RR25 had fractional-grid boundaries. Seven independent raw-count CSD
estimates match saved spectra/gains to roundoff. Twenty candidate windows contain
predicted body-arrival times, requiring packet review. Split-training and band/
model perturbations quantify illustrative sensitivity, not total uncertainty.
The current RESIF service repeats questioned positive pressure poles at RR03/RR22;
34 local pressure responses have such poles. Physical response conventions remain
unverified; no response was edited. Added and tested a calibration-only strict
record-grid guard; original results/masters were preserved. Full suite **98 passed**,
Ruff passed. [Detailed review and uncertainty budget](yv_rayleigh_independent_validation.md).

## Step 23 — optional teleseismic tilt correction

The user specified that teleseismic pressure calibration should be usable without
tilt correction. `--tilt-method none` (standalone) and
`--rayleigh-tilt-method none` (network) now select pressure and the matched
vertical response epochs only. They skip quiet-training reads, horizontal-channel
requirements and rotation. The original vertical is made upward-positive using
its dip and converted to acceleration. Existing timing, response, coherence,
phase, event-count and gain safeguards remain; both tilt methods remain available.
Old-command defaults are retained, so this comparison selects `none` explicitly.

Separate network commands use `--rayleigh-no-master-update` and a new workflow
directory. `--view-workflow` opens those saved comparison plots without depending
on historical master links. Unique station run folders preserve previous results.
No real-data batch or response correction was launched for this change.

Synthetic integration tests cover absent horizontals/quiet data, no TiSKitPy
dependency, upward/downward Z signs, known gain recovery, result provenance,
separate workflow viewing and master preservation. **101 tests pass**, Ruff
check/format pass, and all 2650 historical protected product hashes match the
previous audit. [Run/view commands](yv_rayleigh_calibration.md#comparison-without-tilt-correction).

## Step 24 — downloaded YL five-station pilot preparation

Recognize daily `NET.STA.--.CHA.YYYYMMDD.mseed` files alongside SDS names and
normalize only their filename location identifier for header matching. Preserve
raw filenames, bytes and clocks. Add a preparation script for A02W, A06W, A08W,
A10W and A12W and a deployment-specific USGS catalogue (566 saved events).
The user will run ten-channel PPSD prerequisites and calibration without tilt;
no real YL processing or master registration was launched by the agent.
Document borrowed-model/signed-response assumptions and observed header offsets
requiring timing review. [Commands and limits](yl_five_station_pilot.md).

## Next steps

Latest sampling correction: the user instructed use of all available eligible
earthquakes rather than a three-success target. `event_sampling="all"` is now the
default in standalone/network code; it ignores the old per-level attempt cap and
continues through the magnitude levels. Scientific epoch/magnitude/depth, 24-hour
separation, timing and QC filters remain. Minimum accepted events controls only
availability of a station estimate. Explicit `bounded` mode retains the historical
pilot behavior. Tests verify continued sampling beyond three successes and beyond
the configured attempt cap, without duplicate/cluster reuse. **109 tests pass**,
Ruff check/format pass. No real batch or model replacement was launched.

Model follow-up: cached historical dispersion provenance records disba 0.7.0.
Disba computes phase/group curves from supplied layers and uses Numba for JIT
acceleration; these are not Earth models themselves. At that checkpoint, the
environment read the cached CSV without either package installed. For YL reproduction, prioritize
numeric Lau Basin average dispersion/model from Zha et al. (2014), used by
[Zha & Webb (2016) Appendix A](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012688).
The 2014 study describes a network-average 1-D model; numerical source inputs
still need to be obtained and verified before replacing the provisional model.

Earth-model option follow-up: standalone `--earth-model-file` and network
`--rayleigh-earth-model-file` now accept explicit layers instead of a dispersion
table. Disba 0.7.0 and Numba 0.68.0 are installed, without changing existing
NumPy/SciPy. Generated phase/group curves and model/solver provenance are saved
and included in network reuse requirements. Tests independently check homogeneous
half-space Rayleigh velocities, model guards, known-gain recovery and both API
paths. **117 tests pass**, Ruff check/format pass. Original cached curves reproduce
to numerical precision. No real-data calibration batch or gain application ran.

The original four YL station reciprocals are about 5–6% below Zha's coefficients;
most of the apparent large difference comes from comparing reciprocal correction
directions. A saved-spectrum diagnostic at 0.04–0.1 Hz, holding the original
candidates and historical layers fixed, does not uniformly close that remaining
gap. Model, events, windows, estimator and response conventions need separate
checks. [Model option and comparison](rayleigh_earth_model_options.md).

All-station YL preparation: discover all original station XMLs (currently 25
stations / 50 pressure-vertical channels), register missing metadata once and
cache a full-deployment catalogue (568 events). Reuse completed PPSDs only when
inventory/raw source metadata still match; new or changed products use a fresh
root. A Bash runner stops on errors and performs preparation, required PPSD work
and all-event calibration without calibration master writes. Existing pilot
products remain untouched. **119 tests pass**, Ruff and Bash syntax checks pass.
The agent prepared metadata selectors/catalogue only and did not launch the
all-station registration or batch. [Commands](yl_all_station_run.md).

The user subsequently completed the all-station run: 22 candidates, three
unresolved (A14W/B01W/C10W). A read-only saved-gain comparison with published
Table 1 converts pressure multipliers to reciprocals: 21 stations are consistently
3.4–7.0% lower, while C07W is +86.8%. C07W's already accepted M>=7 subset retains
the discrepancy; its individual event gains span .259–.543. Nominal response
values alone do not explain it. B01W's published reciprocal pressure factor5.26
exceeds the current search upper bound3, and four saved attempts hit a boundary.
No new fit, correction, model/threshold/timing change or score was assigned.
Original workflow/master/result hashes were checked unchanged. Arithmetic/count
checks and Ruff passed. [Comparison](yl_all_zha_comparison.md).

End-of-session checkpoint: README completed steps now extend through 30, including
the saved-result outlier audit. The compact continuation note records the finished
YL workflow and next targeted checks. The comparison is provisional: 21 reciprocal
gains are 3.4–7.0% below Zha (median −5.4%), C07W is +86.8%, and three stations
remain unresolved. All source data, original results, masters and scores stay
preserved locally; code and curated comparison artifacts are tracked on GitHub.

Repository-wide correctness audit (README step 31): 130 tests pass. Read-only
checks cover 222 complete PPSD channels / 980,630 saved segments, 2,219 Rayleigh
spectra across 116 saved workflow runs, the other completed diagnostics, workflow
artifacts and all CSV/Excel cells. No saved arithmetic/provenance mismatch was
found. All 50 YL P/Z products cover their full available intervals. Protected
result hashes and current masters match; 93,784 baseline file stats stayed intact.
Five software corrections address bounded-PPSD reuse, pressure-unit aliases,
explicit-depth preflight, coordinate cache dependencies and auxiliary-fit failures
discarding primary candidates. Primary bounds/QC/model/timing rules and all saved
station gains remain unchanged. Scientific response/window/model validation is
still pending. [Audit and limits](project_audit.md).

Histogram follow-up: 684,714 saved values across 38 channels are outside the
stored −250/+100 dB bin range. RR12 pressure and all four RR50 channels contain
many numerical zero-power floor values near −3076.53 dB. These require quality
screening before accepted noise statistics; histogram edge bins collect the
under/overflow, while original saved arrays remain unchanged. The calibration
fits use raw event windows rather than PPSD amplitudes. Pushed fixes passed CI.

New network preparation (README steps 33–34): isolated 7D and XO masters,
header-named symlink views and complete fetched response epochs. All original
data/metadata/master files remain preserved; incompatible files are logged, never
resampled. 7D finished 30 PPSDs / 15 calibration jobs, with four provisional
estimates and 11 unresolved. No station-matched published gain benchmark verified.
XO's first attempt was interrupted during PPSD processing; the user resumed it.
Recheck live state before reporting calibration results. LA-series XO pairs are
EDH hydrophones / HHZ at 100 Hz, not their 125 Hz HDH APGs. WD pairs are LDH/LHZ
at 1 Hz. Separate regression checks cover identity preservation, mixed-rate
rejection, and refusal to resample APGs. [Commands and limits](7d_xo_runs.md).

YL first-run fix: accept finite nonzero signed response sensitivities while
preserving supplied stage gains/complex responses. The user's five vertical
PPSDs completed; pressure jobs had failed on the old positive-only guard. A
six-hour A02W pressure validation saved three PSD segments with unchanged
raw/XML/master hashes. Retry pressure channels only; timing and physical-model
review remain separate. [Details](yl_five_station_pilot.md).

1. Verify instrument frequency limits and define the six tests' acceptance criteria.
2. Review the initial pressure-level comparison and calculate the remaining
   tests, using the paired foundation; add tilt correction and source-specific
   interpretation before using results for quality classification.
3. Decide whether to add a separate shorter-window calculation for time variations.
4. Set network priorities, extend processing to more stations, and combine accepted
   spectra by depth, season, coastal distance and seafloor type.

**Current status:** the RR52 PPSD calculation, diagnostic plots and initial
expected-pressure comparison are complete.
The initial predicted-infragravity diagnostic is also complete.
The separate screened, full-window Godin pilot is complete.
The nearby-station diagnostic and its network stage are implemented and verified.
The paired spectral foundation and RR52 pilot are implemented and verified.
The teleseismic runner uses actual TiSKitPy and provisional RR52 pressure −1;
arrival isolation, timing policy and model/uncertainty verification remain.
Instrument calibration acceptance and combined ocean-bottom noise models remain pending.
An [initial saved-results review](yv_rayleigh_saved_results_review.md) rederived
all 481 spectral products and 48 candidate-event gains, inspected all 16 combined
candidate plots and sampled four arrival/waveform cases. Arithmetic matches;
packet isolation, response/model validation and uncertainty remain. RR03's sampled
window contains the predicted S arrival; RR57 has a conservative preflight
exclusion. Original statuses, gains, guards and masters were preserved.
The user has completed the network Rayleigh run. Review its saved products;
do not restart it autonomously. The available original chats and four-page WP2
brief were reconciled in [project context](project_context.md), including
multi-network/seasonal/global coverage, sensor weighting not yet adopted,
classification and final noise-model deliverables. The earlier percentage
estimate is not a measured whole-project completion fraction.

[Processing details](rr52_methods.md) · [Exact results and settings](rr52_run_summary.json)
· [References](references.bib)

## Numbered continuation: steps 25–40

The entries below organize the later checkpoints against the README numbering.
Earlier execution observations and test counts above remain historical evidence.

### 25. Accept signed YL sensitivities without changing polarity

Allow finite nonzero signed response sensitivities. Verify PSD sign invariance
and retain complex-response polarity. A bounded pressure pilot passed; original
raw/XML/master hashes and completed vertical PPSDs remained preserved.

### 26. Use all eligible events by default

Separate catalogue sampling from the minimum of three accepted events needed
for a station estimate. Historical attempt/success caps require explicit bounded
sampling. Existing magnitude, depth, clustering and event QC guards remain.

### 27. Add optional disba Earth-layer input

Generate fundamental Rayleigh phase/group dispersion from a layer CSV, retaining
model hashes and solver versions. No verified numeric Zha/Lau model has been
obtained; bundled historical RR52 layers are not Zha's regional model.

### 28. Prepare every downloaded YL station

Discover 25 station XML folders / 50 P/Z channels, cache the deployment catalogue
and reuse only unchanged full-archive PPSD products. Preserve original pilot
results and calibration masters. [Run settings](yl_all_station_run.md).

### 29. Compare the completed YL run with Zha & Webb

YL has 22 provisional gains and three unresolved stations. Compare our reciprocal
1/g to the published coefficient direction. Twenty-one estimates are 3.4–7.0%
lower, median −5.4%; C07W is +86.8%. [Full comparison](yl_all_zha_comparison.md).

### 30. Inspect saved YL outliers and rejection reasons

C07W's accepted M≥7 subset retains the discrepancy. Listed nominal sensitivities
do not explain its factor-of-two event grouping. B01W's published reciprocal
implies g≈5.26, beyond the search upper bound of three. Preserve original fits;
review response, window and policy assumptions before changing parameters.

### 31. Audit software and saved products

Check saved arrays, gains, workflow artifacts, source provenance and CSV/Excel
agreement. Fix bounded-PPSD reuse, unit aliases, explicit-depth preflight,
coordinate-dependent caching and auxiliary-fit failure handling. The historical
audit passed 130 tests and preserved all results. Retain the numerical-floor and
histogram-saturation findings before accepted noise statistics. [Audit](project_audit.md).

### 32. Add a read-only matched-event YL viewer

Compare C07W/C08W/C09W saved event spectra, exact-epoch responses, coherence and
phase, including QC-rejected labels and original arrival/waveform links. The
viewer writes separate review products without refitting or changing scores.
[Instructions](yl_matched_event_review.md).

### 33. Prepare isolated 7D and XO workflows

Use header-named symlinks, fetched full response epochs and separate masters.
The initial native-rate runners log incompatible files without rewriting them;
the later explicit XO processing-copy mode is documented separately.

### 34. Record the completed 7D results

Fifteen stations / 30 PPSD channels completed. J30A/J38A/J45A/J54A have
provisional estimates, each from three accepted events; eleven remain unresolved.
No verified station-matched published gain table was located. [Gains](7d_xo_runs.md).

### 35. Record the original completed XO run

Six stations / 12 P/Z channels completed. WD66/WD70 have provisional DPG gains;
the LA33/LA34/LT02 EDH hydrophone test and WD48 have no accepted events. Identify
the hydrophone/APG channel mismatch before publication comparisons.

### 36. Implement explicit APG processing

Select Paroscientific HDH APGs at nominal 125 Hz with HHZ at 100 Hz. Opt-in
working-copy UTC gridding and anti-alias resampling retain original data and
responses. A recorded 0.1 ppm tolerance supports near-125 Hz encoded rates.
Exact response snapshots remove the need to calculate PPSDs before calibration.
Synthetic checks cover gain/phase, aliasing, gaps, overlaps and epoch bounds.
The separate three-station APG run uses 2 Hz and no tilt.

### 37. Discover all XO gauges by hardware metadata

Prefer APG where identified, otherwise DPG; do not infer hardware from HDH alone.
WD66/WD70 have HDH DPGs as well as 1 Hz LDH DPGs. The all-station runner uses
1 Hz working copies, all eligible events and a new result root without PPSD or
master/gain updates. It opens saved calibration plots after completion.

### 38. Compare the completed APG and all-gauge runs

Both runs completed. The latest XO workflow has four provisional gains:
LA33 g=0.934275 / seven accepted events; LA34 g=0.962680 / eight;
WD66 g=0.909498 / sixteen; WD70 g=0.929889 / twenty-one. LT02 and WD48 remain
unresolved. APG 1/2 Hz gains differ by less than 0.00003%; WD66/WD70 change
+0.199244%/−0.007922% from their original DPG estimates, retaining all previous
accepted events and adding six/five. This demonstrates processing stability, not
external absolute accuracy. [Exact results](calibration_network_summary.json).

### 39. Search for the original 7D published scale factors

Ruan (2014) computes a station scale factor jointly with sediment structure.
Its Table 3 reports delay/thickness/depth, not gains. Publisher/author text,
related papers and indexed public repositories did not yield verified original
numeric scale factors. Request factors and conventions from the authors before
quantitative comparison. Bell/Ruan/Forsyth (2015) corrects the 2014 starting
sediment parameter c to 0.029, not 0.29. No saved model or gain was changed.

### 40. Consolidate the project and verify the current checkout

Four networks have completed calibration attempts on 89 distinct stations:
46 provisional estimates, 41 unresolved outcomes and two reference-band
exclusions. Repeated runs are counted once. The new reporting snapshot checks
saved medians, reciprocals and counts; workflow/result and checked master hashes
remain unchanged. Original and new XO scientific scores are blank. No science
batch, gain application or archive-wide audit was rerun for this checkpoint.

Current validation: **143 tests pass**, Ruff check/format and Bash syntax pass.
Git history and live processes were inspected; no preparation/calibration process
was active. The README, compact and historical handoff headers, session notes and
network report now share the current checkpoint. Code, documentation and curated
summary artifacts are tracked; raw waveforms, inventories, masters and full saved
products remain local under existing exclusions. [Current report](calibration_network_summary.md).

Absolute response, depth/dispersion, polarity, logger timing, event isolation and
held-out uncertainty remain open. The informal research-prototype rating is not
an acceptance score. YL supplies the strongest numerical publication comparison;
C07W remains an outlier. XO's body-wave/band-specific publication is not yet a
matched scalar comparison. No measured whole-project completion percentage exists.


### 41. Review Sentinel/Socrates and close the fresh-export race

Integrate the independent engineering fixes: retain event spectra/QC at a fit
boundary, derive tilt horizontal prefixes from the vertical, refuse fresh master
overwrite and lock cooperating append/update writers. Codex found that fresh
creation bypassed the lock, reproduced the race with a regression, and extended
locking to that path. No scientific response, threshold, gain or timing change.
[Approval and qualifications](reviews/codex_integration_review.md).

### 42. Evaluate all four preserved network runs

Reproduce 416 saved accepted-event fits and 46 station estimates over 89 stations.
Check 18 raw event windows, one B01W end-to-end boundary case and two C07W raw
daily states in isolated review products. 155 tests pass; Ruff and runner syntax
pass. Protected source/result/master hashes stay unchanged within the recorded
bounded scope. No full batch, gain application or score update.

The earlier B01W bound interpretation is incomplete: all four coherent boundary
events have zero phase-valid bins, so larger bounds alone cannot resolve them.
C07W's two relative vertical states are supported; its exact hardware cause is
unconfirmed and peer comparisons need their own QC. RR55 is not established to
share B01W's high-pass mechanism. Aminian (2025) supplies eight published YV
pressure multipliers; the original tilt baseline differs by +0.3–6.1%, median
+3.0%, without establishing a causal glitch bias or absolute acceptance.
Pole reflection can change accepted bins and therefore gains even with unchanged
response amplitudes. No production pole correction was made. IID daily bootstrap
intervals omit serial dependence. [Curated checks](reviews/codex_integration_validation.json).


### 43. Implement optional scientific correction trials

Codex implemented Socrates-informed pole reflection and weighted log fitting as
separate cases on saved no-tilt event windows, alongside guarded C07W state subsets
and diagnostic-only gain/corner fits. Sentinel's merged QC-preservation fixes
retain boundary spectra. Original production defaults and source products remain
intact. Seven-station trials reproduce baselines and isolate each proposed change.
[Scope and commands](socrates_selected_station_trial.md).

### 44. Evaluate all previously tested stations

The user completed all 89 YV/YL/7D/XO trials. YV estimates increase 16 → 31 under
pole reflection; combined trial count is 61 including C07W's three-event normal
subset, with 26 ordinary unresolved cases and two reference exclusions. Weighted
fitting worsens the YL peer median discrepancy from −5.43% to −6.72%; YV's eight
matched INSU comparisons change from +3.20% to +3.50%. C07W state separation
reduces its discrepancy, without confirming the hardware mechanism. All 2,432
protected hashes remain unchanged. [89-row comparison](socrates_all_station_comparison.csv)
and [interpretation](socrates_all_station_comparison.md).

### 45. Integrate reviewed work and record contributions

At the user's explicit request, merge the optional trial tools and comparison
records into main alongside the existing reviewed engineering fixes. Retain the
original fitting objective as the normal calibration default. Software validation:
167 tests passed, Ruff check/format and Bash syntax passed. No further science
batch or gain/response/master/score application. Record the project owner's
scientific leadership and real runs, Sentinel's independent engineering review,
Socrates' scientific review, and Codex's implementation/verification/integration.
[Contribution record](contributions.md). Further physical validation remains open.


### 46. Prepare unused-event frozen-gain validation

At the user's request, add a four-station RR01/RR06/RR48/RR52 validation runner.
Freeze the completed pole-only trial gains and exclude all recorded YV trial
IDs plus a separately documented native development probe. The preflight selects
50 unused Mw≥6 mainshocks under the existing depth/24 h rules. Report fixed-gain
frequency/event errors without refitting, and separate predicted arrival review
from the isolated-window summary. Record inventory metadata and pending physical
hardware identity checks. Native code checks preserve sources; 173 tests pass.
This implements a validation study, not accepted calibration; the user will run
its full four-station workload. [Command and scope](yv_heldout_validation.md).

### 47. Interpret the completed four-station held-out study

The user attempted all 50 selected events at each station with frozen gains.
RR52 has 14 isolated usable events, median prediction error +0.664%; RR06 one
(-1.041%); RR01/RR48 none. RR48 is chiefly limited by its June epoch end.
All 71 protected hashes and 2,208 raw-file stats match. These are conditional
reference-model prediction errors, not independent absolute calibration.
[Completed results](yv_heldout_results.md).

### 48. Review flagged arrivals and test a coverage-based replacement

Codex inspected three bounded native RR01/RR06 P/Z waveform overviews at the
user's request. Strong calibration-band packets follow the predicted body-wave
markers; keep original arrival flags because visual inspection does not prove
mode isolation. Add a separate display-only plotter; no gain refits or application.
Screen four long-epoch candidate stations using headers before new errors.
RR05 has 42 epoch-valid covered windows versus RR48's nine. Its isolated new
50-event run, frozen g=1.12014351665, has zero QC-usable events: 40 coherence
failures and ten input/window refusals. Longer coverage does not overcome all
signal-quality limitations. Preserve all results and thresholds; investigate
frequency-dependent signal/noise and prospective packet placement next.
Ruff check/format, three real-window plots and six held-out tests pass. All
recorded hashes/stats match, including a recheck of the prior four-station run.
[Plots, coverage, result and preservation](yv_heldout_arrival_review.md).

### 49. Close the calibration checkpoint and prepare the microseism stage

The user requested a fresh-chat handoff and chose the second pressure–vertical
check as the next stage, using RR52's saved paired spectra as a pilot. Clarify
that 61/89 optional-trial estimates (46 baseline) are provisional, 7D/XO gains
were computed despite missing comparable published coefficients, and publication
agreement is not absolute accuracy. A read-only same-event check of
`usb000kdb4` finds mean coherence² RR05 0.941675, RR06 0.982187 and RR52 0.098078;
RR52 cannot be assumed a successful control for every earthquake. No refit.

Record the exact existing RR52 paired run, physical units/transfer conventions,
primary-source leads and open identifiability/model/tilt questions. Leave further
Rayleigh troubleshooting deferred. No new microseism calibration, raw processing,
master edits or scores were performed during this documentation-only transition.
[Next-chat brief](microseism_next_chat.md). Local-P checks, classification and
combined low/high ocean-bottom noise models remain later stages.

### 50. Establish the RR52 microseism diagnostic and its identifiability limit

Review Ruan et al. (2014), Bell et al. (2015), Saikia et al. (2020), and the
exact OBSIC quality-assessment description. The microseism observable is the
complex upward-displacement/pressure admittance. Its amplitude jointly contains
pressure gain, vertical gain and coupled water/seafloor structure; published
DPG prefactors assume the vertical response and fix water depth/model constraints.
It is therefore not an independent pressure-gain recipe.

Add a read-only saved-spectrum analyzer. It verifies the designated RR52 pilot's
five artifacts and recorded source/response provenance, exactly reproduces the
weekly aggregate, re-averages the 167 saved hourly spectra into seven daily
products, and writes a fresh derived result. Median squared coherence remains
0.802142 in 0.1–0.2 Hz; 29/361 bins exceed the Ruan 0.95 reference and six of
seven daily means exceed the later Saikia 0.8 reference. Those rules are not
adopted scores because this one-week uncorrected pilot does not reproduce the
published tilt, selection, duration or modal inversion.

Keep RR52's Rayleigh pressure multiplier fixed at 0.6787278140712512. Dividing
the nominal transfer by it changes the weekly median amplitude from
4.9857904e-8 to 7.3457877e-8 m/Pa; no value is fitted or applied. Preserve raw
data, saved paired/Rayleigh/validation products, masters and blank score. Three
new diagnostic unit tests plus the 11 paired-spectra tests pass. The complete
suite passes 176 tests; Ruff check and format checks also pass.
[Theory, result and next constraint](rr52_microseism_diagnostic.md).

### 51. Prepare the baseline-station microseism campaign

At the user's request, extend the bounded diagnostic to every station with a
saved first-method Rayleigh estimate, leaving the user to execute the long
station jobs. Freeze the 46 `previous_g` values and exact P/Z channel pairs in
an attributed manifest: YV 16, YL 22, 7D 4 and XO 4. Exclude pole-trial-only,
unresolved and reference-band-excluded stations from this like-for-like scope.

Add a one-station runner that reuses RR52's designated saved paired spectra and
creates fresh paired/diagnostic directories for other stations. It never enables
master updates. Most pairs use native equal rates; XO LA33/LA34 use the already
audited explicit 125/100 Hz to 1 Hz working-copy resampling, and XO WD66/WD70 use
the same guarded 1 Hz path. All 46 exact saved response pairs reproduce on the
0.10/0.15/0.20 Hz check grid. Publication metadata preserves `g` versus `1/g`
direction and explicitly withholds unmatched 7D/XO scalar claims. No long data
job, gain fit/application, score or master write was performed while preparing
the commands. Five new tests bring the complete suite to 181; Ruff check/format
pass. [Per-station commands and interpretation](microseism_baseline_run.md).

### 52. Compare the completed baseline microseism campaign

The user completed the explicit 46-station run: YV 16, YL 22, 7D 4 and XO 4,
with no failed station. Add a read-only aggregator that requires the exact
campaign log directory, binds every station to the diagnostic path printed in
its log, checks all 462 recorded protected hashes/paths, verifies recorded raw
metadata and response hashes, and reproduces every saved aggregate from the
stored window spectra. It reads no new waveform and changes no prior artifact.

Weekly median coherence² is at least 0.8 at 37/46 stations; 30/46 have at least
one 0.1–0.2 Hz bin above the Ruan 0.95 reference. The Saikia daily mean
reference selects 149/275 station-days at 36/46 stations. These are descriptive
references, not new scores. Coverage is unequal (24–168 hourly segments) and
the automatically selected weeks are not synchronized. The strongest weekly
median is YL.C05W (0.982); the weakest is YV.RR11 (0.583). RR52 retains its
pilot values: weekly median 0.802142, 29/361 Ruan-reference bins and 6/7 selected
days.

All 16 YV phases lie in the near-0° family (network circular mean −3.5°), while
all 30 YL/7D/XO phases lie near 180°. Positive scalar gain cannot produce that
split, so preserve it as an unresolved network sign/response-convention issue
before any common signed modal inversion. Displaying the frozen Rayleigh gains
reduces the within-network geometric amplitude spread for YV (×1.579→×1.342),
YL (×1.629→×1.587) and XO (×1.484→×1.459), but increases it for 7D
(×1.292→×1.362). This is not independent validation because water depth,
structure, vertical response and wavefield remain uncontrolled.

Retain the prior paper comparison honestly: eight YV first-method values differ
from Aminian et al. by +0.35…+6.08% (median +3.20%); 21 YL reciprocal peers
excluding C07W differ from Zha & Webb by −7.04…−3.37% (median −5.43%). No new
microseism coefficient exists to compare. Ruan requires modal reproduction and
Deng's XO body-wave bands are not scalar-equivalent. Five comparison tests bring
the complete suite to 186; Ruff check/format pass.
[Comparison report](microseism_baseline_comparison.md).

### 57. Audit the frozen microseism gain bands and stage table results

Following the integrated method and repeated-week work in README steps 53–56,
audit all seven conditionally repeatable stations on their exact eligible bands.
Reuse saved complex spectra and reproduce the frozen selections: 22 weeks and
58 bands, of which 55 meet the legacy 10-degree phase screen and three RR40
bands are approximate. Retain five inconsistent RR52 daily cases and sparse
7D coverage. Do not retrospectively reselect data or change the frozen rule.

Exact-epoch complex pressure/vertical responses reproduce from saved snapshots.
This checks processing consistency, not absolute sensor calibration. Cached
water/dry model choice shifts 7D band medians by <0.6%, so it alone does not
explain the 4–6% Rayleigh discrepancy. Coherent tilt remains unidentifiable from
P/Z alone. Header-only preflight confirms local YV horizontal channels for the
next bounded pilot; no new waveform processing or calibration batch is run.

Export 89 unique full-epoch staging rows, with seven conditional/unvalidated
gains and all scientific scores blank. Preserve pairing and phase reasons,
unresolved assumptions and full input provenance. No master or original result
is changed. [Audit, evidence and next pilot](microseism_gain_band_audit.md).

### 58. Run a bounded RR52 horizontal-coupling diagnostic

Implement native-rate P/H1/H2/Z cross spectra and conditioning-safe Schur
complement diagnostics. Tests reproduce the existing pair normalization,
check covariance/units and refusals, and demonstrate that physical seismic
motion and finite-sample projection bias can imitate a large H-associated
component without tilt. The full suite has 262 passing tests; Ruff passes.

Two target days in RR52's original week retain 46/48 hourly windows. All retained
P/Z spectra reproduce their saved counterparts. Preserve 0.0024-second midnight
grid refusals and duplicate horizontal azimuth metadata, rather than forcing
alignment or geographic axes. Low-band Z predictability from H is substantial
(67.5%/96.9%) but is not an estimate of tilt power or a pressure correction.
No gain or cleaned data are produced; 2,498 protected inputs match. Reuse saved
four-channel coefficients for any continuation, not another raw run.
[Results, provenance and remaining constraints](rr52_microseism_horizontal_pilot.md).

### 59. Resolve the bounded orientation and separate-training review

Visually verify RR52's row in Scholz et al. (2017), Table 1, and record its
29-degree standardized orientation with separate P/R measurement counts and
reported spreads. Preserve the native/standardized handedness caveat; no XML
azimuth, geographic rotation or physical tilt is inferred from that number.

Select the latest completed historical quiet fit before the first June test,
using timestamps alone, and test its fixed relative-basis weights on the saved
46 hourly four-channel spectra. Retain all five prior fits as declared primary/
sensitivity cases; do not select by gain or event QC. Source epochs, physical
velocity training units, response hashes, original spectra and metadata match.
The June intervals were previously inspected, so no blind-validation claim is made.

Low-band power reductions do not translate into uniformly improved coherence:
highest-band Z power increases ~8-9%, with lower coherence on both days. Do not
adopt a correction or rerun batches to force agreement. Complete the bounded
task with conditional/unresolved results; proceed to local-P/status integration.
Fourteen new tests give 276 passing tests; Ruff passes and 2,545 protected inputs
remain unchanged. [Method, numbers and decision](rr52_microseism_fixed_projection.md).
