# YV Rayleigh: independent raw spectral checks and event/quiet review

The calibration equations and source snapshot are supplied in the
[teleseismic DPG report](teleseismic_dpg_calibration_report.md). This follow-up
checks the completed YV run. **Spectral computation is consistent under the
supplied responses; field calibration acceptance remains pending.** No original
workflow status, fitted gain, master field, waveform, response or PPSD was changed.

## Independent raw estimator

Seven bounded raw-count event windows were read directly with ObsPy, independently
of the production archive helpers. SciPy CSD formed the full count cross-spectral
matrix first; each entry was divided by the corresponding complex response pair,
then the saved TiSKitPy weights projected velocity and acceleration. This differs
from production's response-corrected segment-FFT pathway. Its sign/density/
rotation behavior is tested with a known signal and complex count gains, without
using the water-reference function to generate that test.

Response metadata, evalresp and saved rotation are still shared. The comparison
therefore checks extraction, averaging, spectral density, cross-phase and
projection algebra, not absolute sensor sensitivity or the validity of those
metadata. Maximum representative PSD relative differences were about 5e-15;
phase differences were below 3e-14 degrees. Candidate gains match to roundoff.

| Station | Event | Independent candidate gain | Saved status |
|---|---|---:|---|
| RR03 | usb000ief9 | 1.166216 | candidate |
| RR07 | usc000f40j | 1.151866 | candidate |
| RR09 | usc000f40j | 1.234816 | candidate |
| RR22 | ak0138esnzr | — | polarity_review |
| RR40 | usc000f40j | 0.805504 | candidate |
| RR52 | usc000f40j | 0.671255 | candidate |
| RR55 | usb000g7x7 | — | polarity_review |

## All 48 candidate windows reviewed

All 48 saved broad arrival plots were visually inspected alongside P/S timings
and per-event metadata. **20 windows contain a predicted first P or S arrival**;
every one of the 16 candidate stations has at least one such window. This timing
flag is evidence requiring review, not proof of the degree of contamination.
The remaining 28 windows start after the displayed first arrivals, but body-wave
coda, modal mixtures and local structure still require checks.

Many distant Alaska examples contain a broad, structured wave train; Santa Cruz
examples show a later strong packet, often after the narrow historical group
prediction. Papua and September regional examples overlap early arrival timings.
Pressure and upward acceleration commonly track each other, but high coherence
does not distinguish a pure fundamental Rayleigh packet from other coherent motion.
The current selector has no explicit minimum-distance criterion. Magnitude alone
cannot establish teleseismic suitability.

The original 1800-second window/minimum-eight-200-second-segment choices limit
how much a short regional packet can be isolated. A revised estimator/window
policy requires scientific validation; this review did not shorten windows,
relax thresholds, or select a new window to reproduce the desired gain.
The sidecar records observations for every candidate; original statuses remain.

## Quiet training and native record boundaries

All 48 selected quiet intervals were investigated with independent MiniSEED reads.
**46** passed continuity screening and have twelve 10-minute count-variance metrics
and first-difference screening metrics. Exact aligned duplicate overlaps are
allowed only after verifying identical sample values. These are screening
metrics, not quietness acceptance cutoffs or physical-unit channel comparisons.

Some intervals are nonstationary. For RR52/Papua, horizontal count variance peaks
at about 63/83 times the respective interval median block variance; Z changes much
less. RR07/Santa Cruz has an approximately 12-fold pressure variance peak.
Counts across different channels cannot be interpreted as equal physical noise.
Representative screening plots were saved in the local validation folder.

Two September `usb000jyiv` pressure training intervals have fractional native-grid
boundaries: **RR11, 0.0016 s (0.08 samples)** and **RR25, 0.0011 s (0.055 samples)**
at 50 Hz. Their 0.07 Hz phase bounds are only 0.0403/0.0277 degrees; they cannot
explain large phase disagreements. The previous ObsPy merge returned a continuous
trace despite the boundary offsets. A new **calibration-only strict-grid guard**
now rejects this before merging and records a timing-review reason. It does not
interpolate or alter source records. RR13's initially flagged overlap was verified
as aligned identical duplicates and is accepted by the independent reader.

The generic reader's default remains unchanged; Rayleigh `load_aligned()` opts
into the new guard. Existing PPSD/diagnostic computations were not repeated.
The network cache includes changed scientific code, so a future explicit run can
recompute into separate folders. Do not resume the full batch just to make the
historical snapshot match the new guard.

## Tilt-training sensitivity

The actual TiSKitPy full training and its first/second one-hour halves were
compared for three representative candidates. Same deconvolution/prefilter,
anti-aliasing and constructor bands were used. Gains use the **original accepted
frequency bins**, so these are sensitivity diagnostics, not newly accepted fits.
Full-training reproduction matches saved angles/gains.

| Case | Full-training gain | First-half gain | Second-half gain |
|---|---:|---:|---:|
| RR52/usb000g7x7 | 0.678998 | 0.678992 | 0.679032 |
| RR03/usb000ief9 | 1.166216 | 1.163592 | 1.166051 |
| RR07/usc000f40j | 1.151866 | 1.150880 | Projection guard rejects |

RR52's tested gain is stable to these halves despite the count-variance change.
RR03's fitted axis changes more (about 0.88 degrees full versus 4.05/1.31 degrees
for the halves), while gain changes by about 0.23% at most in this test. These
small gain changes do not certify preservation of every event's Rayleigh motion.

RR07's second-half public-library projection has norm 0.999951 at a fitted angle
of 9.86 degrees and fails BRUIT-FM's unit-projection safeguard. It was neither
normalized nor used for a new gain. The unit-norm policy and TiSKitPy's two-axis
rotation approximation need study; this is not an instrument-failure conclusion.
Other rotation/projection rejections in the original run remain documented.

## Pressure responses and external metadata

At **34 of 43 stations**, the supplied pressure response has a positive-real pole
in a `LAPLACE (RADIANS/SECOND)` stage. For example RR03 has +0.047357 and RR22
+0.099011 radians/s. A conventional stable causal forward instrument transfer
normally has poles in the left half-plane, so these signs require instrument/
metadata verification. Do not silently negate poles: this changes complex phase
and requires a supported forward-response convention.

The original inventory identifies RESIF as its source. Bounded current RESIF
queries for RR03/RR22/RR52 reproduce those pressure poles, sensitivities and
Hydrophone/DPG labels. They are therefore not shown to be a local transcription
error. Archive consistency is not independent physical validation. The first
IRIS queries returned no metadata because RESIF serves these channels.
[RESIF RR03 response](https://ws.resif.fr/fdsnws/station/1/query?net=YV&sta=RR03&cha=BDH&level=response&starttime=2012-10-05&endtime=2013-04-20&format=xml),
[RESIF RR22 response](https://ws.resif.fr/fdsnws/station/1/query?net=YV&sta=RR22&cha=BDH&level=response&starttime=2012-10-05&endtime=2013-04-20&format=xml),
[RESIF RR52 response](https://ws.resif.fr/fdsnws/station/1/query?net=YV&sta=RR52&cha=BDH&level=response&starttime=2012-10-05&endtime=2013-04-20&format=xml).

Large integer FIR coefficient sums are also present at many stations; evalresp
reports normalization. They may be an intentional coefficient convention with
stage gains, rather than failed hardware. Full digitizer transfer/sensitivity and
normalization conventions still need verification. These warnings were retained;
no coefficient, pole, response file or source data was edited.

## Bands, model sensitivity and uncertainty budget

Alternative narrower fit bands and illustrative +/-2% sound-speed and +/-5%
phase-speed perturbations were evaluated from saved spectra, holding original
selected events and phase/coherence masks fixed. This tests sensitivity, not new
field acceptance. A station diagnostic stays blank when fewer than three events
retain at least three accepted bins. A perturbation reaching the resonance guard
is also left unavailable.

| Station | Original median gain | 0.03-0.06 Hz | 0.035-0.065 Hz |
|---|---:|---:|---:|
| RR03 | 1.1224 | Insufficient retained bins | Insufficient retained bins |
| RR07 | 1.0822 | 1.0407 | 1.0537 |
| RR40 | 0.7419 | 0.7538 | 0.7556 |
| RR52 | 0.6790 | 0.6717 | 0.6726 |

No true parameter uncertainties or independent response calibrations are supplied,
so these perturbations cannot be combined into a defensible numerical confidence
interval. The uncertainty budget must include:

| Contribution | Current evidence / missing quantity |
|---|---|
| Event sampling/selection | Only first three successes; scatter omits selection bias; held-out events needed |
| Window and mode | All candidates inspected; 20 body-arrival timing flags; pure-mode isolation unverified |
| Tilt training | Three split-training cases; partial stability evidence; full deployment stationarity unknown |
| Depth | Inventory estimate; +/-5% perturbations move gains roughly -6% to -12% / +7% to +19%; actual error unknown |
| Water properties / phase model | Assumed values; exploratory perturbations saved; appropriate measured/path constraints absent |
| Vertical and pressure response | Absolute sensitivity, pole and FIR conventions unverified |
| Timing / estimator / response shape | Small clock boundaries quantified; seven CSD comparisons agree; scalar frequency independence unverified |

RR41/RR57's conservative upper-band guard and a possible 0.03-0.06 Hz comparison
are discussed in the [initial review](yv_rayleigh_saved_results_review.md).
No replacement band was silently adopted or new station fit invented.

## How to reproduce/view this review

```bash
cd /Users/mohammadamin/Projects/bruit-fm
MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/bruitfm-mpl \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  conda run --no-capture-output -n bruitfm python scripts/validate_yv_rayleigh.py
```

The script writes a new `run_*` subfolder under `results/VALIDATION/YV_rayleigh/`;
it does not write the master. It reproduces raw-count and quiet screening plus
recorded window timing. The separate split-tilt/model and external-metadata
follow-ups are preserved in the curated audit with their methods described above.
The two resonance-excluded stations' original metadata were inspected separately.

Already-saved screening figures from this review can be viewed directly:

```bash
open results/VALIDATION/YV_rayleigh/RR52_usb000g7x7_quiet.png
open results/VALIDATION/YV_rayleigh/RR07_usc000f40j_quiet.png
```

**Verification:** 98 tests passed; Ruff check and format passed. Master CSV/Excel
hashes match before/after the review. Tests confirm software behavior under tested
assumptions; they do not resolve the physical uncertainties above.

Next prioritize original response conventions and event isolation, then held-out
sampling and measured uncertainty before applying gains. Microseism/local-P checks
and final acceptance/classification remain separate unfinished tasks.

The [targeted response/window follow-up](yv_response_window_validation.md) checks
four analog stages explicitly and shows that none of the 20 timing-flagged windows
can retain 1600 seconds while avoiding first P/S arrival instants within its saved
bounds. No response correction or replacement window is adopted.

[Curated numerical audit](yv_rayleigh_independent_validation.json) ·
[Calibration equations/code](teleseismic_dpg_calibration_report.md) ·
[Remaining assumptions](calibration_assumptions.md).
