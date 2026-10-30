# Selected-station scientific correction experiment

The user-approved optional trial tools are integrated into `main`, retaining the
history of `experiment/socrates-selected-stations` and engineering main `0943889`.
Production calibration defaults are unchanged. The user completed the 89-station
trial; use its [saved comparison](socrates_all_station_comparison.md) for review.
The commands below reproduce optional experiments in new result folders.

```bash
cd /Users/mohammadamin/Projects/bruit-fm
PYTHONPATH=src conda run --no-capture-output -n bruitfm python \
  scripts/run_calibration_experiment.py --network both --open-plots
```

This selects YV **RR03, RR22, RR52** and YL **A02W, B01W, C07W, C08W**.
Every stored spectral event at these stations is included; there is no new
three-event cap. This is a fixed-window comparison of the earlier calculations,
not a new all-catalogue discovery or full-network batch. Original timing, data,
clustering and overlap exclusions are retained. Legacy boundary events whose
spectra were not saved are reconstructed with the original no-tilt/raw-window
guards. No PPSD preparation is necessary.

For one network only:

```bash
PYTHONPATH=src conda run --no-capture-output -n bruitfm python \
  scripts/run_calibration_experiment.py --network YV --stations RR03 RR22 RR52
PYTHONPATH=src conda run --no-capture-output -n bruitfm python \
  scripts/run_calibration_experiment.py --network YL --stations A02W B01W C07W C08W
```

## All previously tested stations

```bash
cd /Users/mohammadamin/Projects/bruit-fm
PYTHONPATH=src conda run --no-capture-output -n bruitfm python \
  scripts/run_calibration_experiment.py --network all --all-stations --open-plots
```

This discovers **89 station epochs** from the four preserved workflows:
YV 43, YL 25, 7D 15 and XO 6. Use `--network both --all-stations` for only
YV/YL (68 stations), or select one network with `--network YV`, `YL`, `7D` or
`XO` plus `--all-stations`. The existing small `--network both` default is kept.
RR41/RR57 retain their reference-band exclusions; no guard is relaxed.
XO uses its latest all-gauge APG/DPG baseline and recorded working-copy rate.
Publication comparisons cover all 25 YL coefficients and eight matched YV DPGs.
No verified numeric 7D or directly comparable scalar XO benchmark is inserted.
Each run produces a new result folder and network panels in `comparison.png`.
Masters are loaded once per source folder and reused read-only across stations.

## Options evaluated separately

| Case | Pressure convention | Gain objective |
|---|---|---|
| `nominal_absolute` | Original complex response | Original absolute-unit least squares |
| `pole_absolute` | Reflect real positive YV pressure poles | Original least squares |
| `nominal_weighted_log` | Original response | Coherence-weighted mean ln(T/R) |
| `pole_weighted_log` | Reflect real positive YV pressure poles | Weighted ln(T/R) |

Pole reflection is an explicitly labelled alternative response hypothesis based
on [Stähler et al. (2016), Appendix A, Eq. A5](https://adgeo.copernicus.org/articles/41/43/2016/adgeo-41-43-2016.pdf).
In the small trial it affects RR03/RR22; the all-station trial evaluates every
real positive YV pressure pole. RR52/YL/7D/XO are controls for this option.
The factor R_nominal/R_alternative acts on the pressure–acceleration cross
spectrum, leaving PSDs/coherence unchanged and rechecking phase. No XML or
snapshot is altered. Pressure sign, window, reference model and all guards are
retained. Improved phase alone is not verification of per-instrument hardware.

The weighted log sensitivity fit uses approximate transfer variance
(1−γ²)/(2Kγ²). Coherence is capped at 0.999 to prevent singular weights.
This does not establish estimator-noise correction or confidence intervals;
Hann-bin dependence, channel-noise partition and model errors remain. The
highest accepted bin is also dropped in a separate sensitivity diagnostic.

**C07W:** classify event-window 0.12–0.3 Hz raw-count vertical RMS ratios against
C08W/C09W. Retain a separate normal-state subset and report its median only
with at least the existing three candidate events. The mixed-state median stays
diagnostic. Classification thresholds are explicit trial assumptions, not proof
of a hardware fault; no vertical amplitude is multiplied by two.

**B01W:** jointly fit high-pass gain/corner to amplitude and phase as a diagnostic.
It cannot override phase QC, create an accepted scalar gain or extrapolate a
coefficient to Zha's wider band. Insufficient coherence or poor residuals make
the model unconvincing, even if numerical optimization converges.

INSU periodic-glitch removal, quiet-noise subtraction, new regional models,
density changes and timing-tolerance changes are **not** part of this trial.
They require separate experiments rather than an unexplained combined change.

## Outputs and publication comparison

Each invocation creates a new folder under `results/EXPERIMENTS/SOCRATES/` and
prints its exact path. An existing `--output-dir` is refused. Outputs include:

- `experiment.json`: completion status, input hashes and limitations.
- `comparison.csv` / `comparison.png`: saved-before, replay and each trial case.
- `<network>_<station>/trial.json`: every event, QC/bins, response-shape diagnostics,
  C07W states and per-case station medians.
- Separate recovered NPZ files for legacy boundary events only.

Original baseline/master/inventory/snapshot/NPZ/model hashes are checked before
and after. Raw files actually read by recovery or C07W state diagnostics are also
hashed. This is a bounded preservation check, not a full raw-archive audit.
Original waveforms, responses, masters, saved results and all scores stay intact.
No calibration gain is applied to any waveform, response, PPSD or master.

YL coefficients are compared as **1/g** with Zha & Webb (2016) Table 1. YV
published values are **g** from Aminian et al. (2025) Table 1: RR52 has a matched
benchmark 0.67; RR03/RR22 have no verified numeric benchmark in that table.
The old no-tilt baseline and publication processing differ. Smaller numerical
discrepancies alone do not prove the scientific correction is valid.

Primary tables: [Zha & Webb (2016)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012688)
and [Aminian et al. (2025)](https://www.ipgp.fr/~stutz/2025-aminian-crawford-stutzmann-al-GJI.pdf).
See the [independent integration review](reviews/codex_integration_review.md)
for the earlier merged software fixes and qualifications.


## Local code verification before the user run

Full suite: **162 passed** (77 existing warnings); Ruff check/format passed.
Seven new analytic tests cover response phase/amplitude, pole-unit refusal,
weighting/unit invariance, phase/coherence/overlap precedence, retained bounds
and high-pass recovery without scalar acceptance. CLI dry-run finds all seven
selected station records.

The RR52 control completed end to end in an isolated folder: nominal and
pole-only g=0.6787278141 reproduce the saved result exactly; weighted log fitting
g=0.6840664793 shifts +0.787%, moving farther from the published 0.67. This is
sensitivity evidence, not a reason to tune the fitter toward the publication.
One B01W raw boundary recovery retains zero phase-valid bins; its diagnostic
corner is 0.07514 Hz, without an accepted scalar gain. C07W's minimum-g event
is classified reduced, while the maximum-g probe remains unclassified because
its neighbor-control ratio falls outside the declared trial range. No thresholds
were adjusted to force that event into the normal subset.
These bounded checks do not replace the user's complete selected-station trial.


All-network extension checks: 167 tests pass. A dry-run discovers all 89 epochs;
nominal replay over 1,997 cached spectra reproduces the 416 prior candidates
without station-count mismatches. Isolated real-data XO LA33 and 7D J30A runs
complete with source hashes unchanged. The user subsequently completed the all-station trial, and its optional tools
were integrated into main by explicit request. Physical calibration acceptance
and automatic response correction remain separate, unfinished work.
