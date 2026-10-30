# Paper Figure 5 / thesis Figure 3.13 comparison

The original manuscript Figure 5 (PDF p. 16) and thesis Figure 3.13 (PDF p. 79,
printed p. 57) were rendered and inspected alongside the original
`DPG_Calibration/Pressure_calibration.py` plotting function and ComPy styling.
BRUIT-FM now saves `paper_figure5_comparison.png` in each Rayleigh run.

## Matching the original presentation

Both panels use a logarithmic **0.005–0.1 Hz** x axis. The upper panel retains
`High-Coherence Teleseismic Events`, the `Coherence` label and 0–1 labeled ticks.
Green curves show individual QC-accepted events; blue shows their median;
red dashed lines mark the actual fit band. Coherence here is magnitude-squared,
as in the original code. No extra smoothing has been introduced.

The lower panel retains `Frequency [Hz]`, `Spectral ratio / ρh`, the −1–10
vertical range and the gain-factor title. Green curves show measured event
ratios; blue is their median; purple is the homogeneous-water prediction;
black dashed is the gain-adjusted median ratio. The STIX font family, panel
letters, legend labels/positions and colors follow the original presentation.
A small footer identifies the candidate count, assumed sign and provisional
status. The curve labeled `Corrected P/a ratio` is a plotted comparison only;
it does not modify any pressure data or response.

The figure fits a gain to the **median event spectral ratio**, reproducing the
original aggregate objective. This differs from the station summary's **median
of individually fitted gains**; both quantities are saved. Only QC-accepted
events enter the figure, so the overlapping Santa Cruz packet is excluded.

Display spectra span 0.005–0.1 Hz and are saved as `paper_*` fields in the event
NPZ. Fit/phase/acceptance calculations remain on the explicitly selected band.
The prepared phase-velocity table spans 0.015–0.085 Hz, so purple theory is drawn
only there, also excluding the first acoustic-resonance margin. Its shorter
range is deliberate: no unsupported phase-velocity extrapolation is fabricated
to fill the entire axis. The display outside the fit band is not a calibration
claim. Event segmentation, quiet training and response handling follow the
[current pilot](rayleigh_calibration_review.md), not every historical code choice.
If a response is unusable outside the fit band, the plot retains that event's
fit-band curve and records the display issue; it does not invalidate fit-band QC.

## Magnitude fallback

> This section describes `--event-sampling bounded`, the policy of the saved RR52
> pilot. The current default, `all`, evaluates every eligible event at all three
> magnitude levels and does not stop at `--min-events`; `--max-events` is then
> ignored.

The search starts with **Mw>7**, then expands to **Mw≥6.5**, then **Mw≥6** while
fewer than `--min-events` QC candidates have been found. The default target is
**three**. Catalogue count alone does not satisfy it. The same depth, response,
sample-grid, phase, coherence, coverage and overlap checks apply at every level.
The catalogue must contain the smaller events; the runner does not download a
new catalogue implicitly. The supplied M≥5.5 snapshot is sufficient.

`--max-events 8` now means at most eight new attempts **per magnitude level**
(24 maximum overall). Attempts stop once the candidate target is reached.
Already attempted events are not repeated. A 24-hour separation from previously
attempted mainshocks also applies across levels. Each level's counts and stopping
reason are recorded as `magnitude_search` in `comparison.json`. Exhausting a
budget without enough candidates leaves station status unresolved; it does not
relax the QC criteria. Increase the explicit budget to inspect more events.

Use `--no-magnitude-fallback` to restrict searching to Mw>7. An explicit
`--event-id` processes only that event (including Mw≥6 events) and disables the
fallback search. `--min-events` may be increased, but cannot be less than three.
These are pilot rules; candidate count does not replace window/model validation.

## RR52 comparisons

| Fit band | QC candidates | Median of event gains | Gain fitted to median ratio | Search outcome |
|---|---:|---:|---:|---|
| Existing 0.02–0.06 Hz default | 2 | No station estimate | 0.67386, diagnostic only | Eight attempts at each of the three magnitude levels; target not reached |
| Original figure's 0.03–0.07 Hz | 3 | 0.678998 | 0.674374 | Target reached after four Mw>7 attempts; lower levels not needed |

The paper-band candidates are Alaska (`ak0138esnzr`), the February 8 Santa Cruz
event (`usc000f40j`) and Papua (`usb000g7x7`). The February 6 Santa Cruz example
is still excluded for a predicted overlapping packet. Event-gain sample scatter
is 0.004715; this is not total uncertainty or a confidence interval.
The band comparison was made because the original figure marks 0.03–0.07 Hz;
it was not obtained by scanning bands for a desired gain. Broader window,
mode-isolation, depth/reference and uncertainty review remains necessary.

[Paper-band plot](figures/rr52/rayleigh_paper_figure5.png) ·
[Default-band plot](figures/rr52/rayleigh_paper_figure5_default_band.png) ·
[Validation and exact local run paths](rr52_paper_validation_summary.json).
Previous results, masters and source data remain preserved.

## Run and view

To match the original figure's fit band explicitly:

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python -m bruitfm.rayleigh_calibration \
  --catalog results/RAYLEIGH_REVIEW/catalog_usgs_m5p5_2012_2013.json \
  --dispersion-file src/bruitfm/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --network YV --station RR52 \
  --frequency-min 0.03 --frequency-max 0.07 \
  --min-events 3 --max-events 8 --show
```

Close each event plot to continue; the combined Figure 5 comparison appears at
the end. Omitting the frequency options retains the 0.02–0.06 Hz default.
To view the already generated comparison without recomputing:

```bash
open docs/figures/rr52/rayleigh_paper_figure5.png
open docs/figures/rr52/rayleigh_paper_figure5_default_band.png
```

## Rayleigh-wave package

The fundamental-mode dispersion table was generated with **disba 0.7.0**:
`PhaseDispersion(...)(periods, mode=0, wave="rayleigh")` for the water-column
reference and `GroupDispersion(...)` for the packet timing. BRUIT-FM reads that
prepared table; it does not rerun disba or require it in the runtime environment.
[Official disba documentation](https://keurfonluu.github.io/disba/).
The observed waveforms are actual OBS data, not synthetic Rayleigh waves.
ObsPy TauP predicts P/S arrivals; TiSKitPy handles the quiet-data tilt rotation;
Matplotlib draws the figure. NumPy implements the gain grid search with an
independent analytical least-squares check.
