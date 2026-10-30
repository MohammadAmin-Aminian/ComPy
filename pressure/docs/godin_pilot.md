# Screened full-window pressure pilot

`godin_pilot.py` calculates a separate pressure spectrum directly from raw BDH
counts, following the preprocessing described by
[Godin et al., *Power spectra of infragravity waves in a deep ocean*](https://doi.org/10.1002/grl.50418),
Section 2, paragraphs 8–9. The existing three-hour PPSDs and earlier diagnostic
remain available.

## Calculation

The default selects a seven-day interval near the middle of waveform coverage,
avoiding deployment startup data. An explicit `--start` selects another UTC
interval, and `--days` changes its length. Each 24-hour interval contributes
three overlapping **12.42-hour (44,712 s)** windows: beginning, middle and end.
Their starts are separated by 5.79 hours within an interval; the next day's
first window starts 12.42 hours after the previous day's last start.

Each complete window receives one full-length periodic Hann taper and FFT.
The segment mean is removed before tapering. This detrending rule and the
periodic Hann convention are explicit implementation choices; the source
does not specify them. There are no internal Welch subsegments, FFT padding,
gap interpolation or octave smoothing in this primary product. Gapped,
nonfinite and entirely constant windows are omitted and recorded.

At RR52's 62.5 Hz sample rate, the FFT contains **2,794,500 samples**. Frequency
spacing is approximately **0.0000223654 Hz (0.0223654 mHz)**. Hann equivalent
noise bandwidth is approximately **0.0000335480 Hz**. Frequency spacing does
not imply independent bins, and each spectrum averages 12.42 hours of data.

Screening uses mean **linear counts²/Hz in 0.05–0.06 Hz**, before response
correction. Windows strictly above the 95th percentile are excluded. Threshold
ties remain retained; the excluded fraction need not equal exactly 5% in a
small pilot. This is a noise-based screen, not an earthquake identification.

The full pressure response is evaluated at the exact FFT frequencies and
its squared amplitude divides the raw-count PSD once. Retained spectra are
averaged in **linear power**, then converted to dB for display. Spectra from
positive frequencies through 0.06 Hz are saved, including the screening notch.
Equation 6 is evaluated on this grid, using the existing depth-dependent model.
Only frequencies inside both the selected model and nominal sensor bands
contribute to the residual summary. The sensor-band heuristic and the
inventory-derived depth remain unverified.

Processing is sequential for this single-channel pilot, using bounded reads.
It does not alter the existing four-worker PPSD workflow. Source-file sizes and
modification times are checked before and after waveform processing. The
selected response is checked numerically against the saved PPSD snapshot.

## Run and view

From the project directory, calculate, save and display the pilot:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.godin_pilot
```

To open the already saved three plots without computing again:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.godin_pilot --plot-only
```

For background processing add `--no-show`. Use `--no-update-master` to leave
the master files unchanged. Use `--help` for interval, depth and sensor-band
options. In VS Code, select the `bruitfm` interpreter and run
`src/bruitfm/godin_pilot.py` to calculate and display the pilot.

## Results and figures

The RR52 middle-of-recording week supplied **21 complete windows**. Screening
retained **20** and rejected **one**, an actual rejection fraction of 4.76%.
No requested windows were missing, gapped or entirely constant. These checks
do not establish that the instrument was correctly calibrated.

The retained linear mean is above the model. Its median difference over
**511 eligible FFT bins** is approximately **+17.79 dB**, using the 2,918 m
inventory-derived depth and nominal sensor band. The frequency points are
correlated and do not represent 511 independent tests. This statistic differs
from the earlier +16.13 dB diagnostic in both interval and estimator; they
should not be interpreted as a before/after calibration change.

Three plots distinguish the calculations:

- [Pressure mean and model](figures/rr52/BDH_godin_pilot_pressure_model.png):
  unscreened and retained linear means, retained-window percentile spread,
  the model and residual curve. Gray frequencies are outside the selected
  sensor band. Percentiles describe window spread, not confidence intervals.
- [Screening decisions](figures/rr52/BDH_godin_pilot_screening.png): window
  notch powers over time and their histogram, with threshold and rejections.
- [Same-interval PPSD comparison](figures/rr52/BDH_godin_pilot_ppsd_comparison.png):
  retained Hann spectra and **110 existing PPSD windows** from the same week,
  averaged in matching octave bins. Within each bin, dB values are averaged
  first; their linear powers are averaged across windows. This geometric
  frequency average is a separate diagnostic, not the primary raw FFT mean.
  The two fully eligible bins are approximately **1.79 and 2.32 dB lower**
  than the same-week PPSD mean. Window duration, taper, Welch averaging and
  screening all change together, so this comparison cannot isolate one cause.

The regional model, depth and instrument band still need validation. A short
pilot cannot substitute for the original yearlong average or assign a
calibration correction. **The calibration score remains blank.**

## Saved files and master table

The pressure epoch directory receives a separate folder:

```text
godin_pilot/<interval-start>__<days>days/
  spectra.npz
  windows.csv
  frequency_comparison.csv
  ppsd_comparison.csv
  comparison.json
  pressure_model.png
  screening.png
  ppsd_comparison.png
```

NPZ stores raw-count and response-corrected window PSDs, FFT frequencies,
response powers, window starts, retained flags and both linear means. CSVs
store window decisions, model residuals and matched PPSD comparisons. JSON
records the response epoch, interval, FFT configuration, screening rule,
depth/band assumptions, source metadata and hashes, limitations and paths.
Rerunning the same interval replaces its pilot files. Different intervals
receive different folders; preserve a copy before comparing alternative
settings for the same interval.

The matching master CSV/Excel row receives **`infragravity_pilot_*`** fields
with the interval, window duration, retained/rejected counts, eligible-bin
count, residual and file paths. Both master files are backed up. Existing
fields, scientific scores, earlier diagnostics and extra sheets are preserved.
For RR52, inspect row **13834** on **Master inventory**, in the newly appended
pilot columns. `infragravity_pilot_status` is `comparison_only`.

Curated numerical results and figures are tracked in GitHub. Full spectra and
master files remain local under the existing ignored output directories.

[Settings and summary](rr52_godin_pilot_summary.json) ·
[Window decisions](rr52_godin_pilot_windows.csv) ·
[FFT frequency results](rr52_godin_pilot_frequency.csv) ·
[Matched PPSD results](rr52_godin_pilot_ppsd_comparison.csv).

## Verification

The full suite now passes **47 tests**. Pilot checks independently verify the full
FFT and one-sided normalization for even/odd lengths, Parseval consistency,
window scheduling and sample alignment, screening ties, gap/flat rejection,
response correction and linear averaging. An integration test reads real
MiniSEED files and checks master/extra-sheet preservation and plot-only access.
The RR52 run preserved every existing CSV value and Excel cell, the previous
PPSD and response hashes, and source waveform sizes/mtimes. An independent
NumPy FFT reproduced the first real full-Hann window's saved count spectrum;
response correction and the retained linear mean were checked independently.
The three saved plots were visually reviewed.
