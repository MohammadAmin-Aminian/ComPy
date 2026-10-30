# Depth-dependent infragravity comparison

`infragravity_test.py` implements the third scientific comparison using saved
pressure PPSDs. It requires no nearby stations or waveform reprocessing.

## Model

Reference: Godin, Zabotin, Sheehan, Yang and Collins,
[Power spectra of infragravity waves in a deep ocean](https://doi.org/10.1002/grl.50418),
Equation 6. Its normalization was checked in the
[author-uploaded full text](https://www.researchgate.net/publication/258787403_Power_spectra_of_infragravity_waves_in_a_deep_ocean).

For water depth H, the code solves

```text
x tanh(x) = (2π f)² H / g,   where x = kH
S_p(f,H) = 2000 tanh(x)^(-3/2) cosh(x)^(-2) Pa²/Hz
```

Frequency f is in Hz, depth H in metres, and g defaults to 9.81 m/s².

The figures also draw Eq. 6 without the tanh(kH)^(−3/2) factor,
S_p = 2000 cosh(kH)^(−2), with kH from the same full dispersion relation
(dotted purple). The two agree once tanh(kH) ≈ 1 (above about 10 mHz at
3 km depth) and differ most in the lowest infragravity band; at RR52
(2918 m) the no-tanh form is about 14 dB lower at 1 mHz. Residuals and
summary numbers use the published Eq. 6 only.
Angular frequency enters the dispersion equation; the published coefficient
already has per-Hz units. The model is evaluated in logarithmic form to avoid
overflow. No response correction is reapplied to saved pressure spectra.

## Comparison choices

The default model band is **0.001–0.02 Hz**, a conservative workflow choice.
The sensor band reuses the response-based estimate from the first pressure
comparison. Explicit sensor bounds and water depth can override the defaults.
These inputs still need independent verification.

Only complete smoothing bins contained in both bands are compared. The model's
dB values are averaged over the same FFT frequencies as each saved PPSD bin.
Empirical 5th/50th/95th percentiles describe the observed spectra; the residual
is observed median minus the smoothed prediction. Positive residuals mean
higher observed power density. The `median_minus_model_db` field is populated
only for eligible comparison bins. Separate `diagnostic_*_minus_model_db`
fields retain the median and percentile differences wherever the model can
be averaged, including values outside the selected sensor band. Those
diagnostic values do not contribute to the comparison summary or count as
failed tests. An empty overlap is recorded as `no_comparable_bins`.

The model is a reference prediction, not a lower/upper acceptance envelope.
No tolerance or pass/fail criterion is invented. The original study used
12.42-hour Hann-window spectra with transient screening; this diagnostic
retains the existing three-hour PPSD product. Geographic applicability,
depth, response validity and preprocessing differences need review before
calibration interpretation. Overlapping windows and smoothing are correlated.

The separate [raw-data Godin pilot](godin_pilot.md) now uses full 12.42-hour
Hann spectra and count-based transient screening in a middle-of-recording
week. It saves its results alongside this original diagnostic.

## Run

From the project checkout:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.infragravity_test
```

In VS Code, select the `bruitfm` interpreter and run
`src/bruitfm/infragravity_test.py`. The default station is RR52. Both routes
display and save three figures:

- **Model comparison:** observed percentile spectra in blue and the Godin
  model in orange, including the unsmoothed and PSD-bin-averaged predictions.
- **PPSD with model:** the probability distribution with the orange prediction
  overlaid. Gray frequencies are outside the selected sensor band.
- **Model residuals (Figure 3):** an overview of diagnostic differences and a
  detailed view of the eligible bins. Frequencies are labeled in mHz. The
  horizontal gray segments show the full smoothing-bin width, and vertical
  blue segments show the 5th–95th window spread, not confidence intervals.
  Only blue points contribute to the comparison summary; the two RR52 bins
  overlap strongly and are not independent measurements.

For background use, add `--no-show`. Use `--no-update-master` to generate files
without updating CSV/Excel. Use `--help` for water-depth, sensor-band and
model-band settings. Changing settings replaces this diagnostic's files;
copy its folder first if you want to retain an alternative comparison.

## Files and master results

Each pressure-channel epoch directory receives:

```text
infragravity/
  frequency_comparison.csv
  comparison.json
  model_comparison.png
  ppsd_with_model.png
  model_residuals.png
```

The CSV includes predicted levels, observed percentiles, residuals and reasons
for excluding bins. JSON records depth/band choices, model parameters, source
hashes, counts, limitations and output paths.

The existing master CSV/Excel receive `infragravity_*` fields containing result
paths, depth source, comparison-band bounds, compared-bin count and residual
summaries. The predicted-infragravity test's status becomes `comparison_only`
when its score is blank. **Its score is not assigned.** Existing scores,
other scientific tests, PPSD results, metadata and other Excel sheets are
preserved. Both master files are backed up before replacement.

## RR52 result

Used all **5,249 saved pressure spectra** and the **2,918 m inventory-derived
depth estimate**. The nominal sensor band starts near **0.00857 Hz**. Only
**two** full smoothing bins, centered near **0.01283 and 0.01399 Hz**, fit both
bands. Their median residual is approximately **+16.13 dB**. This limited
comparison cannot establish calibration acceptance or justify a correction.

The [frequency results](rr52_infragravity_comparison.csv),
[settings and summary](rr52_infragravity_comparison.json) and
[figure index](figures/rr52/README.md) preserve the computed diagnostic.

## Recheck of the calculation and Figure 3

Rechecked the equation normalization and units against the source and used an
independent scalar root solver with ObsPy's exact FFT-period masks. The two
retained differences are **+15.8252 dB** and **+16.4279 dB**; their median is
**+16.1266 dB**. The maximum model difference from the independent solver is
below 2×10⁻¹⁴ dB. An independent raw-count PSD calculation for the first saved
window reproduced its stored spectrum exactly at float32 precision.

The original Figure 3 omitted the excluded diagnostic curve and did not show
the large, overlapping frequency-bin widths. Its vertical bars could also be
mistaken for estimation uncertainty. The revised view makes all three
distinctions explicit. The numerical comparison is unchanged; depth,
instrument operating band, geographic applicability and the difference from
the original preprocessing method remain unresolved scientific questions.
