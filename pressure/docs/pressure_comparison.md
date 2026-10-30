# Expected pressure-level comparison

**19 October 2026 correction:** see [the source review](pressure_reference_review.md).
The supplied points attributed to Webb remain unverified as published universal
bounds. The owner confirmed they came from Wayne Crawford and requested all
of them on plots: default figures show all five pairs as Wayne-supplied (Webb).
They remain excluded from numerical comparisons. Historical products below
are preserved, not recomputed. Brown now has a separate approximate full-band
transcription; it is an in-water reference, not a seafloor acceptance rule.

`pressure_test.py` reads the saved pressure PPSD NPZ referenced by the master
CSV. It calculates empirical 5th, 50th and 95th percentiles of the individual
stored spectra, without reading raw waveforms or applying the response again.
The default station is RR52. Seismic channels are excluded.

## Reference bounds

- **Webb, Zhang and Crawford:** [Infragravity waves in the deep ocean](https://doi.org/10.1029/90JC02212).
  The discrete bounds below are transcribed from the supplied project instructions.
  Review of the author-uploaded full text did not establish these values as
  published universal bounds. They remain in metadata for traceability, with
  excluded rows in new numerical diagnostics; all supplied points remain visible
  on the plots at the owner's request.
- **Brown, Ceranna, Prior, Mialle and Le Bras:** [The IDC Seismic, Hydroacoustic and Infrasound Global Low and High Noise Models](https://doi.org/10.1007/s00024-012-0573-6).
  Section 4 uses 5th/95th percentile levels; Section 4.2 and Figure 8 provide
  hydroacoustic comparisons. The bundled curves are approximate digitizations
  of the solid lines in [Figure 8](https://link.springer.com/article/10.1007/s00024-012-0573-6/figures/8),
  sampled only where the solid curves were visually identified. They are not official numeric
  model tables. The source attribution, raster hash, axis coordinates and sampled
  pixel positions are stored with the reference files.

| Frequency (Hz) | Lower bound (dB re 1 Pa²/Hz) | Upper bound |
| ---: | ---: | ---: |
| 0.005 | 10 | 50 |
| 0.05 | −50 | 10 |
| 0.1 | −9 | 20 |
| 0.2 | 10 | 34 |
| 1 | −25 | −2 |

Brown's figure uses **dB re 1 µPa²/Hz**. Subtracting **120 dB** converts these
levels to the saved pressure PSD units, dB re 1 Pa²/Hz. The implementation uses
hydroacoustic curves, not atmospheric infrasound or seismic curves. Brown's
in-water observations provide a reference comparison, not universal acceptance
thresholds for every seafloor differential pressure gauge.

## Frequency selection and calculation

The response snapshot beside the NPZ supplies the complete channel response.
By default, the nominal comparison band is the contiguous interval around its
sensitivity reference frequency whose amplitude attenuation is at most 3 dB.
This is an explicit computational choice, **not independent confirmation of
the instrument's usable band**. Supply both sensor frequency bounds to override
it once suitable limits are established.

Only smoothing bins fully inside the selected band contribute to comparisons.
The unverified handout points remain excluded regardless of center matching.
Brown comparisons use full smoothing-bin support inside approximately
0.01139-78.10 Hz, without extending to the figure's 0.01-100 Hz axis endpoints.
The digitized dB curves are
interpolated over log frequency and averaged over the same FFT frequencies
inside each smoothing bin as the saved spectra.

The frequency CSV contains reference bounds, observed percentiles, median
offsets, window occupancy fractions, and `within`, `below`, `above` or
`excluded` labels. Brown results within 2 dB of a boundary are marked for
digitization review. This tolerance is a review aid, not a validated error bound.
Overlapping windows are correlated; occupancy fractions do not represent
independent trials or statistical confidence. Existing three-hour smoothing
and window settings are retained.

## Run and save

From the project checkout with the `bruitfm` environment active:

```bash
python -m bruitfm.pressure_test
```

Or select that environment in VS Code and run `src/bruitfm/pressure_test.py`.
Both routes display and save three figures. For background use:

```bash
conda run -n bruitfm python -m bruitfm.pressure_test --no-show
```

Each pressure-channel output directory receives:

```text
pressure_expected/
  frequency_comparison.csv   Per-frequency measurements and comparisons
  comparison.json           Settings, source hashes and reference provenance
  comparison.png            Percentile/reference plot and nominal response band
  ppsd_with_limits.png       PPSD probability distribution with reference limits
  test_details.png           Median clearance and percentage of windows in range
```

The existing master CSV and Excel receive `pressure_expected_*` fields for the
result paths, comparison fraction, tested-bin count and selected band. The first
test status becomes `comparison_only` when its score is blank. Existing manual
scores/statuses and the other five tests are preserved. Both master files are
backed up before replacement, using the existing updater's preservation logic.
Use `--no-update-master` to generate files without changing the masters.

Current diagnostic views use **orange for Brown** and **magenta for the
Wayne-supplied (Webb) points**. All five supplied frequency locations are drawn
at their exact frequencies even where numerically excluded. Their numerical
clearance/occupancy column remains excluded. Brown's PSD-bin-averaged reference curves
use dashed/solid lines. The color map shows
occurrence probability. Limits are also shown outside the comparison band for
context; excluded frequencies do not contribute to test results.

The details view separates the references into two columns. The upper panels
show median clearance to the nearest bound: positive means inside both bounds,
negative means outside. The lower panels show the percentage of PSD windows
within both bounds at each compared frequency. Marker shapes distinguish
`within`, `below` and `above`; rings mark Brown digitization review points.
The shaded ±2 dB strip in Brown's panel is a visual review guide.

No final calibration score is assigned by default. An optional `--write-score`
requires explicit minimum/maximum sensor frequencies. It stores the equal-weight
fraction of tested median bins inside the bounds as a **descriptive screening
score**, with status `screened`; this is not a calibration acceptance decision.
This option deliberately replaces the first test's existing score/status.
Skipped frequencies do not count as failures. Source-specific coverage is
recorded; the fraction does not describe a continuous full-band comparison.

## Historical RR52 diagnostic result (before the reference correction)

The initial comparison reused all **5,249 pressure PSD windows**. Its nominal
response band is approximately **0.00857–7.395 Hz**. Nineteen bins contributed:
four supplied Webb points and fifteen Brown bins. Sixteen median levels lie
inside the bounds and three lie above them, giving an in-range fraction of
approximately **0.8421**. Four Brown comparisons lie within the 2 dB review
tolerance of a boundary.

The supplied 0.005 Hz point is excluded because its full smoothing bin is
below the nominal response band. Missing reference coverage and excluded bins
remain visible in the frequency CSV. The result is `comparison_only`, and
RR52's calibration score remains blank.

The [frequency results](rr52_pressure_comparison.csv),
[settings and summary](rr52_pressure_comparison.json), and
[comparison figure](figures/rr52/BDH_pressure_comparison.png) preserve this initial run.

The literature review did not establish the supplied numbers as verified Webb
bounds. Do not reuse the historical in-range fraction as current reference
validation. Sensor bandwidth/in-situ response and applicability remain open;
no acceptance policy was adopted. Use the new saved-spectrum comparison command
in [the source review](pressure_reference_review.md).
