# Pressure-noise reference review (19 October 2026)

## Update, 26 October 2026: limits now in use

The project instructions give the expected-level limits with their sources:
Brown et al. (2014, PAG) at and above 0.2 Hz, and Webb (1998, Rev. Geophys.)
Figure 2 (non-Arctic) at and below 0.2 Hz. The earlier five "Webb" points,
attributed to Webb et al. (1991) and kept out of the comparison as unverified,
are replaced (kept in `pressure_references.json` under `superseded_20261026`).

| Frequency | Limits (dB re 1 Pa²/Hz) | Source | Check |
|---|---|---|---|
| 100 Hz | −56 to −38 | Brown 2014 | digitized curve at 78 Hz: −55.8 to −37.4 |
| 10 Hz | −49 to −25 | Brown 2014 | digitized: −47.8 to −25.2 |
| 1 Hz | −25 to −2 | Brown 2014 | digitized: −26.3 to −2.8 |
| 0.2 Hz | 10 to 50 | Webb 1998 Fig. 2 | Atlantic ~7×10², Pacific ~4×10⁴ Pa²/Hz |
| 0.1 Hz | −20 to 20 | Webb 1998 Fig. 2 | Atlantic ~1, Pacific ~10–100 Pa²/Hz |
| 0.05 Hz | −25 to 5 | Webb 1998 Fig. 2 | ~0.3–0.5 Pa²/Hz in the noise notch |
| 0.005 Hz | 20 to 50 | Webb 1998 Fig. 2 | Atlantic ~3×10³, Pacific ~6×10⁴ Pa²/Hz |

The Brown values agree with the existing digitization within its 2 dB review
tolerance, so the Brown curve is unchanged but compared only at and above
0.2 Hz. The Webb points are compared at their four frequencies, never
interpolated. The 0.2 Hz bin is judged against Webb only. Plots draw Webb in
red and Brown in green.

## Scope and access

This is a focused review of the pressure-spectrum sections, figures and relevant
equations in eight core references below, not a claim to have read every paper
on ocean noise or independently validated every equation in these papers.
Author-uploaded full text was available for the three Webb papers; publisher
full text for Brown, Godin and Aucan; a university-repository published PDF for
Janiszewski; and the publisher PDF for Stahler. Local direct download of the
Janiszewski PDF returned an empty file; its accessible indexed repository PDF
was used instead, not a purported valid local copy. Some scanned equations in
Webb's papers have imperfect OCR: none of the unverified handout numbers was
promoted based on OCR.

## Sources and what they support

1. Webb, Zhang & Crawford (1991), *Infragravity waves in the deep ocean*,
   [DOI](https://doi.org/10.1029/90JC02212),
   [author-uploaded full text](https://www.researchgate.net/publication/229061541_Infragravity_waves_in_the_deep_ocean).
   Focus: deep-ocean pressure spectra, Figs 1-2 and hydrodynamic filtering.
   Infragravity fluctuations dominate below roughly 0.03 Hz at the sites
   discussed; a 0.03-0.1 Hz noise notch and broad 0.1-5 Hz pressure microseisms
   are described. The reported 0.0005-32 Hz gauge capability is specific to
   their DPG design, not a passband for every pressure sensor. The figures are
   site/current-dependent measured spectra, not a universal low/high table.

2. Webb (1998), *Broadband seismology and noise under the ocean*,
   [DOI](https://doi.org/10.1029/97RG02287),
   [author-uploaded full text](https://www.researchgate.net/publication/251426459_Broadband_seismology_and_noise_under_the_ocean).
   Focus: pressure observations/Figs 2-3 and infragravity climatology.
   Pressure spectra differ among basins, depths and instrument types. DPG and
   conventional hydrophone self-noise differ with frequency. These examples
   cannot certify the five handout frequency locations as global bounds.

3. Webb & Crawford (2010), *Shallow-Water Broadband OBS Seismology*,
   [DOI](https://doi.org/10.1785/0120090203),
   [author-uploaded full text](https://www.researchgate.net/publication/234115412_Shallow-Water_Broadband_OBS_Seismology).
   Focus: pressure/depth comparisons in Figs 2 and 4 and associated models.
   Shallower water admits higher-frequency wave loading. Its illustrative
   near-surface spectrum and depth-transformed spectra are not low/high
   acceptance envelopes or an independent pressure gain calibration.

4. Brown et al. (2014; online 2012), *The IDC Seismic, Hydroacoustic and
   Infrasound Global Low and High Noise Models*,
   [open PDF](https://link.springer.com/content/pdf/10.1007/s00024-012-0573-6.pdf).
   Focus: sections 2-4, especially 4.2/Fig. 8, printed p. 371. Low/high sensor
   curves are 5th/95th percentile curves; network envelopes take the minimum
   low and maximum high at each frequency. Hydroacoustic data are in-water IMS
   hydrophones, not seafloor DPGs. Figure 8 has 0.01-100 Hz axes; the solid
   traces do not establish measured numerical support at every axis endpoint.
   Its pressure ordinate is dB re 1 microPa^2/Hz. Subtract 120 dB for Pa^2/Hz.
   Figures 9-10 are atmospheric infrasound and must not be used here.

5. Godin et al. (2013), *Power spectra of infragravity waves in a deep ocean*,
   [full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/grl.50418).
   Focus: observations, theory and Eq. 6 model. The model includes frequency
   and water-depth dependence; its amplitude was derived using MOANA records
   and nominal, not individually calibrated DPG responses. It is not a
   universal upper/lower envelope. The paper's 12.42-hour Hann spectra and
   transient screen differ from our saved 3-hour, octave-smoothed PPSDs.
   No Godin normalization or adopted workflow band was retuned here.

6. Aucan & Ardhuin (2013), *Infragravity waves in the deep ocean: An upward
   revision*, [full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/grl.50321).
   Focus: data preparation, Fig. 1 and spatial/seasonal variability. The
   integrated-wave-height band is 0.00083-0.011 Hz (20-1.5 minutes), not the
   entire sensor band or a universal cutoff. They distinguish low-frequency
   forcing and sheltering/seasonal effects. An integrated surface-wave-height
   measure is not directly a pressure PSD acceptance bound.

7. Janiszewski et al. (2023; online 2022), *Broad-band ocean bottom seismometer
   noise properties*, [published repository PDF](https://escholarship.org/content/qt2zw34048/qt2zw34048_noSplash_a7b96bb4d918fd13f1ac980b22c8562a.pdf),
   [article](https://doi.org/10.1093/gji/ggac450).
   Focus: processing/inclusion, spectra/coherences and limitations. The study
   analyses 551 stations across 18 experiments over 0.001-1 Hz, with strong
   depth/instrument effects. The selected, screened and averaged observations
   are not the distribution of every raw window. Representative microseism
   peaks do not define our pressure sensor's usable band. Their APG treatment
   and units must be reconciled before any numerical pressure-spectrum import.

8. Stahler et al. (2016), *Performance report of the RHUM-RUM ocean bottom
   seismometer network around La Reunion, western Indian Ocean*,
   [open PDF](https://adgeo.copernicus.org/articles/41/43/2016/adgeo-41-43-2016.pdf).
   Focus: instrument descriptions/Fig. 5 and Appendix A2-A3. DEPAS hydrophones
   and INSU DPGs are different hardware; generic nominal corners do not verify
   each in-situ response. Most INSU DPG responses borrowed RR28's calibration,
   with significant uncertainty explicitly acknowledged. RR52's local metadata
   says DPG, 62.5 samples/s; ObsPy's `hydrophone` special handling means pressure
   processing, not verification that the hardware is a conventional hydrophone.

## Corrections implemented

- Preserve all five handout points but mark published-bound attribution
  unverified. At the owner's subsequent request, default figures display
  **all five lower/upper pairs** as "Wayne-supplied (Webb)". Display is separate
  from numerical inclusion; the owner confirmed the instructions came from
  Wayne Crawford. New quantitative diagnostics
  retain excluded rows with reason `unverified_project_instruction_bounds`,
  never count them toward the in-range fraction. No fabricated Webb envelope.
- Keep the original 17-row Brown CSV untouched. A new 37-row full-band CSV
  retains those samples and adds 20 manually inspected samples of the solid
  curves, with all pixel coordinates and the same axis mapping. Actual
  digitized support is approximately **0.01139-78.10 Hz**, not the plot-axis
  limits 0.01-100 Hz. No extrapolation. The approximate 2 dB review flag is
  not a certified digitization uncertainty.
- Compare Brown at centers whose *full smoothing support* lies inside both
  the reference support and the explicitly declared/inferred comparison band.
  No arbitrary 1 Hz numerical restriction remains. This is descriptive, not
  an adopted seafloor pass/fail rule. Plot overlays display reference curves
  for context; they are not certified instrument bands.
- Label approximate environmental regions separately from sensor bandwidth.
  Fixed gains across the full display band and deployment are hypothetical
  extrapolations, not demonstrated broadband/time-invariant calibrations.
  Local-P band-centre shifts remain approximate because the already smoothed
  PSDs span band boundaries; no exact corrected waveform/PSD is claimed.
- Add a gain-only dB plot (20 log10(g)), so near-agreement is visible without
  interpreting closeness to an ambient-noise curve as calibration proof.

## Historical results and reproducibility

Historical CSV/JSON/PNG products retain their old sources, fractions and labels.
For example, RR52's saved 16/19 = 0.8421 included four handout frequencies;
it is a legacy descriptive number, not evidence of verified Webb acceptance.
No saved result, master, raw data, response or gain was rewritten or refitted.
New combined reviews and coverage plans add explicit historical-reference
warnings when reusing old handout-inclusive or restricted-Brown results. They
retain those numbers and do not automatically recompute completed checks.
Do not rerun the completed network/calibration batches after this source update.
New network-code fingerprints may invalidate historical automatic caches;
use this saved-spectrum viewer, not a calibration runner.

Run the latest RR52 display (creates a new timestamped folder and opens the
four-panel comparison plus the gain-shift plot):

```bash
cd /Users/mohammadamin/Projects/bruit-fm
PYTHONPATH=src:. /opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python \
  scripts/plot_rr40_gain_comparison.py --station RR52 --show
```

It uses all 5,249 saved pressure PSD windows, no raw processing or average PSD
curves. The name is historical; `--station RR40` remains available.

## Verification

- Final full suite: **423 passed, 110 existing/synthetic ObsPy warnings**, 106.68 s.
  An additional 60-test review/coverage regression run passed independently.
- New regression cases verify pressure gain uses 20 log10(g), finite positive
  gains, local-P support/no mutation/overlap refusal, the complete pixel-axis
  mapping and exclusion of unverified points from numerical fractions. The
  later display-restoration tests additionally check all ten supplied markers.
- Real RR52 and RR40 figures were generated from saved PSDs only, and the
  RR52 four-panel plot and gain-only figure were visually inspected.
- Preservation audit of the **353** recorded six-check inputs: **352** retain
  their saved SHA-256 hashes. The only changed entry is the authorized
  `src/bruitfm/quality_review.py` code update adding historical-reference
  warnings. No saved data, response, master, result or score input changed.
- The historical Brown CSV is unchanged. New full-band samples remain an
  approximate manual transcription, not certified numerical reference data.
