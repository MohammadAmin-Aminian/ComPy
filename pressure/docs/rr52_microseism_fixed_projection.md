# RR52 orientation and separate-training projection review

## Decision

The requested bounded follow-up is complete. A published RR52 orientation
reference has been recovered, and historically quiet-trained weights have been
tested on temporally disjoint saved June data. **Do not adopt this projection
as an automatic microseism correction.** Lower vertical power alone is not
validation: gain-band coherence does not consistently improve.

Keep the original frozen gains, no-tilt products and all blank scientific scores.
The [89-row epoch-bound staging export](microseism_table_staging.csv) remains
unchanged. Its conditional/unresolved statuses are the appropriate software
outcome at this stage; an accepted absolute calibration is not established.

## Orientation source review

The local XML identifies Trillium240OBS sensors at RR52, with both horizontal
dips and azimuths set to zero, no station comments and no sensor serial number.
No deployment/orientation log was found in the scoped local RR52/metadata
folders. This was a targeted search, not an archive-wide records audit.

[Scholz et al. (2017)](https://doi.org/10.1093/gji/ggw426), Table 1, printed
page 1284, reports:

| RR52 estimate | Measurements | Circular-mean orientation | Reported 2-standard-deviation spread |
|---|---:|---:|---:|
| P-wave polarization | 6 | 29 degrees | 10 degrees |
| Rayleigh polarization | 27 | 29 degrees | 17 degrees |

The paper defines this angle clockwise from north to BH1 in its standardized
frame. Its supporting-information description flags different native component
handedness for INSU-IPGP and DEPAS/GEOMAR. The exact conversion to this local
archive remains unverified. Therefore **no 29/119-degree azimuth pair was
inserted into XML**, and no geographic rotation was applied. Horizontal yaw is
also not a measurement of instrument tilt or absolute vertical sensitivity.

The relevant methods/table pages were extracted and the full table was visually
checked using the PDF-reading workflow. The [reference record](rr52_orientation_reference.json)
contains the source URL, local preserved PDF path and SHA-256. The supplement's
detailed mapping was not obtained; its description is not a substitute for
verifying that mapping. A geographic azimuth is unnecessary for the relative-
basis operator test below because training and testing retain identical channel
IDs, responses and native signs.

## Fixed selection and separation

No new fit, waveform read or catalogue query was performed. The saved
[TiSKitPy validation summary](rr52_tiskitpy_validation_summary.json) reproduces
the original `run_f0fa112ac43c/comparison.json` event records. The source
uses physical velocity traces, upward Z, CleanRotator 2.3.1, a catalogue-screened
two-hour interval and its documented low-frequency training bands. A constant
real projection commutes with conversion from velocity to displacement.

Rule `saved-quiet-projection-v1` selects the **latest completed quiet interval
before the earliest test sample**, without using target amplitudes, fitted gains
or earthquake-event QC. This selects 18 April 2013, 03:05:52.650–05:05:52.650 UTC,
about 53/57 days before the 10/14 June tests. The five prior available fits are
reported, with the other four retained only as sensitivities. The July fit is
excluded from the primary/prior set. No best-performing alternative is chosen.

The primary saved H1/H2/upward-Z coefficients are:

```text
[-0.002037953145410086, -0.006681646428432743, 0.9999756008145222]
```

The earthquake associated with that fit was `rejected_qc`; its label remains
unchanged. Event QC is not evidence for or against the separate quiet fit, and
the historical catalogue quiet screen is not proof of an uncontaminated quiet
interval. This separation and the long training-test gap are explicit limitations.
The June intervals were already inspected in earlier diagnostics, so this is
**not preregistered blind validation** or an independent accuracy certificate.

## Result on saved test spectra

The 46 retained hourly four-channel coefficients and original frozen band edges
are reused without raw processing. Exact response/inventory bindings and
pressure/vertical channel epochs match the training source. Coherence and
transfer are computed only after averaging complex cross spectra.

| Test day | Band (Hz) | Projected/original Z power | Transfer-amplitude change | Original coherence squared | Projected coherence squared |
|---|---|---:|---:|---:|---:|
| 10 June | 0.03469-0.04688 | 0.370 | -50.13% | 0.063 | 0.035 |
| 10 June | 0.04688-0.06334 | 0.778 | -19.63% | 0.202 | 0.188 |
| 10 June | 0.06334-0.08559 | 1.093 | -1.48% | 0.481 | 0.435 |
| 14 June | 0.03469-0.04688 | 0.050 | -77.42% | 0.114 | 0.091 |
| 14 June | 0.04688-0.06334 | 0.830 | -1.15% | 0.388 | 0.489 |
| 14 June | 0.06334-0.08559 | 1.080 | +0.92% | 0.865 | 0.798 |

Entries are medians of per-bin ratios/products; a ratio of band medians can
differ. None is a revised pressure multiplier. The lowest band has zero bins
with nominal coherence squared >=0.5 on either day; its large amplitude changes
must not be interpreted as a calibrated correction. In the highest band the
primary projection increases vertical power by about 8-9% and lowers coherence
on both days. Therefore variance reduction at low frequency does not establish
a generally beneficial target-band correction.

Across all five prior fits, highest-band transfer ratios span 0.985-0.990 on
10 June and 1.007-1.012 on 14 June. These are sensitivity ranges, not confidence
intervals or independent estimates. Changing the training choice does not turn
this into validation of an absolute pressure gain.

## Signal-distortion checks

An illustrative pure Rayleigh wave with quadrature horizontal motion was
tested over all relative azimuths and H/Z amplitude ratios 0.5, 1, 2 and 5.
For the primary weights, the maximum amplitude perturbation on that grid is
0.059%, and maximum phase change is 2.01 degrees. This is a conditional toy
model, not an empirical field-distortion bound. Tests also show that an in-phase
physical horizontal component can produce first-order signal loss, invalidating
any general claim of signal preservation from the idealized example alone.

No pressure regressor, frequency-dependent subtraction, new gain fit or response
repair is used. The operator is applied only to temporary Fourier copies;
no corrected waveform or replacement input spectrum is saved. This follows the
existing public-`apply` extraction of saved CleanRotator weights, not its
unvalidated `tfs()` interface; see [TiSKitPy documentation](https://tiskitpy.readthedocs.io/latest/classes/clean_rotator.html).

## Continuation boundary

The microseism research workflow now includes the conditional estimator, frozen
week/band rule, exact-band phase audit, horizontal-coupling pilot, separate-
training sensitivity test and an epoch-bound table export. The result of this
last test is **not to promote an automatic correction**. There is no need to
rerun completed calibration batches or keep trying rotations until they match
Rayleigh. A broader RR40/7D correction campaign is not justified by this test.

Proceed with the separate local-P implementation/review and six-check table
integration. Keep microseism results conditional/unresolved unless new independent
response/model evidence and an agreed scoring policy support acceptance. This
completes the current bounded validation task, not absolute calibration of every
station. A future physical tilt/geometry study should have its own explicit
protocol and evidence, rather than retrospectively tuning this frozen rule.

## Provenance and reproduction

[Machine-readable results](rr52_microseism_fixed_projection.json) contain all
30 band/day/prior-fit comparisons and the full local review's path/hash.
The read-only review verifies 2,545 protected files, including existing staging
results, spectra, responses, masters and historical raw sources. Zero protected
hashes changed. Original calibration gains, scores and response XMLs remain
unchanged.

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python -m scripts.review_microseism_fixed_projection \
  --pilot docs/rr52_microseism_horizontal_pilot.json \
  --quiet-summary docs/rr52_tiskitpy_validation_summary.json
```

This reuses saved products and writes a fresh diagnostic folder. Do not rerun
the raw horizontal pilot. Original snapshots and the completed numerical review
remain sufficient for continuation.
