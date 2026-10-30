# Compare pressure spectra with nearby stations at similar depth

`src/bruitfm/nearby_test.py` implements the second spectral check in the
master table. It reuses response-corrected pressure PPSD NPZ files and their
saved StationXML snapshots. It reads no raw waveforms and applies no second
response correction. The network runner discovers targets automatically;
references can come from any network already registered with saved pressure
spectra in the same master.

## Starting criteria

These are **proposed screening settings**, configurable for discussion with
the supervisor. They are not published universal calibration thresholds and
are not fitted to the YV station spacing.

| Criterion | Preferred comparison | Exploratory comparison |
| --- | --- | --- |
| WGS84 geodesic separation | At most 50 km | At most 100 km |
| Relative water-depth difference | At most 5% | At most 10% |
| Declared sensor type/model | Matching nonempty type and matching model labels | Different or unknown labels allowed, explicitly flagged |
| Common qualifying UTC days | At least 30 | At least 30 |
| Distinct reference stations supporting a day/frequency | At least 2 | At least 2 |

Depth difference is `abs(H_target - H_reference) / mean(H_target, H_reference)`.
Depth comes from negative inventory elevation; it is an estimate rather than
verified bathymetry. Sensor labels are normalized for case/whitespace. The
channel code does not identify the hardware. In particular, generic YV
`Hydrophone` inventory labels do not independently confirm APG, hydrophone or
DPG hardware. Review these labels before interpreting a preferred comparison.

A reference with a recorded score of zero in any of the six checks is excluded.
Blank scores allow a provisional reference; they do not certify calibration.
Channels at the target's own `network.station` are excluded even if they have
different locations or response epochs. Single-pair spectra can be inspected
as context, but cannot produce an ensemble comparison result.

## Paired estimator and frequency support

1. Reload and validate each pressure NPZ, checking identity, physical units,
   sample rate, finite spectra and unique timestamps. Hash the NPZ and saved
   response for provenance. Convert the saved spectra to daily summaries once
   and reuse those arrays across targets.
2. Use complete PPSD windows contained in each UTC day. Exclude windows crossing
   midnight. Require at least **12 hours of union coverage per station/day**;
   overlapping windows do not double-count coverage. Twelve hours is another
   configurable operational choice. The daily statistic is the median of
   saved **dB** spectra, not a mean of linear power.
3. Pair common qualifying UTC days. This is a daily-resolution comparison:
   within-day window start times and exact shared hours are not synchronized.
   Sparse or strongly different within-day coverage needs review even when
   both days meet the coverage floor. No whole-deployment median is compared
   against a different station's unrelated recording interval.
4. Require matching PPSD window length, overlap and smoothing width. Sample
   rates and period-bin centers may differ. Use only full smoothing bins inside
   both nominal response bands and both actual FFT-supported bands. The response
   band uses the existing 3 dB heuristic and is not an independently verified
   sensor cutoff. Interpolate the already smoothed reference dB spectrum in
   log frequency onto target centers, without extrapolation; this is an
   approximation, not exact raw-FFT rebinning.
5. Choose one epoch/location per distinct `network.station`, ordered by
   preferred tier, most paired days, shortest distance, then stable identifiers.
   Separate epochs from a station are not merged, and multiple channels from
   that station do not count as independent references.
6. First try the preferred reference group. If it has no supported comparison
   bins, use the allowed wider cohort and label the result `exploratory_only`.
   Each daily/frequency reference is the equally weighted station median.
   Require at least two simultaneous reference stations and at least 30
   supported days **per retained frequency**. Disjoint pairs covering different
   months cannot falsely satisfy simultaneous reference coverage.
7. Subtract the same-day peer median from the target daily median. Save the
   median and 5th/95th percentiles of these paired daily differences at each
   frequency. The scalar summary is the median of the retained per-frequency
   median differences. A positive difference means higher target pressure PSD.

Default PPSD windows remain three hours with 90-minute spacing. The comparison
time plot has **one-day** resolution. Window overlap, spectral smoothing and
correlated neighboring stations prevent interpreting the plotted spread as a
confidence interval or the number of days as independent measurements.
Coast distance, seafloor type, weather and depth-dependent wave propagation
are not corrected by this check. A spectral offset alone cannot establish a
gain error, correction factor or failed calibration.

## Outputs and master files

Every median difference is accompanied by its root mean square,
√(mean d²) over the paired days (`difference_rms_db` in the frequency and
pair tables, `rms_difference_db` in the daily table, and a dashed black
curve in the difference plot). Unlike the median, the RMS combines the
systematic offset and the day-to-day scatter, so it is never smaller than
the absolute median offset.

Each target pressure epoch gets a `nearby/` folder beside its NPZ:

| File | Content |
| --- | --- |
| `candidates.csv` | Every available reference epoch, distance, depth difference, label match, common days, selection tier and exclusion reason |
| `frequency_comparison.csv` | Retained frequencies, paired days, observed/reference levels and paired residual percentiles |
| `daily_comparison.csv` | Daily residual summaries, frequency support and target coverage |
| `pair_comparisons.csv` | Per-reference per-frequency paired differences, including single-pair context |
| `comparison.json` | Criteria, source hashes, nominal bands, estimator definitions, counts, status and limitations |
| `spectra_comparison.png` | Matched target/peer spectra and paired residual spread; zero marks equality, not an acceptance limit |
| `daily_residuals.png` | Time versus frequency with red positive and blue negative residuals; missing/unsupported cells remain blank |
| `neighbor_selection.png` | Local separation/azimuth and distance/depth screening; all exclusions remain available in the CSV |

The CSV and Excel masters receive matching `nearby_*` file paths and numerical
diagnostics, including `nearby_reference_stations`, `nearby_joint_common_days`,
`nearby_compared_frequency_bins` and `nearby_median_difference_db`. Network
processing records `nearby_error` separately from scientific status.

The second test's status is updated when its score is blank. Its score remains
blank, and existing manual scores/statuses, unrelated fields, Excel formatting
and extra sheets are preserved. Both masters are backed up before replacement.

| Status | Meaning |
| --- | --- |
| `comparison_only` | Preferred group supports a diagnostic comparison |
| `exploratory_only` | Comparison needs the wider or less certain reference cohort |
| `no_suitable_neighbors` | No candidate provides a usable pair under the selected settings |
| `insufficient_common_days` | Otherwise suitable candidates lack the required common days |
| `insufficient_reference_stations` | Valid pairs exist but fewer than the required distinct stations |
| `insufficient_joint_coverage` | Enough separate pairs exist but not enough jointly supported daily/frequency cells |
| `no_comparable_frequencies` | No full smoothing bins survive the shared band checks |

These are completed diagnostic outcomes, not pass/fail scores. Unreadable target
products or invalid target metadata are processing errors; the network runner
records them and continues other epochs. Reference errors are audited and
excluded. Empty results still save CSV headers, JSON and informative plots.

## What to run

**Wait for the current network PPSD calculation to finish.** Run one master
writer at a time. The nearby stage checks the PPSD workflow summary and refuses
master updates while that batch is unfinished. If a batch was interrupted,
resume it first. This guard is a workflow check, not a general file lock.

Run from Terminal:

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python -m bruitfm.network nearby
```

This processes discovered pressure targets across all available networks.
To process YV targets only, keeping the full registered reference pool:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network nearby --network YV
```

Open one station's three saved figures:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR44 --view-kind nearby
```

For an interactive single-station run using saved NPZ paths directly:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.nearby_test --network YV --station RR44
```

Use `--no-show` to save without opening figures. `--no-update-master` leaves
the CSV/Excel untouched while still saving diagnostic files beside the NPZ;
use an isolated export with copied NPZs for a fully isolated pilot.

Every criterion is configurable on either processing entry point:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network nearby \
  --preferred-distance-km 50 --maximum-distance-km 100 \
  --preferred-depth-fraction 0.05 --maximum-depth-fraction 0.10 \
  --minimum-common-days 30 --minimum-daily-coverage-hours 12 \
  --minimum-reference-stations 2 --response-loss-db 3
```

Network comparisons and plots run sequentially using cached daily arrays; no
additional waveform FFT calculation is performed. Saved-stage reuse depends on
the code, settings and the entire ready pressure-reference pool's metadata,
scores, NPZ hashes and response hashes. New or changed peers therefore invalidate
previous unresolved results. Stage progress is in `_workflow/nearby.json`.

## Reference context and verification

[Doran et al., *Calibration of Differential Pressure Gauges Through in Situ
Testing*](https://doi.org/10.1029/2019EA000783) examined a compact array with
approximately 5 km aperture and independent in-situ/post-deployment calibration.
Its [author PDF](https://igppweb.ucsd.edu/~gabi/papers/ess.dec19.pdf) provides
calibration context; it does **not** prescribe the 50/100 km screening settings
or a broadband neighbor pass/fail tolerance.

[Godin et al., *Power spectra of infragravity waves in a deep
ocean*](https://doi.org/10.1002/grl.50418) motivates attention to water depth and
regional wavefields. It also does **not** provide a universal nearby-station
calibration radius. The present code implements an operational paired spectral
comparison, not either study's instrument-calibration procedure.

Automated checks recover a known 6.0206 dB PSD offset after doubling synthetic
counts, verify geographic/depth masks and log-frequency interpolation, reject
duplicate-station and disjoint-coverage references, exercise cache invalidation,
and preserve CSV/Excel scores and extra sheets.

Real-data verification uses isolated copies of completed RR43, RR44, RR47,
RR48 and RR52 pressure products with the unchanged default criteria. The
[verification summary](nearby_validation_summary.json) records the actual
diagnostic outcomes; it does not claim a whole-network execution or acceptance.
Example plots: [spectra](figures/nearby_validation/RR44_spectra_comparison.png),
[daily differences](figures/nearby_validation/RR44_daily_residuals.png),
[selection](figures/nearby_validation/RR44_neighbor_selection.png).
