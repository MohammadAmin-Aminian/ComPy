# Paired pressure–vertical spectra

`paired_spectra.py` reads synchronized raw waveforms and saves pressure and
upward vertical-displacement spectra, complex cross-spectra, squared coherence,
and a pressure-to-displacement transfer function. This is a shared diagnostic
foundation; it does not assign a calibration score or correct an instrument.

## Estimator and conventions

The default pilot is seven interior days, split into nonoverlapping 3,600-second
segments. Each segment has its mean removed and receives a periodic Hann taper.
There is no padding, octave smoothing, gap filling or implicit resampling.
One-sided FFT densities use `1 / (sample_rate * sum(taper²))`, doubled away from
DC and the even-length Nyquist bin. DC is omitted. The display/save interval is
0.01–0.5 Hz, configurable and explicitly **not a verified sensor band**.

Pressure counts are divided by the full complex pressure response (`DEF`);
vertical counts by the full complex displacement response (`DISP`), once each.
The vertical axis is upward: inventory dip −90° keeps the sign; +90° reverses it.
Both channels must have equal sample rates, contiguous finite nonflat samples
and aligned sample grids within one microsecond. Each retained pair stays inside
the requested interval and shared response epoch. Missing, gapped, conflicting,
flat or unaligned segments are audited and omitted.

For corrected Fourier coefficients P and D, define `S_dp = mean(D * conj(P))`.
Auto-power and **complex** cross-power are averaged in linear units across all
retained segments before computing:

- `H_dp = S_dp / S_pp`, in m/Pa, predicting upward displacement from pressure.
- `coherence² = abs(S_dp)² / (S_pp * S_dd)`.
- Transfer phase is `arg(H_dp)`, in degrees, wrapped to ±180°.

At least two retained segments are required. Segment/bin correlation prevents
interpreting counts as independent degrees of freedom. Coherence has no adopted
pass threshold or confidence interval; high coherence alone is not calibration.
Pressure and displacement PSD units are Pa²/Hz and m²/Hz, respectively. Invalid
response/power bins are saved with an explicit Boolean mask; placeholder zeros
in those bins must not be interpreted as measurements.

[Ruan et al.](https://doi.org/10.1002/2014JB011162) motivates a displacement/pressure
transfer measurement in the 0.1–0.2 Hz microseism band. Their workflow includes
tilt correction. This stage uses uncorrected vertical data and does not reproduce
the full published method, infer sediment properties, test event arrivals or
establish a gain correction. The plot shades that band for inspection only.

## Pair selection and provenance

Pairing uses pressure input units, seismic input units, exact inventory/source
keys and near-vertical orientation (absolute dip within 0.1° of 90°), rather
than channel suffixes alone. Default partners are at the same network, station
and location. An explicit vertical SEED ID can select another location at the
same station. Multiple overlapping candidates/vertical epochs remain
`ambiguous_pair`; no partner is selected silently or epoch stitched. Each
pressure epoch is handled independently by the network stage.

Both channels require saved PPSD response snapshots. These are compared with
live complex responses on the exact saved frequency grid, and orientation is
checked too. Saved PPSD magnitudes are neither loaded nor used for cross-power.
Raw header coverage is clipped to both response epochs. Default interval
selection centers on these bounds; it does not guarantee continuity and does
not search the archive for the most favorable coherence.

JSON records both channel epochs, physical units, orientation, estimator,
source sizes/mtimes, inventory/snapshot/code hashes, retained counts and limits.
Source metadata and inventory/snapshot hashes are checked after processing.
Clock correctness, actual hardware bands and water depth remain unverified.

## Commands and preservation

A single-station pilot defaults to **no master updates**:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.paired_spectra \
  --network YV --station RR52 --start 2013-06-10 --days 7
```

Every standalone calculation creates a new folder under the pressure epoch's
`paired_spectra/run_<unique-id>/`. Existing runs, Godin products and PPSDs survive.
Outputs are `spectra.npz`, `spectra.csv`, `windows.csv`, `comparison.json`, and
`spectra_coherence_transfer.png`. NPZ includes individual auto/cross-spectra,
linear means, complex responses, transfer, coherence and retained start times.

Open a saved plot without another waveform calculation:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.paired_spectra \
  --plot-only /absolute/path/to/paired_spectra/run_ID/comparison.json
```

`--update-master` appends only `paired_*` fields with CSV/Excel backups. It preserves
all six scientific scores/statuses and existing sheets/styles. An unfinished
saved workflow blocks these standalone writes. The guard is a checkpoint check,
not a general process lock: run production updates only when other batches are idle.

The network command discovers paired channels, isolates errors, saves checkpoint
reports and **does update the master**. Start with a selected station:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network paired \
  --network YV --station RR52 --paired-start 2013-06-10 --paired-days 7
conda run --no-capture-output -n bruitfm python -m bruitfm.network view \
  --network YV --station RR52 --view-kind paired
```

Use `--station` selection rather than restricting discovery to the pressure
`--seed-id`, which would exclude its raw partner. Network options also expose
segment length, frequency limits and explicit vertical ID. Cached results depend
on both raw channels, selected settings, station epoch availability, live
inventories, saved response snapshots and code. Interrupted paired batches may
resume; unfinished other workflows block master writes. `all` now includes
`paired`; it should not be used merely to retry RR12 or RR42.

Coverage/pairing outcomes are `insufficient_common_coverage`, `no_vertical_pair`,
`ambiguous_pair`, `no_complete_pairs` and `no_valid_frequencies`. These are saved
unresolved diagnostics. Unexpected exceptions are processing failures. Neither
kind is a failed calibration score.

## RR52 validation

For 2013-06-10 through 2013-06-17 UTC, 167 of 168 hourly pairs were retained.
The last pair's sample alignment would extend beyond the requested interval.
All 1,765 saved bins are finite; 361 are in 0.1–0.2 Hz. Median squared coherence
there is 0.80214, median transfer amplitude is 4.98579e−8 m/Pa, and median wrapped
phase is +4.119°. These are correlated-bin summaries of an uncorrected pilot,
not thresholds, uncertainty estimates or instrument corrections.

Independent SciPy periodograms and `csd` reproduce the first pressure,
displacement and complex cross-spectrum after response correction. Independent
means reproduce the saved transfer/coherence. The PNG was visually reviewed.
Synthetic tests verify even/odd normalization and Parseval consistency, known
complex gain, phase/delay correction, polarity and independent noise. Real
MiniSEED tests cover gaps, flat segments, unaligned grids, pairing ambiguity,
unequal rates, cache invalidation, restricted master writes, extra-sheet
preservation and saved-plot viewing. Production CSV/Excel hashes are unchanged.
Network integration was tested in temporary exports; no full YV paired batch ran.

[Validation summary](rr52_paired_validation_summary.json) ·
[Saved figure](figures/rr52/paired_spectra_coherence_transfer.png).

The subsequent microseism stage uses the current master-linked
`run_bfbef17f2c13`. Its `spectra.npz` is byte-identical to the older validation
run above. A new read-only analyzer verifies the saved aggregate from its 167
window spectra, reports daily stability and displays RR52's fixed Rayleigh
pressure multiplier without fitting it. The result is a non-identifiable
admittance diagnostic, not a calibration or score. See the
[source-supported RR52 diagnostic](rr52_microseism_diagnostic.md).

## Existing Godin pilot errors

RR12's 21 selected 12.42-hour windows are continuous but entirely constant at
−8,268,016 counts. Four additional interior windows sampled in March, April,
June and July are also flat. This explains the default rejection, but does not
establish that the entire archive is flat or identify the hardware/root cause.
An alternate week should not be selected to hide this evidence. Check original
pressure acquisition/export and response metadata before interpreting it.

RR42's pressure archive has one file covering 2013-05-08 23:58:42.300 UTC through
2013-05-10 00:00:01.060 UTC. Seven days cannot fit. Additional source data, or an
explicitly separate shorter-duration experiment, would be needed. The 12.42-hour
scientific window should remain unchanged. Original pilot/master statuses remain
preserved; neither station was rerun or classified.

[Bounded-read investigation](pilot_error_investigation.json).
