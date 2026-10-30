# RR52 bounded horizontal-coupling pilot

## Decision

The P/H1/H2/Z diagnostic is implemented and exercised on two days in RR52's
original base week. It **does not establish that tilt is negligible**, and it
does not produce a cleaned waveform, corrected gain or scientific score.
The lowest gain band has substantial horizontal-associated vertical motion.
The physical source cannot be identified from this projection alone.

Retain all frozen estimates as conditional research results. Do not populate
pass scores or apply gains. The [89-row staging export](microseism_table_staging.csv)
and earlier [phase audit](microseism_gain_band_audit.md) remain unchanged.

## Real-data result

Two nominal days were requested: 10 and flagged 14 June 2013. Each retained
23 hourly windows. The 23:00 hour on each day was refused because pressure
files crossing midnight have a 0.0024-second native-grid offset. No interpolation,
clock adjustment or forced merge was made. Every retained P/Z hourly spectrum
reproduces its saved counterpart within relative tolerance 1e-10.

| Day | Gain band (Hz) | P/Z coherence squared | Z power predictable from H | Partial P/Z coherence squared given H | Residual/original transfer amplitude |
|---|---|---:|---:|---:|---:|
| 10 June | 0.03469–0.04688 | 0.063 | 67.5% | 0.036 | 0.447 |
| 10 June | 0.04688–0.06334 | 0.202 | 37.0% | 0.267 | 0.858 |
| 10 June | 0.06334–0.08559 | 0.481 | 8.4% | 0.503 | 0.926 |
| 14 June | 0.03469–0.04688 | 0.114 | 96.9% | 0.043 | 0.093 |
| 14 June | 0.04688–0.06334 | 0.388 | 38.0% | 0.563 | 0.885 |
| 14 June | 0.06334–0.08559 | 0.865 | 9.4% | 0.881 | 0.914 |

All entries are band medians, not estimates of instrument error or literal
percentages of tilt. Low P/Z coherence in the lowest band makes it especially
unsafe to interpret an amplitude as a calibrated pressure response on these
days. A strong weekly aggregate does not imply every day has the same wavefield.
No flagged days/bands were removed from the original frozen estimate.

These are targeted diagnostics, not independent held-out validation. The
full-day run extends an initial six-hour pilot and is not an independent repeat.
Finite-sample projection bias is important: tests with independent synthetic
channels still show apparent horizontal predictability. Frequency bins are
correlated, so the table is not an uncertainty/confidence calculation.

## Method and limitations

The four channels are nominal pressure BDH and displacement-corrected
BH1/BH2/BHZ, with metadata-based upward Z. Native 62.5 Hz samples are processed
without resampling or gap filling. A periodic Hann taper, segment-mean removal,
full complex responses and one-sided density normalization reproduce the
existing P/Z convention. The complete four-channel cross matrix is averaged
before forming any coherence or transfer. Individual coefficients are saved
to allow further diagnostics without rereading waveforms.

The H block is normalized by channel powers before conditioning checks, avoiding
unit-dependent numerical inversion. Its Schur complement gives conditional
P/Z coherence. The last column is the magnitude of the Z/P transfer after
statistically projecting Z onto H and subtracting that projection, divided by
the original transfer. It is **not** a pressure-multiplier estimate and must
never be applied to data. Ill-conditioned horizontal blocks and vanishing
residual powers are masked rather than assigned plausible values.

A synthetic no-tilt common seismic mode loses over 99% of its transfer in the
same projection. Thus a large change cannot, by itself, be attributed to tilt:
real Rayleigh motion also has coherent horizontal components. This is why no
automatic subtraction is adopted. [Bell, Forsyth & Ruan (2015)](https://doi.org/10.1785/0120140054)
describe horizontal-channel tilt removal and signal-distortion considerations;
the present calculation is a diagnostic, not a reproduction of their correction.

Both RR52 horizontal azimuths are recorded as **0 degrees**, while both dips are
0 degrees. The channels are statistically usable as two observed signals, but
this metadata does not establish two geographic axes or a physical tilt angle.
The ambiguity is preserved, not repaired by assuming BH1/BH2 are north/east.
Absolute vertical response, sensor orientation and clock remain unverified.
The old P/Z constructor used a less strict cross-file grid merge; this pilot's
two midnight refusals flag a processing limitation, not a retrospective change
to or rejection of the original weekly products.

## What remains before acceptance

**Completed follow-up:** the [orientation/separate-training review](rr52_microseism_fixed_projection.md)
now addresses items 1-2 within a bounded diagnostic scope. It records a published
orientation reference without repairing native metadata and tests fixed prior
quiet weights on separate spectra. The outcome is not to adopt an automatic
correction. Retain conditional results and proceed to local-P/status integration;
the original broader acceptance requirements below are not claimed satisfied.

1. Verify horizontal-axis metadata from deployment/instrument records, or adopt
   an explicitly reviewed relative-axis procedure without claiming geographic
   orientation. Do not edit response snapshots from a guess.
2. Design a quiet-trained, signal-preserving tilt diagnostic and test its
   distortion on independent periods. The current in-sample H projection is not
   an acceptable calibration correction. Use saved coefficients first.
3. Check RR40 separately and obtain horizontal data for the 7D controls if
   available; those configured inputs currently contain P/Z only.
4. Retain independent vertical-response/model constraints and agree scientific
   scoring criteria before accepting absolute gains. The 7D offset is unresolved.

Further blind calibration batches would not resolve these missing constraints.
The code can already export conditional/unresolved outcomes for the reviewed
stations, but that is not a completed absolute-calibration validation.

## Provenance and execution

The [curated numeric snapshot](rr52_microseism_horizontal_pilot.json) identifies
the full local `pilot.json` and its SHA-256. Full cross matrices, hourly
coefficients, window decisions, exact response/input hashes and band statistics
are under:

`results/REVIEWS/MICROSEISM_HORIZONTAL/run_20261015T155523Z_d3564085/`

All **2,498 protected-file hashes** match after the run. Original raw data,
responses, masters, scores, calibration results and frozen selections are
unchanged. Preliminary/unfinished development outputs remain separate and are
not used in the curated result.

The completed pilot command was:

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python -m scripts.pilot_microseism_horizontal \
  --gain-audit docs/microseism_gain_band_audit.json \
  --master-csv exports/master_inventory.csv \
  --data-dir /Users/mohammadamin/Desktop/Data_Bruit_FM/Data/YV/RR52 \
  --station YV.RR52 --start 2013-06-10 --start 2013-06-14 --hours 24 \
  --horizontal-ids YV.RR52.00.BH1 YV.RR52.00.BH2
```

Do not rerun it simply to continue analysis; use the saved coefficients.
