# Pressure atlas: global seafloor pressure-noise models (`compy_pressure.atlas`)

`compy-pressure-atlas` turns an archive of seafloor pressure gauges (differential pressure gauges, absolute
pressure gauges and hydrophones) into station-balanced noise models from 0.02 to 1260 s, stratified by water
depth, distance to the coast and to the shelf break, sensor type and brand, basin, latitude and season. It
replaces the analysis scripts of the November 2026 atlas (branch `analysis/pressure-atlas`, kept as an archive)
with one tested package.

## Three commands

```bash
# 1. native-rate spectra of every channel (one sequential reader of the archive, N workers; resumable)
compy-pressure-atlas native plan.csv --out native/ --cache /fast/disk/sds_cache --sds /Volumes/DATA/SDS --workers 12

# 2. distance of every gauge to the coastline and to the shelf break (Natural Earth 1:10 m, needs cartopy)
compy-pressure-atlas distances channels.csv --out channels_with_distances.csv

# 3. selection, screens, models and figures
compy-pressure-atlas build channels_with_distances.csv --out atlas/ --config atlas.json --workers 10
```

`plan.csv` has one row per channel: `seed_id`, `fs`, `xml` (StationXML with the response) and optionally `files`
(day files separated by `;`; otherwise `<sds>/<year>/<NET>/<STA>/<CHA>.D/*`). The channel table of `build` is
described in `compy_pressure/atlas/pipeline.py`; its required columns are the position, depth, sampling rate,
sensor type and brand, and the path of the native-rate file; optional columns carry the gain correction, the
usable band, upstream exclusions with their reason, the P-Z calibration level and the distances. A missing
optional column switches its step off. Every threshold is in `compy_pressure/atlas/config.py`; a JSON or YAML
file overrides any of them.

## One estimator from 0.02 to 1260 s (`native`)

Every channel is read at its own sampling rate. Windows of 3 h with 50 % overlap start on a fixed 1.5-h grid; a
window with a gap is not used. The PSD of a window is the Welch average of 2700-s Hann segments with 75 % overlap
(the ObsPy PPSD convention for that window), divided by the squared modulus of the full StationXML response and
averaged as power over 1/4-octave bands on the 1/8-octave grid of the atlas (nodes 1060 / 2^(j/8) s).

Periods shorter than 4 s come from the native samples; periods of 4 s and longer from the same window decimated to
1 Hz (zero-phase polyphase filter, flat to 0.4 Hz) and high-passed at 0.2 mHz over 6 h of context on each side,
which removes the tide before it can leak into the lowest Welch bins of absolute gauges. Above about 470 s a
1/4-octave band is narrower than the Welch resolution; there the level is interpolated between Welch bins.

Checks (`tests/test_pressure_atlas.py`):

- white noise of known level is recovered within 1.5 dB from the Nyquist band to 1260 s at 1 and 50 Hz (and
  within 1 dB at 62.5 and 250 Hz in the checks of the November 2026 run);
- a tide of 10^5 counts on 10-count white noise changes no period above 300 s by more than 1.5 dB;
- the vectorised Welch equals `scipy.signal.welch` to 1e-9.

On real data the product equals the 1.25 Hz product of ComPy stage A within 0.2 dB at 3-200 s (YV.RR01, 6360
windows), continues past 667 s where that product is cut by the stage A pre-filter, and does not suffer the
1-octave smoothing of the ObsPy PPSDs (which filled the 15-50 s notch by up to 7 dB).

**1-Hz channels** (most 7D LDH): the data fall by 34 dB between 2.9 and 2.3 s, a decimation filter that their
StationXML does not describe. They are used from 4 s (`population.one_hz_min_period_s`).

## Selection and screens (`build`)

1. **Channel cascade**, each step with its reason: no product; moored or unverified deployment; logger or unknown
   channel; input units not pressure; upstream exclusions (for example metadata known to be placeholders); fewer
   than 7 days; a second channel of the same gauge (the ?DH with the widest band is kept); a level more than 20 dB
   from the basin median at 1-2 s or 2.5-4 s (gross metadata errors); a P-Z calibration level more than 20 dB from
   the DPGs and APGs of the same depth class, or a network and sensor type whose median fails.
2. **Window screens**: no level, clipped (99 % of 24-bit full scale, or 32-bit where most windows exceed 24 bits),
   dead (30 dB below the channel median), deployment edges (days at either end whose 1-10 s level departs by more
   than 15 dB).
3. **Day screen**: day levels against the same-day median of network peers of the same depth class (20 dB after the
   gauge's usual offset); deep gauges without peers against their own median (30 dB).
4. **Band screen**, gauges at 1000 m or deeper: a band median more than 25 dB from peers removes that band.

## Models

Each gauge counts once. Its histogram of window levels, shifted by its gain, is normalised at each period and
weighted by its usable band (cos² taper over the outer octave). Percentiles of the sum give the group's 1-99 %
curves; a continuous piecewise-linear fit in log period with the fewest segments reaching 1 dB RMS gives the low
(5 %), median and high (95 %) models in the form of Peterson (1993), P = A + B log10 T. Curves are drawn only where
8 gauges count at full weight (3 for thin groups, flagged), and fitted only where a quarter of the group's gauges
count. Groups: all, deeper than 200 m, sensor type, depth class, sensor x depth, brand, basin, distance to the
coast, distance to the shelf break, latitude (each also restricted to gauges deeper than 1 km), season, hemisphere
and season, month by hemisphere, and the paired winter-minus-summer difference of each gauge.

## Results

The archive of November 2026 (975 channels, 35 networks) recomputed at native rate, with the report on what
seafloor pressure gauges record: see the run of 18 November 2026 in the ComPy workspace
(`runs/pressure_2026-11-18`) and the summary in the README.
