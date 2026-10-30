# Station-coloured network PPSD overlays

**Reference correction, 19 October 2026:** at the owner's request, all five
lower/upper pairs from Wayne's instructions are shown as "Wayne-supplied (Webb)".
Their display does not establish published universal acceptance bounds or put
them into numerical screening. Brown uses the approximate full-band solid
curves (about 0.01139-78.10 Hz actual digitized support), labelled as in-water
context, not universal seafloor limits. Old PNGs remain historical.
See [the literature review and RR52 command](pressure_reference_review.md).

## Pressure probability heat map (current requested view)

```bash
cd /Users/mohammadamin/Projects/bruit-fm
PYTHONPATH=src:. conda run --no-capture-output -n bruitfm python \
  -m bruitfm.plot_network_ppsd --network YV --channels BDH --probability --show
```

The default map is ObsPy's classic `pqlx` PPSD palette (white near zero,
blue/green/yellow/red toward higher occurrence). Installed modern ObsPy's
current default is actually viridis; `--cmap obspy-default` selects it exactly.
The default colour scale is 0–30%, matching ObsPy's non-cumulative convention.
For the present pooled YV distribution, `--max-percentage 15` gives stronger
contrast. Optional `--min-percentage 0.1` hides bins below 0.1% **only in the
display**; the title/manifest disclose it, saved counts and probability
normalization are unchanged, and no data selection or renormalization occurs.
That cutoff applies to both detail and full-range PNGs; all counts remain in
the histogram NPZ. Omit it to show all nonzero bins.

The saved YV spectra have one-octave frequency smoothing and 1/8-octave center
spacing. These are original processing settings, not heat-map blur. Neither a
colour map nor smaller drawn pixels can undo that smoothing. A genuinely
narrower smoothing bandwidth requires a separately authorized new PSD run;
this plotting command does not reprocess raw data or modify the saved spectra.
The network's vertical distribution also contains real between-window and
between-station variability; its width must not be cosmetically compressed.

Frequency on x, pressure PSD on y, colour = occurrence percentage in each
1-dB bin at each frequency. Every saved window contributes one count. Windows
are pooled before normalization, so stations with longer records contribute
more; this is not an equal-station probability or independent-sample count
(the saved PPSD windows overlap). No mean, median or percentile is calculated
or drawn. Different frequency grids are combined on the union of their native
display-bin edges: each saved bin is piecewise constant inside its original
edges, without interpolating PSD amplitudes. At each frequency the denominator
is the number of windows actually covering it, not the whole-network count.
Extreme levels are retained in the saved count grid and full-range PNG;
the detail plot does not renormalize after zooming. Pressure reference bounds
remain overlaid. The new histogram NPZ is derived only; sources stay unchanged.

## Pressure only: every individual spectrum, one colour, no summary lines

```bash
cd /Users/mohammadamin/Projects/bruit-fm
PYTHONPATH=src:. conda run --no-capture-output -n bruitfm python \
  -m bruitfm.plot_network_ppsd --network YV --channels BDH \
  --curves-only --single-color --show
```

This draws all 181,575 saved YV pressure PSD windows, without subsampling or
calculating/overlaying any mean, median or percentile. All observed curves are
blue; reference bounds remain black. It can take time and substantial memory.
These are the individual smoothed spectra already stored by the original PPSD
processing, not newly computed or additionally averaged waveforms/spectra.
Both a labelled pressure detail view and a full-range PNG retaining RR12's
extreme recorded level are saved. Nothing is dropped from the source dataset.

## Optional station-coloured summary mode

Reuse the response-corrected saved PPSD NPZs registered in the master. No raw
data are read, no response/gain is changed, and no master/score is written.
All 172 available YV products are registered: 43 each of BDH, BHZ, BH1 and BH2.

```bash
cd /Users/mohammadamin/Projects/bruit-fm
PYTHONPATH=src:. conda run --no-capture-output -n bruitfm python \
  -m bruitfm.plot_network_ppsd --network YV --channels BDH BHZ BH1 BH2
```

One frequency-domain PNG per channel, in a new timestamped folder under
`results/REVIEWS/PPSD_NETWORK/`. Each station keeps the same colour across
figures. Solid curves are exact per-epoch medians; thin dashed curves are
5th/95th percentiles over **all** saved PSD windows, not a window subsample or
histogram-bin approximation. Multiple epochs remain separate, in the same
station colour. This is a station-coloured comparison of PPSD summaries, not
an occurrence-probability histogram: probability cannot also be the colour.

Pressure additionally gets a `_full_range.png` retaining all extremes. Its main
PNG is explicitly labelled a −100 to +120 dB detail view (customize with
`--pressure-limits LOW_DB HIGH_DB`). This is a display zoom, not a station or
window filter. The actual YV RR12 pressure spectra include an extreme level
near −3,077 dB that would otherwise compress the entire useful display; RR12
is retained in the full-range PNG and manifest. No hardware cause is inferred.

For pressure only, use `--channels BDH`. To additionally draw every individual
saved window, append `--all-curves` (potentially very slow and memory-heavy).
Append `--show` to open figures; headless plotting saves without opening windows.
An explicit `--output-dir` must not already exist. Missing registered NPZs are
reported, not fabricated; malformed or unit-mismatched inputs fail the command.

Pressure is in dB re 1 Pa²/Hz. Brown's low/high hydroacoustic curves are the
repository's approximate Figure 8 digitization, converted from µPa²/Hz to
Pa²/Hz. Current digitized support is about 0.01139-78.10 Hz; no extrapolation
is made. The original 1-31.8 Hz transcription is preserved separately. The
discrete user-supplied points attributed to Webb are not independently verified
published bounds; they are displayed as Wayne-supplied reference points.
Provenance: `src/bruitfm/reference_data/pressure_references.json` and
`brown_hydroacoustic.csv`, DOI `10.1007/s00024-012-0573-6`.

Seismic channels use acceleration PSD, dB re 1 (m/s²)²/Hz, and ObsPy's Peterson
NLNM/NHNM, DOI `10.3133/ofr93322`. These seismic curves must not be overlaid on
pressure axes. The full saved frequency range is displayed; this does not
certify the instrument's response/passband at every frequency. No unverified
pressure-to-acceleration conversion is used.

The station percentiles describe the available observations. They do **not**
constitute a new accepted YV/global low/high noise model: station acceptance,
response uncertainties, deployment/seasonal coverage and environmental grouping
remain separate research decisions. No quality-based station selection is made.
The manifest records exact source paths, full epoch keys, source hashes, time
coverage and window counts. Source NPZ and master CSV hashes are checked before
and after plotting. Original PPSD settings and temporal gaps remain unchanged.
