# RR52: processing settings and results

This document describes the implemented pilot
and its verified outputs. The numbered development history is in the
[research log](research_log.md). Exact values are also preserved in the
[run summary](rr52_run_summary.json).

## Scope and data

The pilot uses the local RR52 waveform subset and the matching `YV.xml`
StationXML inventory. The master table contains 16,784 source channel epochs;
this pilot processed the four `YV.RR52.00.*` channels from the larger research archive.

The inventory locates RR52 at latitude -20.4723° and longitude 68.1094°. Its
negative station elevation gives an estimated depth of 2,918 m below sea level.
This estimate has not been checked against bathymetry.

| SEED ID | Component | Sensor described in inventory | Input units | PSD reference units |
| --- | --- | --- | --- | --- |
| `YV.RR52.00.BDH` | Pressure | DPG | Pa | Pa²/Hz |
| `YV.RR52.00.BH1` | Horizontal 1 | Trillium240OBS | m/s | m²/s⁴/Hz |
| `YV.RR52.00.BH2` | Horizontal 2 | Trillium240OBS | m/s | m²/s⁴/Hz |
| `YV.RR52.00.BHZ` | Vertical | Trillium240OBS | m/s | m²/s⁴/Hz |

Horizontal labels do not establish cardinal orientation. No component rotation
was applied. All channels have a sampling rate of 62.5 Hz and a Nyquist
frequency of 31.25 Hz. The selected metadata epoch is recorded in the
[run summary](rr52_run_summary.json).

There are 329 daily MiniSEED files per channel, 1,316 total, occupying approximately
9.33 GB (decimal units). Actual file-header bounds determine waveform selection,
including samples outside nominal daily filename boundaries. The processing
interval uses an exclusive end one sample after the final waveform sample;
exact bounds are retained in the run summary and local processing records.

The StationXML inventory SHA-256 is:

```text
58afc2d6f0056e52cf954bdd5b9fb80e92f9572d44d067b30966a160d2044494
```

The source archive's formal dataset citation/accession is not established in
the current project record and still needs to be confirmed.

## Waveform selection and response correction

Files are indexed by their actual MiniSEED trace headers. Identifier consistency,
time bounds, and sampling rate are checked. Reads select bounded half-open UTC
intervals and merge adjacent/overlapping data using ObsPy method 0 at the
declared sampling rate. Gaps and conflicting overlaps are masked and split;
they are not interpolated or filled with zeros for spectral estimation.
Source data files are not rewritten or resampled.

The master row retrieves the matching full StationXML channel epoch. Processing
requires response stages, finite nonzero signed sensitivity, recognized physical
input units, digital-count output units, and a matching waveform sample rate.
The tabulated sensitivity alone is not substituted for the complete response.

Waveform samples remain in counts when passed to ObsPy PPSD. Instrument-response
correction occurs internally once. Pressure units select
`special_handling="hydrophone"`, including the DPG, avoiding the seismic
differentiation step. Normal seismic handling produces acceleration power
spectral density. Stored spectra are `10 log10(PSD)` in dB relative to the
reference units listed above. This is nominal response correction, not an
independent verification of instrument calibration.

ObsPy's spectral estimator is related to the probabilistic noise methodology of
[McNamara and Buland (2004)](https://doi.org/10.1785/012003001).
The implementation details are described in the
[ObsPy PPSD documentation](https://docs.obspy.org/packages/autogen/obspy.signal.spectral_estimation.PPSD.html).

## Spectral parameters and temporal interpretation

| Parameter | RR52 value |
| --- | --- |
| Outer PSD window | 10,800 s = 3 h; 675,000 samples |
| Outer overlap | 50% |
| Nominal spacing between PSD starts | 5,400 s = 90 min |
| FFT size | 131,072 samples |
| FFT section duration | 2,097.152 s ≈ 34 min 57 s |
| Internal FFT overlap | 98,304 samples = 75% |
| Full-window Welch sections | 17 |
| Lowest nonzero FFT frequency | 0.000476837158203125 Hz |
| Period smoothing width | 1 octave |
| Period-center step | 1/8 octave |
| Requested lower frequency | 0.001 Hz |
| Actual lowest saved frequency center | 0.00095367431640625 Hz |
| Actual highest saved frequency center | 31.25 Hz |
| Saved frequency/period centers | 121 |
| Histogram amplitude range | -250 to 100 dB |
| Histogram bin width | 1 dB |
| Gap handling | `skip_on_gaps=True`; masked gaps split before PPSD |

For this implementation, ObsPy uses Welch averaging with linear detrending and
a cosine taper within each FFT section. The stored spectra then use octave
smoothing. The 3-hour duration is the averaging support of each saved PSD; the
approximately 35-minute internal sections do not create extra time frames.

A spectrogram column is placed at the saved PSD start and spans one nominal
processing step. It represents a 3-hour window despite the 90-minute spacing.
Adjacent windows overlap and are statistically correlated. Increasing overlap
gives more densely sampled starts while retaining the same averaging duration.

The reader handles daily-file fractional-sample phases by advancing the next
read boundary relative to the accepted PSD start. Source timestamps are not
changed on disk; in-memory merges follow ObsPy's sample-grid rules. Verification
found no unintended whole-window gaps at read boundaries.

FFT support is checked against the full lowest-frequency smoothing bin, not
just the requested plotting axis. The maximum required period-bin edge is
approximately 1,482.9104 s, within the 2,097.152 s FFT support. This computational
support does not establish an instrument's usable low-frequency cutoff.

## Execution, restart, and preservation

The final full run used **four worker processes**, one per channel. Each worker
processed its time chunks sequentially. The approximate read size is one day
plus the overlap needed for complete PSD windows at chunk boundaries. Explicit
`OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` were used during the final run.

NPZ/provenance checkpoints are saved every ten attempted chunks and at completion.
Each output directory contains a `ppsd.npz`, an epoch-specific `inventory.xml`
response snapshot, and `processing.json`. The latter records configuration,
ObsPy version, source inventory hash, waveform paths/sizes/mtimes, completed
chunks, warnings, errors, and summary statistics. Waveform content is not hashed.

Restart constructs a PPSD with the correct channel statistics and inventory,
then imports existing spectra with `add_npz()`. This was regression-tested after
the alternative `load_npz(metadata=...)` path gave a response-lookup failure in
the installed ObsPy 1.5.1 environment. Changed settings or source-file attributes
require a deliberate recalculation or a separate output location.

Only the parent updates the four master rows, matched by SEED ID, channel epoch,
and source inventory. Result fields include paths, counts, units, configuration,
frequency bounds, histogram diagnostics, dates, processing status, and errors.
The updater preserves existing calibration scores, other rows, Excel sheets,
and formatting, and backs up both master files before replacement.

## Verified results

| Channel | Daily files | Saved PSDs | Values below/above histogram bounds | Recorded chunk errors |
| --- | ---: | ---: | ---: | ---: |
| BDH | 329 | 5,249 | 0 | 0 |
| BH1 | 329 | 5,249 | 0 | 0 |
| BH2 | 329 | 5,249 | 0 | 0 |
| BHZ | 329 | 5,249 | 0 | 0 |
| Total | 1,316 | 20,996 | 0 | 0 |

Each saved spectral array has shape **5,249 × 121**. All stored values are finite;
PSD start timestamps are unique. The complete-window bounds are recorded in
the run summary. The remaining waveform tail cannot form another
complete 3-hour window at the next nominal start. Final manifests contain no
recorded processing warnings or chunk errors.

The initial histogram floor of -200 dB did not contain the lower seismic tail.
The complete calculation was rerun with a -250 dB floor. All individual PSD
values lie within the final histogram bounds.

Four NPZs reload successfully and their histograms use all 5,249 saved windows.
The original master CSV values and Excel cells were compared with the original
backups. Waveform source sizes and modification times matched their recorded
values after processing; the StationXML hash was unchanged. These checks are
not a full bytewise hash verification of the waveform archive.

The original processing suite passed 14 tests covering metadata/export preservation,
gap handling, response scaling, FFT support, epoch boundaries, checkpoint
extension without duplicate windows, failure/retry behavior, fractional file
timestamps, and parallel channel updates. Ruff lint/format checks passed.
Real-data plot verification confirmed all curves, exact spectral values,
frequency/time orientation, and seismic-only model selection.
The additional expected-pressure comparison and its checks are described in
[pressure comparison](pressure_comparison.md). The subsequent
[infragravity comparison](infragravity_comparison.md) adds model diagnostics;
the expanded suite passes 33 tests.

## Figures and reference models

The [twelve diagnostic figures](figures/rr52/README.md) provide three views of
each channel. Distribution colors show occurrence probability; spectrogram
colors show PSD amplitude in dB. The line views draw every individual saved
smoothed spectrum, not a subsample or a single averaged curve.

Seismic distributions and curve plots show the NLNM and NHNM from
[Peterson (1993)](https://doi.org/10.3133/ofr93322), accessed through ObsPy. The
models are drawn only over their available frequencies within the data band;
they are not extrapolated to Nyquist. Pressure views retain pressure units.
Model overlays are references for visual comparison, not implemented pass/fail
criteria or ocean-bottom pressure-noise models.

Suggested caption for a seismic distribution:

> Probabilistic acceleration power spectral density for YV.RR52.00.BHZ, computed
> from 5,249 overlapping 3-hour windows with 50% overlap. Spectra are corrected
> using the matching StationXML response and smoothed over one octave at
> 1/8-octave steps. Color represents occurrence probability. The NLNM and NHNM
> are shown as seismic reference models (Peterson, 1993).

Suggested caption for the time view:

> Temporal evolution of response-corrected power spectral density at RR52.
> Each column is located at the UTC start of a 3-hour PSD window, with nominal
> 90-minute spacing; color shows amplitude in dB. Adjacent windows overlap.

Suggested caption for the individual curves:

> All 5,249 individual response-corrected, octave-smoothed PSD curves for the
> selected RR52 channel. Each curve represents one 3-hour window; transparency
> shows overlapping curves without assigning a probability color scale.

## Processing summary

Daily MiniSEED waveforms for four channels of station YV.RR52 were indexed using
trace-header identifiers and time bounds and read in bounded intervals. The
matching StationXML channel epoch supplied the complete instrument response.
Probabilistic power spectral densities were calculated with ObsPy from raw
counts, with response correction applied internally once. Pressure spectra
were expressed in dB relative to 1 Pa²/Hz, and seismic spectra in dB relative to
1 m²/s⁴/Hz. Calculations used 3-hour windows with 50% overlap, one-octave
smoothing, and 1/8-octave period steps, with a target band extending from
approximately 0.001 Hz to the 31.25 Hz Nyquist frequency. Gaps were excluded
rather than interpolated or zero-filled. Four channel processes produced
5,249 complete spectra per channel, saved as reusable NPZ data with response
snapshots and processing provenance. Diagnostic views included probability
distributions, time spectrograms, and all individual PSD curves; seismic
distributions and curve plots included the Peterson NLNM/NHNM references.
Instrument usable-band verification and calibration acceptance tests remain
pending, so these results constitute a processing pilot.

## Reproduce the workflow

Use Python 3.13 or newer and install the project with `pip install -e '.[dev]'`.
On the current machine, the working conda environment is named `bruitfm`.

```bash
# Create a master inventory for a new run; keep existing scored exports separate.
python -m bruitfm.inventory \
  --inventory-dir /path/to/Data_Bruit_FM/INVENTORIES --output-dir exports

# Calculate the default full-band product.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m bruitfm.ppsd \
  --data-dir /path/to/Data_Bruit_FM/RR52 --workers 4

# Generate the twelve figures from saved spectra without opening GUI windows.
python -m bruitfm.plot_ppsd --no-show

# Run automated checks.
python -m pytest
ruff check .
ruff format --check .
```

The measured environment was Python 3.13.15, ObsPy 1.5.1, NumPy 2.5.3,
pandas 3.0.6, Matplotlib 3.11.2, openpyxl 3.1.5, pytest 9.1.1, and Ruff 0.16.9,
on macOS 27.0.1 ARM64. The dependency list in `pyproject.toml` is not a version
lock; these versions describe the environment used for this pilot.

## Remaining questions and temporal-resolution alternatives

No earthquake/event rejection, six-test calibration scoring, pressure–vertical
cross-spectral analysis, coastal-distance/seafloor classification, or combined
accepted-station noise model has been completed. Instrument cutoffs and
inventory calibration accuracy require verification. These gaps must remain
explicit when interpreting the pilot.

A separate short-window product could improve the time view while preserving
the current full-band spectra. Checked examples at 62.5 Hz and the current
smoothing are 1 hour/0.003 Hz, 30 minutes/0.01 Hz, and 15 minutes/0.02 Hz, with
50% overlap. They were not used for the current results. Two-hour windows fail
the current 0.001 Hz full-bin support check. Shorter-window spectra require
reprocessing raw waveforms; they cannot be recovered from the 3-hour averages.

## References

- McNamara, D. E., and Buland, R. P. (2004). *Ambient Noise Levels in the
  Continental United States*. BSSA, 94(4), 1517–1527.
  [doi:10.1785/012003001](https://doi.org/10.1785/012003001).
- Peterson, J. R. (1993). *Observations and modeling of seismic background
  noise*. USGS Open-File Report 93-322.
  [doi:10.3133/ofr93322](https://doi.org/10.3133/ofr93322).
- Beyreuther, M., Barsch, R., Krischer, L., Megies, T., Behr, Y., and
  Wassermann, J. (2010). *ObsPy: A Python Toolbox for Seismology*. SRL, 81(3),
  530–533. [doi:10.1785/gssrl.81.3.530](https://doi.org/10.1785/gssrl.81.3.530).
- Krischer, L., Megies, T., Barsch, R., Beyreuther, M., Lecocq, T., Caudron, C.,
  and Wassermann, J. (2015). *ObsPy: a bridge for seismology into the scientific
  Python ecosystem*. Computational Science & Discovery, 8(1), 014003.
  [doi:10.1088/1749-4699/8/1/014003](https://doi.org/10.1088/1749-4699/8/1/014003).

Bibliographic metadata were checked against the
[ObsPy citation list](https://docs.obspy.org/citations.html) and the
[USGS publication records](https://pubs.usgs.gov/). BibTeX entries are in
[references.bib](references.bib).
