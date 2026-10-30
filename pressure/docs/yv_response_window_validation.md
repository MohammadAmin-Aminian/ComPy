# YV response conventions and saved-window feasibility

This bounded follow-up to handoff `833cc05` reads saved metadata only. It does
not rerun the YV batch, read raw waveforms, fit gains or change event decisions.
The production report under review remains unchanged. Scientific acceptance is
still pending.

## Explicit analog response check

[FDSN's forward-response convention](https://docs.fdsn.org/projects/stationxml/en/latest/response.html)
uses output/input and evaluates the analog pole-zero expansion at
`s = i*2*pi*f` for radians/second stages. Here the explicit calculation is
`G*A0*product(s-zero)/product(s-pole)`, including the stage gain `G`.
Saved epoch-bound pressure snapshots were checked at 0.03, 0.05 and 0.07 Hz.
RR03/RR22 stage 1 and RR52 stages 1/2 match stage-limited evalresp to a maximum
complex relative error of **4.48e-16**. This checks metadata interpretation;
it does not independently measure sensor sensitivity or hardware polarity.

An algebra-only comparison replaces each positive-real pole with its reflection
into the left half-plane, keeping zeros and gain factors fixed. For these
single-real-pole stages, the amplitude ratio is unity; the phase is different:

| Original / reflected forward response | 0.03 Hz | 0.05 Hz | 0.07 Hz |
|---|---:|---:|---:|
| RR03 phase difference | −28.21° | −17.14° | −12.29° |
| RR22 phase difference | −55.42° | −34.99° | −25.37° |

These are response phase differences, **not newly computed event residuals**.
Replacing the denominator during deconvolution would change recovered pressure
phase by the opposite sign. A pressure-polarity multiplier contributes 180°;
it cannot generally reproduce this frequency-dependent change. Thus matching
amplitudes alone cannot settle the pole signs, and the discrepancy can matter
relative to the existing ±20° phase mask. No reflected response file is created.
RR52 has no positive-real analog poles in its two checked stages.

Stage-limited evalresp emitted sensitivity-mismatch warnings; these were retained.
The isolated stages are not the full counts/Pa chain, so these warnings do not
establish a new full-chain sensitivity failure. Digitizer/FIR normalization and
absolute pressure/Z calibration remain unverified. Obtain the original sensor
and electronics transfer specification before proposing any response correction.

## Can the saved windows be trimmed under the current estimator?

For all 48 candidate windows, the saved first P/S arrival instants partition
each original 1800-second interval. The largest partition gives an upper bound
on the duration available while avoiding those instants. An interval starting
strictly after the last predicted first body arrival was also checked.
Neither calculation excludes body-wave duration or coda.

| Event | Candidates | First-body-arrival flags | Longest arrival-free span across stations |
|---|---:|---:|---:|
| `ak0138esnzr` | 15 | 0 | 1800 s |
| `usb000g7x7` | 15 | 15 | 1301.6–1538.9 s |
| `usb000ief9` | 2 | 2 | 1591.5–1596.9 s |
| `usb000jyiv` | 3 | 3 | 1286.5–1359.8 s |
| `usc000f40j` | 13 | 0 | 1800 s |

**All 20 flagged windows lack a 1600-second arrival-free span within their
original bounds.** All 20 also lack 1600 seconds after the last predicted first
body arrival within those bounds. Trimming alone cannot both avoid those instants
and retain the current eight nonoverlapping 200-second segments. Changing bounds
outside the original window or changing the estimator would require new waveform
and scientific validation; no replacement bounds or shorter estimator are adopted.

The saved predicted group-arrival ranges span about 38–130 seconds. These ranges
describe model dispersion across the display band, not observed packet duration.
They cannot justify a 1600-second isolated-mode interval. The other 28 windows
remain candidates for waveform-based isolation review, not certified pure modes.

## Reproduction and preservation

Run the bounded audit from the project root:

```bash
OPENBLAS_NUM_THREADS=1 /opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python scripts/validate_yv_conventions.py
```

Each invocation writes a new `conventions_*` folder under
`results/VALIDATION/YV_rayleigh/`. The initial output is
`conventions_22bf4b45f61d/validation.json`; its curated numerical record is
[yv_response_window_validation.json](yv_response_window_validation.json).
SHA-256 checks preserved all 2650 selected saved result files, the Rayleigh
workflow, pressure snapshots and both masters. Raw records and the original
inventory were not written. The bounded PAZ assertions and Ruff check/format
passed; no production processing changed.

Next inspect an existing unflagged distant event such as RR52/`ak0138esnzr` for
packet duration and coda, using its saved overview before choosing any new bounds.
In parallel scientific work, resolve the original RR03/RR22 analog transfer
specifications and digitizer conventions. Held-out fits remain premature until
the response and window assumptions are defensible.
