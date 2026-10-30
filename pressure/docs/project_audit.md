# Repository and saved-result audit

The completed calculations are internally consistent under their recorded inputs.
The audit found and corrected **five software issues**, with **130 tests passing**
after the changes. This establishes computational checks, not final physical
calibration acceptance. No completed science batch was rerun, no gain was applied,
and no saved station result or scientific score was changed.

## What was checked

The starting checkout was clean at `c41d0ad`; its latest three GitHub Checks runs
had passed. Review covered the numbered development history, continuation notes,
scientific method/assumption documents, processing and calibration code, tests,
runner scripts, package/CI configuration and saved workflow state. Source review
focused on mathematical calculations, response/epoch selection, time/grid guards,
restart dependencies and writes that could affect protected products.

A read-only audit checked:

| Check | Scope | Outcome |
|---|---:|---|
| Tracked Python / JSON / Markdown parsing | 34 / 23 / 28 files at the starting checkpoint | Parsed; no broken local document links |
| Workflow artifact references | 14,174 references, including repeat references | All present; saved selected-scope status counts match |
| Master-linked complete PPSDs | 222 channels; 980,630 saved segments | Finite arrays, unique ordered starts, matching IDs/rates/counts and requested bounds |
| Saved Rayleigh station runs | 116 runs across original YV, no-tilt YV, YL pilot and all-YL workflows | Counts and station medians agree |
| Rayleigh event spectral products | 2,219 saved spectra; 408 candidate occurrences | Power ratios, complex transfers, squared coherence, phase, masks and closed-form gains agree |
| Other saved diagnostics | 43 pressure-bound comparisons, 43 Godin comparisons, 41 raw Godin pilots, one paired-spectrum run | Percentiles, bound fractions, residuals, response-power division and linear/complex averages agree |
| CSV / Excel agreement | 16,884 rows; 2,363,760 compared cells | Agree, allowing numeric formatting and split long text |
| Waveform provenance | 61,694 distinct referenced source files | Recorded sizes/mtimes match |
| Current preservation baseline | 93,784 production/source files | No size/mtime change during the audit |
| Historical YV protection | 2,648 original result/workflow/snapshot hashes | Match; the two old master hashes changed with authorized YL registration |
| Completed YL comparison protection | All 28 protected hashes | Match, including current masters |
| YL full-interval PPSD check | All 50 complete P/Z products | Header-derived available bounds match full requested intervals |

Run counts include separate saved calculations of the same station/events; they
are not counts of independent observations. Fifty additional YL master rows are
horizontal LH1/LH2 metadata, not missing P/Z jobs. All six score columns are blank.
The two original YV Godin pilot failures remain recorded: RR12's selected pressure
windows are flat; RR42 lacks enough coverage for seven days.

The array audit uses existing spectral inputs and nominal responses. It does not
repeat raw response removal or establish absolute response accuracy. Raw waveform
contents were not exhaustively rehashed: preservation uses source size/mtime
provenance, baseline stats, and hashes of masters/inventories and protected saved
products. Plot verification checked referenced files and sampled the saved A02W,
C07W and B01W comparison figures; it did not visually reinspect every plot.

## Five corrections and their implications

1. **Short PPSD pilot reuse.** A completed bounded interval could pass the YL
   full-archive reuse check. Preparation now checks the requested start against
   the full origin and the requested end against all source headers, clipped to
   the response epoch. Regression cases cover a late start and an early end.
   The existing 50 YL P/Z products already cover their full available intervals;
   this issue did not truncate the completed YL run.
2. **Pressure-unit selection.** PPSD accepted `PA`, `Pa`, `PASCAL` and `PASCALS`,
   while Rayleigh selection required exactly `PA`. Selection now shares one
   pressure-unit predicate, including case/space normalization. Tests cover both
   standalone selection and the network path. Current saved YV/YL pressure rows
   use `PA`, so this does not change their results.
3. **Explicit depth preflight.** The network's resonance preflight used master
   depth even when an explicit depth was supplied to calibration. It now uses the
   same effective depth as the fit. Existing runs did not supply that override.
4. **Coordinate-dependent restart.** Editing master latitude/longitude did not
   invalidate cached Rayleigh geometry. Coordinates now enter the dependency
   fingerprint; regression cases cover both latitude and longitude.
5. **Auxiliary fit failure isolation.** A failed H1 or unrotated-Z diagnostic fit
   could discard a valid primary P/a candidate. A failed median-curve display fit
   could also prevent the whole station result from being saved. These failures
   now retain explicit diagnostic reasons while preserving eligible primary
   results. Primary coherence, phase, overlap, timing and gain-bound safeguards
   remain in force. Synthetic failure-injection tests cover all three diagnostics.
   This does not retroactively accept previously unresolved events or enlarge
   B01W's primary search bounds.

Core response correction, water-column equation, primary least-squares objective,
gain bounds, event selection and timing tolerances were not changed. Existing
saved plots can show a gain fitted to the median curve, whereas the station table
uses the median of event gains; those are documented distinct estimators.

## Scientific questions still open

A follow-up checked the recorded histogram-range flags: **684,714 values across
38 channels** lie outside −250 to +100 dB. Stored spectral arrays retain these
values, while ObsPy collects histogram under/overflow in the edge bins. RR12
pressure and all four RR50 channels include many values near **−3076.53 dB**, the
numerical zero-power floor. These are not measurements of extremely quiet ocean
noise. Other pressure channels have some values above +100 dB. Original spectra
and statuses remain preserved; zero-power/transient screening and histogram
interpretation need review before constructing accepted noise statistics.
This does not alter the saved raw-event calibration fits, which do not use PPSD
amplitudes. The additional counts are retained in the numerical audit sidecar.

The primary [Zha & Webb paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012688)
was rechecked against the transcribed Table 1 and Appendix A. Reciprocal comparison
aligns correction direction; it does not establish identical fitting objectives
or nominal responses. Our 0.03–0.07 Hz borrowed model and event/averaging policy
still differ from their 0.04–0.1 Hz regional-model procedure. The available YL
comparison remains 21 reciprocals 3.4–7.0% lower, C07W +86.8%, three unresolved.
See [the station comparison](yl_all_zha_comparison.md).

The [Godin reference](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/grl.50418)
provides an empirical infragravity model, not an instrument acceptance envelope.
Webb handout bounds and Brown digitization/nominal sensor-band applicability still
need source and field review; no acceptance rule was inferred from those plots.

Remaining priorities are C07W's response/matched-event/window investigation,
B01W's primary 0.1–3 gain-search limitation, A14W/C10W timing/coherence/coverage,
and YV positive-pole/absolute-response and mode-isolation questions. Across the
saved workflow candidate occurrences, 214 windows contain a saved predicted first
P/S arrival; occurrences repeat across runs. These are review flags, not proof
of contamination or automatically rejected data. The original YV audit already
identified 20 such windows among its 48 candidate events. Hardware, bathymetry,
clock accuracy, held-out validation and total uncertainty remain unverified.
No gain should be presented as a physically validated calibration on the basis
of this software audit alone.

## Reproduce the saved-product audit

From the checkout, with the local products present:

```bash
conda run --no-capture-output -n bruitfm python scripts/audit_saved_project.py
```

This writes a new `results/VALIDATION/PROJECT_AUDIT/run_*/audit.json`, reads saved
products/provenance and exits nonzero if checks fail. It does not invoke waveform
processing, update masters or rerun YV/YL. The initial audit is `run_5af6e0fc9900`;
the expanded diagnostic audit is `run_843492622eee`. The header-only YL coverage
check is retained beside the initial audit. [Curated numerical record](project_audit.json).

Software verification: **130 tests passed**, Ruff check/format, Bash syntax and
Git whitespace checks passed. The retained test warnings describe ObsPy's inability
to map pressure units onto seismic motion units; response handling tests confirm
pressure follows its declared pressure response without seismic differentiation.
The pushed fixes also passed [GitHub Checks](https://github.com/MohammadAmin-Aminian/bruit-fm/actions/runs/37952188324).
