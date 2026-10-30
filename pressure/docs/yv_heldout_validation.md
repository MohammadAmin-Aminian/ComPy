# Frozen-gain validation on unused YV earthquakes

Run the four-station validation:

```bash
cd /Users/mohammadamin/Projects/bruit-fm
PYTHONPATH=src:. conda run --no-capture-output -n bruitfm python \
  scripts/run_yv_heldout_validation.py --open-plots
```

Default stations: **RR01, RR06, RR48**, with **RR52 as the control**. Frozen
pressure gains are respectively 1.04137996, 1.13967879, 1.38671034 and 0.67872781,
from the completed pole-only/original-fit trial. They are not refitted on the
validation events or applied to any waveform, XML, PPSD, master or score.

For a metadata/selection check without reading waveform samples:

```bash
PYTHONPATH=src:. conda run --no-capture-output -n bruitfm python \
  scripts/run_yv_heldout_validation.py --dry-run
```

## Event independence and selection

Use the existing cached USGS catalogue. Exclude every recorded event from the
entire YV trial, including QC failures and selection exclusions, plus the native
software-development probe recorded in
[the development exclusion list](yv_heldout_development_events.json).
Also exclude origins within 24 h of those events, preventing alternate catalogue
IDs or nearby shocks from leaking into the validation set.

Retain Mw-type magnitude ≥6, depth <150 km and the 24 h mainshock rule. Choose
largest events first for cluster resolution, then evaluate **every** selected
unused mainshock: no three-event cap, early success stopping or amplitude-based
event selection. `--minimum-magnitude 6.5` is an explicit stronger-event variant;
changing it defines a separate study. Selection is saved before sample reads.

## What is checked

- Verify the frozen gain against its training candidates and experiment summary.
  Refuse changed masters, nominal response snapshots, inventories or the reference
  dispersion model. Record sensor labels, serial fields, pressure units and the
  nominal/alternative poles in `response_audit.json`.
- Hardware metadata are an audit trail, not independent manufacturer calibration
  certificates. Physical identity verification remains explicitly pending.
- Keep the original group-velocity packet definition, no tilt, native sample-grid,
  epoch, continuity, coherence and phase rules. Pressure-pole reflection changes
  cross-spectral phase in memory only; nominal PSDs and coherence are unchanged.
- Evaluate `g_frozen * measured_P/a / predicted_P/a - 1` by frequency and event.
  Large amplitude errors remain in usable data instead of being removed or refitted.
- Flag other predicted earthquake packets and IASP91 P/S-family arrivals inside
  the window. The primary summary excludes these overlaps; a separate QC-only
  summary retains body-arrival review cases. Visual mode isolation remains needed.
- Report event error spread and, with at least five usable events, an explicitly
  conditional IID-event bootstrap interval for median prediction error. This is
  not an absolute calibration confidence interval: response/model systematics,
  shared estimator assumptions and serial dependence remain unresolved.

## Outputs and preservation

Every invocation creates a fresh `results/VALIDATION/YV_HELDOUT/run_*` folder.
An existing explicit `--output-dir` is refused. Outputs include:

- `selection.json`: all selected and excluded IDs/reasons.
- `heldout_summary.csv` / `.png`: frozen-gain prediction errors by station/event.
- `heldout_events.csv`: every attempt and its QC/isolation/refusal status.
- Station-level `validation.json`, `response_audit.json`, new spectra NPZs and
  per-event review plots.
- Root `validation.json`: completion status, metadata/source hashes, raw-file
  size/mtime records and limitations.

Source metadata/response/master hashes and indexed raw-file size/mtime are
rechecked after processing. This is a bounded preservation check; raw byte hashes
are not claimed. Original results and all scores are untouched. No catalogue
download, PPSD stage, complete YV batch or instrument correction is performed.

Fewer than three usable isolated events produces `insufficient_heldout_events`,
not a pass. Successful script completion means that the study ran; it does not
certify calibration. RR55 is deliberately outside the default station set until
its hardware/response identity is resolved.

The native RR01 development probe produced spectra and failed coherence QC;
it is preserved separately under `YV_HELDOUT_CODE_CHECKS` and excluded from the
formal event set. It did not change any frozen gain or scientific threshold.


## Software checks before the user run

173 tests passed, including six independent frozen-gain/leakage tests. Ruff
check/format passed. Native RR01 samples were processed and saved to a separate
code-check folder with source hashes and raw-file stats unchanged. The formal
preflight excludes 42 prior YV events plus the development probe and selects
**50 unused mainshocks** with the default Mw≥6 policy, without a count cap.
The complete four-station validation is left for the user to run.
