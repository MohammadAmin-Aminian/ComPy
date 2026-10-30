# Baseline-station microseism commands

This campaign covers the **46 stations with a saved first-method Rayleigh
estimate**: YV 16, YL 22, 7D 4 and XO 4. It intentionally excludes the 15
pole-trial-only estimates, 26 unresolved stations and two reference-band
exclusions. The frozen values are copied exactly from `previous_g` in
`socrates_all_station_comparison.csv`.

Each invocation processes one station. It selects the exact pressure and
vertical channels from
`src/bruitfm/reference_data/microseism_baseline_stations.csv`, writes fresh
paired spectra and a fresh diagnostic, and never requests a master update.
RR52 is the exception: its designated saved paired spectra are reused by
default. No command fits a microseism gain or applies the fixed Rayleigh gain to
waveforms, responses or saved spectra.

To run all 46 stations sequentially and continue past individual failures:

```bash
bash scripts/run_microseism_baseline_all.sh
```

The wrapper writes per-station terminal logs and a `summary.tsv` under a fresh
`results/DIAGNOSTICS/MICROSEISM_BASELINE_LOGS/run_<UTC>/` directory. It exits
nonzero if any station fails, but only after attempting every station.

From the repository root, define this convenience function once in `zsh`:

```bash
ms () {
  PYTHONPATH=src:. conda run --no-capture-output -n bruitfm \
    python -m bruitfm.microseism_campaign --station "$1"
}
```

Then run stations individually. The exact channel pair and frozen multiplier
are printed by `python -m bruitfm.microseism_campaign --list`.

## YV: `BDH / BHZ`

```bash
ms YV.RR03
ms YV.RR07
ms YV.RR08
ms YV.RR09
ms YV.RR11
ms YV.RR13
ms YV.RR14
ms YV.RR25
ms YV.RR28
ms YV.RR29
ms YV.RR34
ms YV.RR36
ms YV.RR38
ms YV.RR40
ms YV.RR50
ms YV.RR52
```

## YL: `LDH / LHZ`

```bash
ms YL.A02W
ms YL.A06W
ms YL.A08W
ms YL.A10W
ms YL.A12W
ms YL.B04W
ms YL.B09W
ms YL.B11W
ms YL.B14W
ms YL.C01W
ms YL.C03W
ms YL.C05W
ms YL.C07W
ms YL.C08W
ms YL.C09W
ms YL.C11W
ms YL.C16W
ms YL.C17W
ms YL.F02W
ms YL.N02W
ms YL.N03W
ms YL.S01W
```

## 7D: `LDH / LHZ`

```bash
ms 7D.J30A
ms 7D.J38A
ms 7D.J45A
ms 7D.J54A
```

## XO

```bash
ms XO.LA33  # HDH APG / HHZ; explicit 125/100 Hz -> 1 Hz working copies
ms XO.LA34  # HDH APG / HHZ; explicit 125/100 Hz -> 1 Hz working copies
ms XO.WD66  # LDH DPG / LHZ; explicit 1 Hz working copies
ms XO.WD70  # LDH DPG / LHZ; explicit 1 Hz working copies
```

The XO resampling is the existing audited polyphase/Lanczos path. It changes
only in-memory working copies, uses native complex responses, fills no gaps and
does not correct logger clocks or rewrite headers.

## Outputs and later comparison

Each successful command prints JSON containing `paired_comparison` and
`diagnostic` paths. Preserve those paths for the combined table. Diagnostic
JSON records the fixed first-method coefficient and the applicable publication
reference:

- YV matched values use Aminian et al. (2025) Table 1 in the same `g`
  direction where available.
- YL uses Zha & Webb (2016) Table 1 in the reciprocal `1/g` direction.
- Ruan et al. (2014) 7D scale factors require the station-specific coupled
  modal reproduction; they are not inserted as guessed scalar targets.
- Deng et al. (2022) XO results are frequency-dependent body-wave values in
  different bands and are not treated as directly comparable scalars.

The cross-station comparison can rank coherence, phase and transfer stability,
and can display the frozen first-method/published references. It cannot call the
microseism result an independent gain until vertical response, water depth and
station-specific structure constrain
`H_nominal = (g_pressure/g_vertical) * eta`.

## Completed comparison

The user completed all 46 jobs in the explicit
`run_20261014T170400Z` campaign. The verified cross-station result, including
the full CSV and figure, is in
[the baseline comparison report](microseism_baseline_comparison.md). No
microseism gain was fitted or applied.
