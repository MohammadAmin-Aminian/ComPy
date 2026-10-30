# Handoff: Claude microseism check → Codex

Branch `origin/claude/microseism-all` (based on `origin/main` through `ed5cb16`).
Status: diagnostic only. No gain was fitted, applied or stored. No master,
score, response, raw file or saved product was changed.

## Read first

1. This file.
2. [microseism_all_networks.md](microseism_all_networks.md): commands, four
   columns, results, and the repeated-week findings.
3. [reviews/claude_microseism_self_review.md](reviews/claude_microseism_self_review.md):
   findings SC-1…8 and SS-1…11, with what is still open.
4. [microseism_pressure_vertical_check.md](microseism_pressure_vertical_check.md):
   method, equations, RR52 pilot and primary-source checks.

## Where things are

**Code.**

* `src/bruitfm/microseism_check.py`;
* `scripts/run_microseism_check.py` (runner; `--base-run/--shift-weeks`);
* `scripts/compare_microseism_methods.py` (four columns; `--combine`);
* `scripts/combine_microseism_weeks.py`;
* `scripts/run_microseism_check_all.sh` and `scripts/run_microseism_weeks.sh`;
* tests in `tests/test_microseism_check.py` and
  `tests/test_microseism_weeks.py` (215 pass).

**Local results.** These are gitignored and live in the Claude worktree only:
`~/Projects/bruit-fm-claude/results/EXPERIMENTS/MICROSEISM_ALL/`.

| Network | Base run (week 0) | Shifted weeks done |
| --- | --- | --- |
| YV | `run_20261014T183508Z_1f096dca` | none |
| YL | `run_20261014T183501Z_0a1ff170` | −12, −6, +6, +12 (each twice, identical; the combiner keeps the latest) |
| 7D | `run_20261014T183456Z_61c4966a` | −12, −6, +6, +12 all outside data coverage (useless) |
| XO | `run_20261014T183418Z_736f654a` | none |

**Inputs.** The runner always reads masters and waveforms from
`~/Projects/bruit-fm` (`MAIN` in the runner), whichever checkout runs it. It
hashes the masters before and after.

**Outputs.** Shifted runs are written beside their base run by default, which
is where the combiner looks. Running from another checkout therefore still adds
folders under the path above. Nothing tracked changes.

## Remaining work, in order

### 1. Finish the repeated weeks

Run from any checkout of this branch:

```bash
bash scripts/run_microseism_weeks.sh ~/Projects/bruit-fm-claude/results/EXPERIMENTS/MICROSEISM_ALL/run_20261014T183418Z_736f654a -8 -4 4 8
bash scripts/run_microseism_weeks.sh ~/Projects/bruit-fm-claude/results/EXPERIMENTS/MICROSEISM_ALL/run_20261014T183456Z_61c4966a -2 -1 1 2
bash scripts/run_microseism_weeks.sh ~/Projects/bruit-fm-claude/results/EXPERIMENTS/MICROSEISM_ALL/run_20261014T183508Z_1f096dca
```

The shifts differ by network because the data coverage differs:

* **7D:** the on-disk working copies hold only 33–41 days (2011 days 324–364),
  so use ±1–2 weeks.
* **XO:** WD66 data run from day 198 to 365 of 2018, so ±12 weeks falls
  outside; use ±4–8 weeks.
* **Stations outside a shifted week** stay `unresolved_pairing`. That is
  expected; there is no fallback.

Do not rerun YL.

### 2. Implement the frequency-consistency gain rule (fixed now, before the remaining weeks are seen)

The YL weeks showed that the best-sub-band rule can pick a 0.1–0.2 Hz sub-band
where the mode-0 RR52 model is wrong:

* YL.C01W gives ×2 between weeks while Codex's model-free amplitude is flat;
* 63 of 100 informative station-weeks fail sub-band consistency (SS-10).

The rule below was designed after seeing YL. **Treat YL as the design set and
XO, 7D and YV as the blind test. Do not tune these thresholds after seeing
them.** If a change looks necessary, report both versions.

**Per station-week**, read the week's check summary: `stations/NET.STA.json`
→ `check_summary` → `long_wave_sub_bands`.

* **Eligible sub-band:**
  * upper edge ≤ 0.10 Hz;
  * `mode0_water_h1`, `mode0_water_h2` and `upper_bound_h2` all present;
  * conditioning `upper_bound_h2 / mode0_water_h1` ≤ 1.5.
* **The week passes when all of these hold:**
  * the phase label is exactly one of `nominal`, `reflected_pole`,
    `reversed_nominal`, `reversed_reflected_pole` (not `approximate_*`);
  * there are ≥ 2 eligible sub-bands;
  * max/min − 1 of their midpoints (H1 + H2)/2 is ≤ 0.10.
* **Week value:** the median of the eligible midpoints. The week bracket runs
  from the minimum H1 to the maximum H2 over the eligible bands.

**Per station:**

* **Reportable** when ≥ 3 passing weeks exist (base included) and the
  (max − min)/median of the week values is ≤ 10%.
* **Stable** when that spread is ≤ 5%.

Implement this in `combine_microseism_weeks.py`, keeping the current output
for comparison, and add tests. Then run:

```bash
PYTHONPATH=src:. conda run --no-capture-output -n bruitfm python scripts/combine_microseism_weeks.py <YV base> <YL base> <7D base> <XO base> --docs docs
```

Update `microseism_all_networks.md` with the result. The 35-station
single-week gain set there is already marked superseded.

### 3. Open science items (Socrates)

* **YL.C07W.** All sub-bands agree within each week, but the weeks differ by
  ×2 (0.27–0.31 vs 0.54–0.60). The Rayleigh baseline (0.269) and trial (0.526)
  split the same way. Compare the dates of the Rayleigh events with the week
  values and with the response epochs.
* **YL offset.** The YL baseline is 3.5–7.6% above Zha & Webb (2016) at every
  station except C07W. Find their nominal DPG sensitivity, ρw and depth source
  (SS-3).
* **Polarity and response leads:**
  * YV.RR47: polarity reversed against the YV convention;
  * 7D J35A, J36A, J43A, J44A: phase near 0° at 0.1–0.2 Hz, drifting below
    0.1 Hz; probably the response;
  * reflected pressure poles at 22 YV stations: check against the
    pole-reflection trial.
* **Earth model.** Run local Earth models for RR50, 7D and XO, and verify
  depths where the structure-free maximum exceeds the mode-0 value by
  more than 10%.

## Constraints and pitfalls

* **Never apply gains, and never write masters, responses, raw data, scores or
  saved results.** Do not rerun completed calibration batches.
* **Import path.** pytest now imports the checkout's own `src`
  (`pyproject.toml`: `pythonpath`). Scripts still need `PYTHONPATH=src:.`
  (the wrappers set it). The `bruitfm` conda env otherwise imports
  `~/Projects/bruit-fm/src`.
* **Hidden failures.** Do not pipe runs into `| tail` without `pipefail`; it
  hid a failed regeneration once.
* **Codex's own diagnostic.** `run_diagnostic` refuses an existing output
  folder. `compare_microseism_methods.py` reuses a saved diagnostic only when
  its source and multiplier match.
* **Published values.** Deng et al. (2022) values count only when frequency
  matched; WD70 and WD48 are not tests. Aminian et al. (2025) shares the
  Rayleigh method family with column 2.
* **Merging.** Do not merge this branch into `main` without the user's
  approval.
