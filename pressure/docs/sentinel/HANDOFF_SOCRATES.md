# Sentinel → Socrates handoff

Socrates, I am **Sentinel**, the software-engineering reviewer of BRUIT-FM
(architecture, bugs, performance, reliability, testing). This branch holds
everything I reviewed and changed. Scientific validity is yours to judge: treat
my scientific findings as leads to verify or overturn, not as conclusions.

## 1. Getting the branch

```bash
git fetch origin
git switch review/sentinel        # based on faa39ff; merges cleanly onto main 4387fbe
```

- Do not run tests in the shared conda env `bruitfm`: it imports `bruitfm` from
  Codex's checkout, not from this branch. Use a private environment with the
  same versions (ObsPy 1.5.1, NumPy 2.5.3, SciPy 1.18.1, pandas 3.0.6,
  TiSKitPy 2.3.1, disba 0.7.0), e.g. `python -m venv .venv` then
  `pip install -e '.[dev,rayleigh,dispersion]'`.
- Run tests with `HOME` pointed at a scratch folder. Several defaults resolve
  under `~/Desktop/Data_Bruit_FM`, so this keeps a test from touching real data.
- `pytest -q`: 141 passed. `ruff check` and `ruff format --check` passed.

## 2. What is on the branch

| Kind | Commits | Files |
|---|---|---|
| Software fixes, each with a regression test that fails on the previous code | `2d4e8d6` SENT-002, `41a86f5` SENT-003, `04350d5` SENT-013, `63ef818` SENT-014 | `src/bruitfm/rayleigh_calibration.py`, `src/bruitfm/inventory.py`, `tests/` |
| Independent tests of conventions and guards (pass on unchanged code) | `3f663c5` SENT-011 | `tests/` |
| Documentation corrections (warnings, no numbers changed) | `ddb5c35` | `README.md`, several `docs/*.md` |
| Audit record and evidence | `ce41f54`, `da8932a`, `41cafa1`, `e142d81`, latest | `docs/SENTINEL_AUDIT.md`, `docs/sentinel/*`, `scripts/sentinel_yl_review.py` |

None of the fixes changes an acceptance rule, a bound, a sign, a gain or a
saved result.

- SENT-002 keeps the spectra and true QC status when a fit hits the grid
  boundary.
- SENT-003 pairs tilt-correction horizontals by band code.
- SENT-013 stops the inventory export from overwriting an existing master.
- SENT-014 locks master updates across processes.

## 3. Scientific items for your review

Each item has my evidence in `docs/SENTINEL_AUDIT.md` §3. All are open and
none is acted on.

| ID | Claim to test | Where my evidence is weakest |
|---|---|---|
| SENT-001 | YL C07W's low gains come from an intermittent ≈0.5× vertical (LHZ) amplitude state, not the DPG. | Cross-station ratios only; the physical cause and the deployment logs are unknown. |
| SENT-004 | YL B01W and YV RR55 have a frequency-dependent pressure response (P/a phase +45° to +76°, R/T rising with frequency), so a scalar gain is inadequate. | The ≈0.07 Hz extra high-pass is a hypothesis. |
| SENT-005 | The 1 µs between-channel timing tolerance excludes 515 YL windows over 2–200 µs offsets (≤0.005° phase at 0.07 Hz). | Q330 per-channel timestamp behaviour is not documented here. |
| SENT-006 | The −5.4% median offset from Zha & Webb (2016) lies within method sensitivity (density, response treatment). | Their Appendix A processing parameters were not accessible. |
| SENT-007 | ObsPy's −3076.5 dB zero-power floor enters pressure diagnostics percentiles; RR12's "0% in range" is an artefact. | — |
| SENT-010 | Differentiating after tapering biases g by about +0.65% on real events. | The magnitude is from 48 events at six stations. |
| SENT-012 | The network CLI's default pressure sign of −1 (from YV) applies to every network. | — |
| RR22, RR56 | Not analysed. Their positive-real poles predict 55–57° of phase at 0.03 Hz, so they need the same R/T-against-frequency check. | — |

### Conventions I verified, and how

Please re-derive these independently rather than rely on them:

- Water-column reference P/a_up = ρ·tan(rH)/r with r = ω·√(1/α² − 1/c²), for z down, P(0)=0 and ∂P/∂z(H)=ρ·a_up. Checked by derivation and by series expansion.
- The numpy e^{−iωt} FFT convention is consistent with evalresp's s = 2πif; upward Z negates channels with dip > 0. Checked with an analytic time-domain derivative through production `event_spectra`: ratio exact to 1e-9, phase 0 to 1e-6°.
- Gain direction: g multiplies pressure, so 1/g multiplies counts/Pa, and comparison with Zha uses 1/g.
- PPSD: hydrophones have no ω² term; seismic VEL × ω² gives acceleration; ObsPy averages dB within period bins. Checked numerically.
- An independent estimator (SciPy Welch CSD, cross-spectral response division, my own water formula) reproduces 16 saved YL gains within 0.07%. It shares responses, the dispersion table and event windows with production, so it is not a full independent calibration.

### Sources I could not read

- Zha & Webb (2016), full text.
- Zha's thesis.
- Godin (2013) Eq. 6 normalisation.
- The primary source for the Webb bounds.
- The Brown source raster.

Table 1 values come from the repository's transcription.

## 4. Working on the same branch without collisions

- **Sentinel writes only** `docs/SENTINEL_AUDIT.md`, `docs/sentinel/**` and
  `scripts/sentinel_yl_review.py`.
- **Suggested for Socrates:** your own files, e.g. `docs/SOCRATES_REVIEW.md` and
  `docs/socrates/**`. If you correct one of my scientific statements, write it
  in your report and I will update the audit record.
- **Code changes go through Codex.** Neither of us should change `src/` or
  `tests/` on this branch without telling the other in a commit message or
  report. Propose code changes in your report.
- Pull before committing (`git pull --ff-only`). Never force-push or rewrite
  history.

## 5. Standing rules from the user

- Never apply calibration gains to original observations or responses.
- Never overwrite raw waveforms, StationXML, masters, scores or saved results;
  never delete or reorganize research data. Put every experiment's output in a
  new directory.
- Never rerun completed YV/YL batches to establish a baseline. Ask the user
  before any expensive or network-wide real-data calculation.
- No provisional gain is an accepted correction. Do not widen bounds to
  reproduce published coefficients. Do not call guessed layers or the RR52
  model "Zha's model".
- Tilt correction is optional. Do not require horizontals or quiet training
  when the user selects no tilt.
- Passing tests do not make the project scientifically validated.

Decisions pending with the user:

1. Integration of this branch.
2. C07W treatment.
3. B01W/RR55 handling.
4. Timing tolerance.
5. Zero-power screening.
6. Estimator ordering.
7. The network default sign.
8. Master flush cadence and backup retention.
