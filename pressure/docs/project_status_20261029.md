# Project status, 29 October 2026

Start here next session, together with
[`integration/README.md`](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/integration/README.md).

**30 October.** BRUIT-FM 1.3.0 is the last release: the owner moved
everything into ComPy (3.0.0, PR #10), and this repository is kept as a
backup with the dated reports. See `integration/README.md`. The notes below
are from 29 October.

**Evening update.** The owner asked for the integration without waiting for
the supervisor:

- BRUIT-FM is now 1.3.0, uncommitted on `main` (not tagged).
- Its core is ComPy's `compy_pressure`, on the uncommitted ComPy branch
  `pressure-integration` (worktree `~/Desktop/BRUIT-FM/ComPy-pressure-integration`).
- All tests pass on Python 3.10 and 3.13.
- In the reproduction check, every QC and noise-model output matches
  BRUIT-FM and the archived 28 October results.
- What remains (commit, rebase onto ComPy 2.4, PR, the supervisor's answers):
  `integration/README.md`.

Nothing is running.

## Today

| Item | Result | Where |
|---|---|---|
| Report to the supervisor | 9-page PDF: gain-corrected noise model (map, all groupings as heat maps and medians, before/after gain), selection criteria, provisional decisions A–G and questions 1–7; sent by the owner with a short email | `docs/reports/2026-10-29_noise_model_report/` |
| Shallow-water peak explained | 0.05–0.1 Hz peak in the < 1 km class = ocean swell at 13 Cascadia shelf stations (56–173 m); passes QC test 1 because the swell band is 2 of 31–54 bins | report section 3.3 |
| Godin with and without tanh | with tanh within ±4 dB of the depth-class medians at 5 mHz; without, 1–18 dB too low | report appendix |
| Integration preparation | core separated (22 modules), manifest, boundary test, export tool (225 passed, 3 skipped), comparison tool, rehearsal (12/12 noise-model files identical) | `integration/` |

| Godin with a different B<sub>p</sub> | "bp" is B<sub>p</sub>, the amplitude of Eq. 6 (2.0·10³ Pa²/Hz). Fitted per depth class: 2.3–5.0·10³ (+0.7 to +4.0 dB); fitting removes only 0.1–1.3 dB RMS. The model does not apply below 200 m water depth | report v2, Figure 15 |

The paper is now in `~/Desktop/BRUIT-FM/Papers/Godin_2013_GRL.pdf`. Report
v2 (10 pages) adds Figure 15 and asks about B<sub>p</sub> in question 7; v1 is
kept as sent.

## Owner's decisions today

- BRUIT-FM's core replaces ComPy's `Pressure_calibration.py`.
- Tilt is not moved; ComPy has its own.
- Do not touch ComPy until the owner asks; integrate after the supervisor's
  answers.

## Waiting for

1. The supervisor's answers (record them in `integration/README.md`). Apply
   them in BRUIT-FM, export again, rerun the reproduction check, then release
   v2.0.0 as the frozen campaign record.
2. The owner's go-ahead to commit both repositories and rebase the ComPy
   branch onto `main` after the 2.4 merge.

Earlier status: [`project_status_20261028.md`](project_status_20261028.md).
