# BRUIT-FM contribution and review record

This records the people and agents responsible for the completed work through
README step 56. Repository commit authorship is preserved; the role attribution
below documents work performed through the user's coding-agent sessions.

| Contributor | Role and completed work |
|---|---|
| **Mohammad Amin Aminian — project owner** | Scientific direction and original research context; waveform/inventory acquisition; executing the YV, YL, 7D and XO runs and scientific trials; supplying published-method references; reviewing results and authorizing integration. |
| **Codex — development and integration agent** | Original BRUIT-FM implementation with the owner; inventories, processing, diagnostics and network workflows; calibration and APG support; independent regression/raw-window checks; closing the fresh-master lock race; implementing optional scientific experiments; comparing results, maintaining numbered reports and integrating approved branches. |
| **Sentinel — independent Claude engineering reviewer** | Separate code/test/reliability audit; fixes retaining spectra and true QC at gain boundaries, matching tilt horizontal band codes, refusing master overwrite and serializing master updates; independent convention tests, station diagnostics and evidence preservation. |
| **Socrates — independent Claude scientific reviewer** | Reviewing Sentinel's evidence against physics and primary publications; response-pole, periodic-glitch, estimator/model/uncertainty and response-shape analysis; read-only V1–V11 diagnostics and qualified proposals. Socrates' review supplied scientific hypotheses; it did not directly implement production corrections. |
| **Claude — independent microseism implementation** | Water-layer admittance, RR52 pilot, modal calculations, H1/H2 conditional gains, all-network and repeated-week tools; identified the best-band selection flaw and specified the frozen frequency-consistency rule in the handoff. Its Sentinel/Socrates-style self-review is explicitly a self-review. |

For the microseism integration, Codex retained Claude's implementation and
branch history, implemented the handoff rule with provenance/independence
guards, qualified tilt and published-band comparisons, and evaluated saved
weeks without another calibration batch. See the
[integration review](reviews/microseism_integration.md).

## Review and integration history

- Sentinel/Socrates branch `review/sentinel` at `6f074c9` was integrated with
  Codex's additional locking correction in `0943889`.
  [Engineering review and qualifications](reviews/codex_integration_review.md).
- `a751592`: Codex implemented the isolated selected-station scientific trials.
- `e79f773`: extend those trials to all 89 previously tested station epochs.
- `c1513b4`: record the user's completed run and the full station comparison.
  [Results and publication discrepancies](socrates_all_station_comparison.md).
- The project owner then explicitly requested merging, pushing, contributor
  attribution and a final handoff. The final merge retains both branch histories;
  inspect Git for its identifier.

## What the integration establishes

The software suite passes **167 tests**; style checks pass. Numerical trials
produce 61 provisional estimates versus 46 baseline estimates, including C07W's
three-event normal-state subset. That is processing/QC coverage, not a scientific
accuracy percentage. The original fitter stays the production default; alternative
pole conventions and weighted fitting remain in the explicit trial runner.

All original waveforms, input XMLs, response snapshots, masters, saved results
and scores are preserved. No fitted gain is applied. Glitch removal, automatic
hardware/response repair, final calibration acceptance and instrument scoring
are not claimed as completed contributions. The original independent reports
remain attributable to Sentinel/Socrates, with Codex qualifications clearly marked.
