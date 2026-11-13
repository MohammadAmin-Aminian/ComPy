# ComPy 3: a new default Vs inversion, with its audit

**Vs results change. Compliance estimates do not.** Not for merging before the author has read `docs/inversion_overhaul_2026-11/REPORT.md` (PDF in the same folder).

## Why

An audit of the inversion against the primary literature and against tests that use no earlier ComPy result found:

- **Errors.** The 1 % independent floor per frequency was five times the real scatter of a stacked curve (0.2-0.4 %, correlated between bins) and too strict on the slope of the curve. Median chi-square per datum was 0.38; it is 1.08 now.
- **Prior.** Flat layers with second-difference smoothing pushed unresolved deep layers to their bounds.
- **Vp and density.** Compliance measures stiffness, `(1 - nu)/mu` for a half-space. One free Vp/Vs settled near 2.1 (lower-crust Vp above 8.5 km/s) and lowered Vs by about 10 % at 2-5 km.
- **Sediment.** A thin soft layer adds `k h/(rho Vp^2)` to compliance whatever its Vs. Density was fixed at 2.0 g/cc and thickness limited to 1 km.
- **Band.** At shelf stations the estimate collapses where the pressure signal no longer reaches the seafloor, and the saved curve stopped at 40 mHz.

The forward model (4e-7 against a 50-digit solver), the gravity corrections, the sampler and the H1 estimator were checked and are right.

## What

- `compy.inversion.auto` (`invert`, `SeafloorModel`, `compliance_data`), `compy.inversion.errors`, the `compy-invert` command, `compy.plotting.auto.station_figure`.
- Stage C saves the bootstrap covariance of each stack; stage E runs the new inversion by default; figures and report stages read either result.
- The 2.x inversion is unchanged and repeats an archived 2.8 run draw for draw.
- `docs/INVERSION.md` rewritten; user guide, README, example 4 and changelog updated.
- 18 new tests; 940 pass.

## Evidence

- Calibration: on synthetic data from the prior, 50 % and 90 % intervals hold the truth at those rates on four station types; with an Earth the layers cannot represent, intervals stay honest below 700 m.
- Archive rerun (287 stations): tier 1 135 good, 20 usable, 17 poor. Vs against 2.8: +6 % at 2 km, +7 % at 3 km, +11 % at 5 km (medians of 155 stations); sediment thickness x0.76.
- Independent map: soft-sediment thickness against GlobSed, rank correlation 0.66 over 92 thin-sediment stations.

## Open (section 10 of the report)

The slope allowance and prior widths are judgements (prior widths change the intervals, not the medians); Vs between 0.5 and 3 km depends on the upper limit of the band, most at thick-sediment stations; the fast lower crust at four YV stations; sources not obtained in the original (Crawford 1991/1998, Ruan 2014, Bell 2015, Brocher 2005).
