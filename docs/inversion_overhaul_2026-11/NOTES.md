# Working notes (survive context loss; raw material for REPORT.md)

## Literature status (10 Nov 2026)

Read in full or in the relevant sections (text in `~/Desktop/BRUIT-FM/compliance_yv_2026-10-28/literature/papers/txt` and `inv2026-11/`):

| Source | Read | Taken from it |
|---|---|---|
| Crawford (1994) thesis, ch. 2 and 6 | yes (ch. 3-4 are placeholders in this PDF: they are the 1991 and 1998 papers) | Eq. 2.19: half-space compliance xi = -(lambda+2mu)/(2 mu k (lambda+mu)); normalized eta = -k xi. Eq. 2.20-2.22: d eta/dVs is 2-9 times d eta/dVp for 1.8 < Vp/Vs < 2.4. Forward code: Gomberg & Masters (1988) minor vectors. Ch. 6: Vp fixed from refraction, start Vs = Vp/1.87, density from Carlson & Raskin (1984): rho = 3.50 - 3.79/Vp (Vp < 6.65 km/s), rho = 3.81 - 5.99/Vp (Vp >= 6.65). Occam inversion (Constable et al. 1987), second-derivative smoothing. Sensitivity near lambda/6. |
| Crawford (2004) GJI 157 | yes | Eq. 1: eta = k_w u/tau_zz. Eq. 5: vertical half-space compliance (lambda+2mu)/(2mu(lambda+mu)). "Density has almost no effect on the seafloor compliance", Vp less than Vs, less still for Poisson > 0.4. Sensitive at depths between lambda/2pi and lambda/4. Eq. 7-8: eta = k gamma |S_a|/(omega^2 |S_p|); error = sqrt(1-gamma^2)/(|gamma| sqrt(2 n_d)) |eta|, n_d data windows. Minimum structure: layers 0.05 km growing by 1.1, second differences of Vs, target misfit; Bayesian: few layers, shear modulus and interface depths, prior sigmas (Jackson & Matsu'ura 1985). Full elastodynamic equations always solved. |
| Zha (2014 thesis; = Zha et al. 2014, Zha & Webb 2016) | yes, ch. 2 and 4 | Gravity: a_w = 2 pi G rho_w exp(-kH) h_w (Crawford 1998 eq. 6 without an errant 2 pi); after removal the actual compliance is slightly HIGHER than apparent below 5 mHz; free-air term negligible above 1 mHz. Inversion: DRAM (Haario 2006), likelihood exp(-(chi2 + alpha R)/2); sediment thickness 0-2000 m in 50 m steps with fixed Vs(z) from Ruan et al. (2014) and Vp(z) from Hamilton (1979), density 2.0; crust layers 100 m doubling; Moho fixed from refraction (7-15 km); mantle Vs 3.8-4.2; Vp/Vs 1.85 in crust and mantle; rho = 1.85 + 0.165 Vp (Christensen & Shaw 1970); 10 chains x 10000; keep chi2 < max(2 chi2_min, chi2_min + 0.25) (Shen et al. 2013). Sediment Vp/Vs up to 15 at the seafloor, 2.5 at 1 km. Hamilton instead of Ruan changes sediment thickness by about 200 m at thick-sediment sites. DPG gains from Rayleigh waves 1.3-2 (one 0.2). |
| Hulme et al. (2003, 2005), Crawford et al. (2002), Crawford (2005) | on disk, sections consulted as needed | 2-D effects. |
| Mosher et al. (2021) GJI 227 | sections | Bernstein-polynomial Vs(z), monotonic, constant Vp 6 km/s and constant density; mixture density network. |
| Janiszewski et al. (2019) GJI 217 | sections | Sediment density 2.0 g/cc "following Ruan et al. (2014)". |
| Aminian et al. (2025), Aminian (2024) thesis | yes (earlier audits) | The method under test. |

NOT obtained (publisher and host block automated download; closed access): Crawford et al. (1991, 1998) JGR; Ruan et al. (2014) JGR; Bell et al. (2015) SRL; Doran & Laske (2019); Zheng et al. (2024). Consequence: the Ruan sediment relation is NOT used (a search snippet of a later paper gave a = 0.15608, b = 1.2198, c = 0.49473 without the functional form: unverified, do not use). The 1998 gravity equation is taken from Zha's quotation of it and re-derived below.

## Physics checks done

1. **Gravity terms (data side), re-derived.** z up; under a crest p > 0 the seafloor goes down, u = -xi p. A seismometer reads specific force: a_meas = d2u/dt2 - dg, with dg the change of the upward gravitational acceleration. Wave crest: dg = +a_w = (2 pi G/g) cosh(kH) exp(-kH) p. Moving down by |u| increases gravity: dg = (2g/R) u. So a_meas/p = (omega^2 + 2g/R) xi - a_w/p, hence xi = (omega^2 |Z/P| + a_w/p)/(omega^2 + 2g/R). `gravity_corrected_compliance` implements exactly this. The correction RAISES compliance (as Zha says). The GJI paper's sentence that a_w "can increase the measured compliance by up to 5 per cent" has the sign and the size wrong: at RR52, 7.3 mHz, a_w/p is about 19 % of omega^2 |Z/P|. The neglected Bouguer term of the displaced seafloor, 2 pi G (rho_rock - rho_w) u, is about 20 % of the free-air term, itself 0.3 % at 5 mHz: negligible.
2. **What compliance constrains.** Quasi-static half-space: eta = (lambda+2mu)/(2mu(lambda+mu)) = (1 - nu)/mu. The data constrain M = mu/(1 - nu) = 2 rho Vs^2 (1 - Vs^2/Vp^2). Vs needs rho and Vp/Vs from elsewhere. Between Vp/Vs 1.73 and 2.4 the factor r^2/(r^2-1) goes from 1.50 to 1.21: up to 11 % in Vs for the same data. A 10 % density error is 5 % in Vs. This is why earlier results moved with the Vp rule.

## Code reading (part A) findings

- `compy/compliance.py`, `numerics.py`: definitions consistent (k from full dispersion relation, g = 9.79329, Hann 50 % Welch with 0.947 dof factor; error formula = Crawford 2004 eq. 8 with effective averages).
- Stage C: 35-min windows, nseg 1024 at 1.25 Hz (df 1.22 mHz), 4 Welch segments per window. Adjacent Hann bins are correlated (not in the covariance).
- `dram.py`: DRAM acceptance, Gibbs smoothing step and tempering swaps read as correct. Smoothing acts on second differences in LAYER INDEX (layers double in thickness), i.e. roughly in log depth.
- **Error model (lead confirmed on RR52):** the covariance is random + monthly/n (shrunk 50 %) + 1 % floor as INDEPENDENT noise per bin; sampled gain. At RR52 the marginal sigma is 1.1-1.3 %, almost all from the floor, correlation between bins about 0.1. The actual bin-to-bin noise of the stack is about 0.3 %. Real systematic errors (response, estimator bias, 2-D structure) are smooth in frequency, not white. So white 1 % is too large for scatter (chi2/N 0.39) and too strict against smooth distortions of the curve.
- RR52: 9 data, 10 parameters + smoothing weight. The problem is prior-dominated at depth.

## E1 forward checks (script `scripts/e01_forward_checks.py`, output `14_inversion_overhaul/e01_forward/`)

- **Verified by test.** A global-matrix solver written for this audit (50-digit mpmath, all layer amplitudes in one linear system) agrees with `dram.Forward` to 3.6e-7 and with `impedance` to 2.8e-8 over 60 random models (half with slow sediment, water 150-5000 m). Sign positive. The forward model is correct.
- Half-space: eta/((1-nu)/mu) = 1.001-1.003 for Vs 3 km/s, 1.005-1.012 for 1.5 km/s, 1.16-1.51 for Vs 300 m/s (inertial terms matter for slow half-spaces, as Crawford 2004 says).
- Gradient layer split in n sublayers: error 16 % (2), 4.9 % (4), 1.3 % (8), 0.34 % (16), 0.09 % (32): second-order convergence. A gradient needs about 16 sublayers for 0.3 %.
- **Sensitivities (log derivatives), bare rock under 3000 m:** summed over layers d ln eta/d ln Vs = -0.80 to -0.96, d ln eta/d ln Vp = -1.02 to -1.17, d ln eta/d ln rho = -0.99 (exact scaling). In RELATIVE terms compliance is as sensitive to Vp as to Vs. The literature's "2 to 9 times more sensitive to Vs" (Crawford 1994 eq. 2.21) is per m/s, i.e. the same thing divided by Vp/Vs. Half-space: d ln eta/d ln Vs = -2 + 2/(r^2-1), d ln eta/d ln Vp = -2/(r^2-1), r = Vp/Vs.
- **Thin soft layers act through their P modulus.** For kh << 1 a layer of thickness h bonded to stiffer rock is in uniaxial strain and adds eta = k h/(rho Vp^2), independent of Vs. Check with the mpmath solver, 3000 m water, 100 m (Vp 1560, Vs 200) + 200 m (Vp 1750, Vs 450) of sediment on rock: the sediment adds 58 % (4 mHz) to 80 % (6 mHz) and 133 % (10 mHz) to the bare-rock compliance; added/(k h/(rho Vp^2)) = 0.97, 0.98, 1.03 at 4, 6, 10 mHz (1.25 at 15 mHz, 1.76 at 20 mHz where shear enters). Halving the sediment Vs changes eta by 1 %, 3 %, 16 % at 4, 6, 10 mHz. So at the low frequencies that carry the deep-crust information, sediment enters only through thickness/(rho Vp^2); the sediment Vs matters where kh > 0.15-0.2. The same holds for a thin layer 2A.
- Consequences for the inversion: (a) sediment thickness and sediment P modulus trade directly with deep crustal Vs; (b) a fixed sediment Vp(z) and density are strong assumptions, not harmless ones; (c) a prior on sediment thickness from an independent map (GlobSed) is valuable; (d) statements that compliance is insensitive to Vp at high Vp/Vs hold for a half-space, not for a thin layer.

## E2 error model on the archive (script `scripts/e02_error_model.py`, output `14_inversion_overhaul/e02_errors/error_model.csv`)

285 of 287 stations (XO.LT02, LT15 have no band). Restack equals the archive to 1.5e-14. Tier 1 (172 stations), medians [10 %-90 %], all relative to compliance:

| Quantity | Value |
|---|---|
| Bendat & Piersol error stored by stage C | 0.15 % [0.09-0.38] |
| Bootstrap over windows (independent) | 0.23 % [0.10-0.55]; 1.36 x B&P |
| Correlation of neighbouring bins (Hann, 50 % overlap) | 0.48; lag 2: 0.11 |
| Bootstrap over 30-day blocks | 0.38 % [0.14-2.2]; 1.7 x independent; neighbouring-bin correlation 0.80 |
| Share of the 30-day covariance that is a common level / a tilt in log f | 59 % / 11 % |
| First half against second half of the record | level 0.67 %, tilt 0.83 %, rms 0.93 % [0.36-5.9] |
| H2/H1 = 1/coherence^2 of the stack | 1.09 band median; 1.15 at the low edge, 1.20 at the high edge |
| eta of the loudest fifth of windows over the quietest fifth | 0.985 [0.88-1.003] |

Reading:
- The old covariance gave about 1.2 % INDEPENDENT error per bin (1 % floor). Real scatter is 0.2-0.4 % and strongly correlated between bins. That alone explains chi2/N = 0.39.
- What is uncertain is the LEVEL (gain) and the SLOPE of the curve, not its bin-to-bin shape. The old model was too strict on slope errors and too loose on wiggles.
- Pressure-channel noise would make H1 rise with the pressure level. It falls slightly instead (-1.5 % median), so pressure noise is not what limits coherence, H1 is the right estimator and H2 (+9 %) is not a bracket worth carrying. The small fall with level is in the direction of the coherence-gate selection bias found in the 2.9 audit.
- Two independent calibrations (stage D Rayleigh against BRUIT-FM microseism, 36 stations): robust scatter of the log ratio 2.0 %; all 73 pairs: 2.7 % robust, 16-84 % range -5.4 to +2.6 %, with outliers. A Gaussian gain prior of 5 % (1 sigma) is supported; uniform +-5 % is too narrow in the tails.

## Design decisions taken from E1-E2 (to be tested, D49-D5x)

1. Covariance = bootstrap covariance of the stack (independent windows, with its bin correlation) + slow level and tilt variance measured from 30-day blocks + a tilt term for the pressure-response shape + a smooth (correlated in log f) model-error term. No white floor. Gain sampled with a Gaussian prior.
2. Prior: ln Vs of each crustal layer Gaussian about a reference model (CRUST1.0 or a user model), sigma 0.3, correlation 0.7 between neighbouring layers; proper, so unresolved depths return to the reference with the prior width instead of extrapolating a gradient to the bounds. The old flat prior with hierarchical second-difference smoothing stays as a sensitivity check.
3. Vp follows Vs by a rock relation with a sampled scatter (upper and lower crust separately); density follows Vp (Carlson & Raskin 1984). Sediment: thickness, Vs scale and P-modulus factor free.
4. Brocher's regression (eq. 9) excludes mafic rocks, gabbros and serpentinites by its own scope (secondary source: R package hydrogeo, which lists eq. 6 as "all except Ca-rich and mafic, gabbros and serpentinites" and eq. 8, the mafic line Vs = 2.88 + 0.52 (Vp - 5.25) for 5.25 <= Vp <= 7.25 km/s). The original paper was not obtained. For oceanic crust use the mafic line (Vp/Vs 1.82-1.85), which agrees with Zha & Webb's 1.85.

## New code (10 Nov, afternoon)

- `compy/inversion/errors.py`: `stack_covariance` (bootstrap over windows + slow level/tilt from 30-day blocks), `default_covariance` (same structure from plain errors), `usable_bins` (trim noisy band edges), `compliance_from_sums`, `attraction_ratio`.
- `compy/inversion/auto.py`: `SeafloorModel` (soft sediment unit: thickness, Vs scale, P-modulus factor; crust grid 100 m doubling with ln Vs Gaussian prior about the reference, sigma 0.3 (0.5 where CRUST1.0 sediment >= 500 m, whose layers then join the reference), correlation 0.7; free Moho depth cutting the grid, Gaussian prior about the reference, 30 %; Vp = `oceanic_vp`(Vs) x sampled factor (upper, lower crust, sigma 7 %); density = `crust_density`(Vp) x sampled factor (3 %); mantle Vs uniform 3.9-4.9 with the Moho rule), `compliance_data` (covariance = bootstrap + tilt + smooth 0.5 %; no white floor; gain not in covariance), `invert` (chains from prior draws, PT ladder 1.6^0..5, rounds doubling until R-hat <= 1.05 and ESS >= 400; exact prior draws), `InversionResult` (profile with resolved flag, `parameter_resolution`, quality label, save).
- `dram.py`: `ComplianceData.build` gained `random_cov`, `tilt_rel`, `smooth_rel`; `_Target.log_prior` (model prior + Gaussian gain prior `gain_sigma`); `_Dram` carries the untempered log prior. Defaults unchanged, old tests pass.
- Tests: `tests/test_inversion_auto.py` (17 tests).
- Archive scripts in `scripts/`: `archive.py` (inputs: windows from the 13_audit cache, tide-based response correction where the 2.9 audit has a leak factor, tilt = hypot(0.02, half the tide-phase difference)), `r01_run_station.py`, `r02_batch.py`, `r03_collect.py`, `e03_sbc.py`.

## Pilot 1 (t02_pilot, 15 stations; before the band trim and the thick-sediment rule)

- 90 s per station on 6 cores when one round suffices. YV: RR52, RR28, RR29, RR34, RR36, RR40 good; RR38 (chi2 24.5/6) and RR50 (15.3/6) misfit. 7B.NN04 good, 7D.J44A good (9.5/10), XE.CC05 marginal (16.5/8), 7D.J35A gain pulled 3.4 sigma (0.83), 7K.AN3 fits (11.8/13) but not converged after 3 rounds (R-hat 1.06), 7D.J09B fails (276/33).
- chi2_min for N = 5-9 data: 0.3-9.8, median about N-2: the new covariance is about right; no white floor needed at these stations.
- **Per-parameter resolution (posterior/prior width), deep stations:** sediment thickness 0.02-0.24 (well determined, given the sediment P modulus); crust layers 1-3 (top 700 m) 0.55-0.98 (hardly resolved one by one); layers 4-6 (0.7-6.3 km) 0.16-0.7; mantle 0.35-1.0; Moho depth 0.5-0.9; Vp factors, density, gain about 1.0 (pure prior, as they should be).
- **Compared with the 2.8 run:** new Vs is 10-30 % lower at 500 m at five YV stations and 3-23 % higher at 2-5 km at most. At 2-5 km the YV medians are 3.3-4.6 km/s (RR38 4.58, RR40 4.0-4.5, RR50 4.1-4.4, RR52 4.1, RR36 4.4 at 5 km): stiffer than gabbro. About +5 % comes from Vp/Vs (2.8 run posterior Vp/Vs 2.1 for all layers against 1.84 +- 7 % now) and +4.5 % from density (Christensen-Shaw with that Vp gave 3.29 g/cc against 3.0). A 3 % slope error of the curve moves Vs at 2-5 km by about 6 %: the tilt allowance decides how firmly the deep Vs is known.
- 7D.J09B: 2 km of Hamilton sediment fits the curve to 10 % between 15 and 34 mHz; above 35 mHz the data collapse (bootstrap sd 1-6 %, eta of loud windows 20-28 % below quiet ones): pressure signal gone (kH > 1.4). The stage C coherence gate (0.8) does not catch it; `usable_bins` does.
- 7K.AN3 (CRUST1.0: 6.5 km of sediment): the posterior used the "crust" layers at 450-1700 m/s, 3-6 sigma from an igneous reference. Hence the thick-sediment rule.

## E3/E4 calibration (144 synthetic data sets each; `scripts/e03_sbc.py`, `e04_sbc_summary.py`; output `14_inversion_overhaul/e03_sbc/`, `summary.json`)

Expected coverage for 144 replicates (95 % range): 0.42-0.58 for 50 % intervals, 0.85-0.94 for 90 % intervals.

- **Truth from the prior (calibration proper), templates YV.RR52 (deep, 9 data), 7B.NN04 (abyssal, 8), 7K.AN3 (slope, thick sediment, 13), 7D.J09B (shelf, 23):** 67 parameter rows; 50 % coverage 0.42-0.60, 90 % coverage 0.83-0.95; five rows marginally outside (density 0.83-0.84 twice, vs_mantle 0.59, vs_crust_6 0.60, sediment_vs_scale 0.95), as expected by chance for this many rows. Vs at eight depths: 50 %: 0.39-0.60, 90 %: 0.83-0.94; median bias within 3 %. Median chi2/N 0.82-0.95. The inversion is calibrated.
- **Truth the model cannot represent (smooth gradients on 60 layers, Vp/Vs 1.76-2.5 varying with depth, Christensen-Shaw density, own sediment and Moho), RR52 and NN04:** below 700 m the 90 % intervals contain the truth 97-100 % of the time and the 50 % intervals 51-79 % (conservative), median error 2-9 %, bias -5 to +8 %. In the top 300 m Vs is overestimated by 16-26 % (median) and 90 % coverage falls to 0.66-0.83 at 50-150 m: the truths have a slow layer 2A (0.9-2.6 km/s at the seafloor), the CRUST1.0 reference has 2.54 km/s, and the data do not resolve those depths. This is the regime the result flags as unresolved.
- **Unmodelled tilt of 4 % (twice the allowance for a corrected response), RR52:** Vs coverage stays 0.82-0.92 (90 %); bias +2 % at 3 km, +4 % at 5 km, +7 % at 7 km.
- **Gain error of 10 % (twice the prior sigma), RR52:** Vs 90 % coverage 0.78-0.88; the gain's own coverage collapses (0.56), as it must.
- **Legacy identity (C1):** variant `final` of ComPy 2.8 on RR52 reproduces the archived draws exactly (max difference 0.0) with the ComPy 3 sampler code.

## Sensitivity runs (43 benchmark stations; `scripts/r02_batch.py`; output `14_inversion_overhaul/r02_sensitivity/`)

Change of the median Vs against the standard run, stations where the choice differs and neither run is poor.

- **Prior width of ln Vs 0.5 instead of 0.3 (17 stations; the 24 with thick reference sediment already use 0.5 and repeat the standard run exactly):** median change below 1 % at every depth, largest 11 % (5 km); intervals 1.2-1.3 times wider down to 3 km.
- **Layer correlation 0.9 instead of 0.7 (38):** median change at most 1.5 %; intervals 0.7-0.9 of their width; single stations up to 34 % at 1 km.
- **No response correction, tilt 4 % (34):** median change below 1 %; 16-84 % range -7 to +2 % at 5 km; six labels change (four worse, two better).
- **Vp factor 15 % instead of 7 % (38):** median change below 1 %; intervals 1.1 times wider from 2 km down.
- **Band cut at 0.6 f_c instead of 0.8 (35; two stations left with fewer than 3 bins):** three data fewer (median); median change +3.5, +4.2, -0.6, +6.4, -0.9 % at 0.5, 1, 2, 3, 5 km; 20 % or more somewhere between 0.5 and 3 km at 14 stations (7K.AN2-AN6, 7D.M11B, ZU.NZ19 largest); new median outside the standard 16-84 % interval at 14-37 % of stations.
- **The "resolved" label depends on the prior width:** thin-sediment stations resolved at 1 km: 7 of 19 with 0.3, 15 with 0.5. With correlation 0.9, 0.5 km resolved at 31 of 43 instead of 22.
