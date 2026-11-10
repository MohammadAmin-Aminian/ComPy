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
