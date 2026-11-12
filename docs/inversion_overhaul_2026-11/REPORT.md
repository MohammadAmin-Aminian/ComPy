# ComPy inversion overhaul, November 2026

Audit and rebuild of the shear-velocity inversion of ComPy, against the primary literature and against tests that do not use any earlier ComPy result. Branch `inversion-overhaul`, proposed as ComPy 3.0.0. Brief: [PROMPT.md](PROMPT.md). Working notes: [NOTES.md](NOTES.md).

Every statement below is marked **[test]** (verified by a script in `scripts/`, output in `~/Desktop/ComPy/14_inversion_overhaul/`), **[source]** (taken from a publication that was read) or **[inferred]**.

## 1. Summary

**What was wrong in the inversion of ComPy 2.3 to 2.8**

1. **The error model.** Each frequency bin carried an independent 1 % error. The real scatter of a stacked compliance curve is 0.2 to 0.4 % and strongly correlated between bins; what is uncertain is the level of the curve (pressure gain) and its slope (pressure response at low frequency). The old model was five times too loose on the shape of the curve, which explains chi-square/N near 0.4, and too strict on its slope, which is what carries the information about the deep crust. [test]
2. **The prior.** Layers had a flat prior with a penalty on second differences. In the layers the data do not see, that penalty continues the shallow gradient and drives Vs to its bounds. This is the "deep layers press against their limits" of the earlier reports. [test]
3. **Vp and density.** Compliance measures stiffness, not Vs. One free Vp/Vs for all crustal layers settled near 2.1, which gives lower-crust Vp above 8.5 km/s and, with the density rule tied to Vp, a density of 3.3 g/cc. Both lower Vs: by about 10 % at 2 to 5 km. [test]
4. **Sediment.** A soft layer much thinner than the wavelength adds compliance through its P modulus, `k h / (rho Vp^2)`, whatever its Vs. The model fixed sediment density at 2.0 g/cc (soft sediment is nearer 1.6) and its Vp, and limited thickness to 1 km. Sediment thickness, and through it deep Vs, inherited those choices. [test]
5. **Fixed Moho** from a 1-degree model, and **uniform gain bounds** of 5 %, where two independent calibrations of the same gauges differ by more than that in one case out of six. [test]
6. **Band.** At shelf stations the estimate collapses above the frequency where the pressure signal still reaches the seafloor; the coherence gate does not detect it, and the saved curve was also cut at 40 mHz. Shelf station 7D.J09B never fitted for these reasons. [test]

**What was right.** The forward model (agreement with an independent 50-digit solution to 4e-7), the gravity corrections of the data (re-derived; the 2025 paper states their sign and size wrongly, the code is right), the sampler (acceptance rules read and tested; a ComPy 2.8 run repeats draw for draw), and the H1 estimator. [test]

**What changed.** A new default inversion, `compy.inversion.auto`, with errors measured from each station's own windows, a proper prior about a reference model, free Moho depth, Vp and density uncertain, sediment with its P modulus, a Gaussian gain prior, automatic band limits, sampling that repeats until it converges, and a result that marks the depths the data resolve. One command: `compy-invert`. Method: [../INVERSION.md](../INVERSION.md).

**Is it calibrated?** On synthetic data drawn from the prior, on four station types, the 50 % and 90 % intervals contain the truth at those rates (section 5). When the true Earth is one the model cannot represent, intervals stay honest below 700 m and Vs in the top 300 m follows the reference, not the data; the result flags those depths as unresolved. [test]

**What the data resolve.** At a deep-ocean station: the thickness of soft sediment, and broad averages of Vs between about 0.7 and 6 km. Not the top few hundred metres of basement, not thin layers, and nothing about Vp/Vs or density. [test]

RESULTS_SUMMARY_PLACEHOLDER

## 2. Literature

| Source | Read | What was taken |
|---|---|---|
| Crawford (1994), PhD thesis, UCSD | ch. 2 and 6 (ch. 3 and 4 are placeholders in the PDF) | Half-space compliance, eq. 2.19; sensitivity ratios, eqs 2.20 to 2.22; Vp fixed from refraction; start Vs = Vp/1.87; density of Carlson & Raskin (1984), eq. 6.1; Occam inversion; sensitivity near lambda/6 |
| Crawford (2004), GJI 157 | in full | Definition, eq. 1; `eta = (lambda + 2 mu)/(2 mu (lambda + mu))`, eq. 5; "density has almost no effect"; depth between lambda/2pi and lambda/4; error formula, eq. 8; minimum-structure and Bayesian inversions; full equations of motion always solved |
| Zha (2014), PhD thesis, Columbia (holds Zha et al. 2014 and Zha & Webb 2016) | ch. 2 and 4 | Gravity correction as Crawford et al. (1998) eq. 6 without an errant 2 pi, sign of its effect; DRAM; sediment 0 to 2000 m with fixed Vs(z) (Ruan) and Vp(z) (Hamilton), density 2.0; crust layers 100 m doubling; Vp/Vs 1.85; density 1.85 + 0.165 Vp; mantle 3.8 to 4.2 km/s; Moho from refraction |
| Hulme et al. (2003, 2005); Crawford et al. (2002); Crawford (2005) | sections | 2-D effects, as a limit |
| Mosher et al. (2021), GJI 227 | sections | Another parameterization (Bernstein polynomials, constant Vp and density) |
| Janiszewski et al. (2019), GJI 217 | one statement | Sediment density 2.0 "following Ruan et al. (2014)" |
| Tsekhmistrenko et al. (2026), Seismica 5(2) | sections | Instruments of network 8J (for the 8J paper) |
| Aminian et al. (2025), GJI 242; Aminian (2024) thesis | earlier audits | The method under test |

**Not obtained** (publisher or host refused automated access, or closed): Crawford et al. (1991, 1998); Ruan et al. (2014); Bell et al. (2015); Doran & Laske (2019); Zheng et al. (2024); Brocher (2005), Hamilton (1976, 1979), Carlson & Raskin (1984), Christensen & Shaw (1970) in the original. Consequences: the Ruan sediment relation is not used; Brocher's relations are taken from a published implementation (R package hydrogeo) and marked so in the code; the 1998 gravity formula is taken from Zha's quotation and re-derived (section 3.1). Please check these against the originals.

## 3. Physics checks

### 3.1 Gravity terms of the data

With z up, a pressure p pushes the seafloor down, `u = -xi p`. A seismometer reads specific force, `a = d2u/dt2 - dg`. A wave crest attracts the mass upward, `dg = a_w = (2 pi G / g) cosh(kH) exp(-kH) p`; moving down by |u| raises gravity by the free-air gradient, `dg = (2g/R) u`. Hence `a/p = (omega^2 + 2g/R) xi - a_w/p`, and

`eta = k (omega^2 |Z/P| + a_w/p) / (omega^2 + 2g/R)`,

which is what `compy.compliance.gravity_corrected_compliance` computes. The correction raises compliance, as Zha et al. (2014) say. At YV.RR52 it is 18 % of compliance at 7.3 mHz and 6 % at 17 mHz. The 2025 paper says the attraction "can increase the measured compliance by up to 5 per cent": sign and size are wrong in the text, not in the code. The attraction of the displaced seafloor itself is a fifth of the free-air term, which is 0.3 % at 5 mHz: negligible. [inferred from first principles; formula of a_w from source]

### 3.2 Forward model

A global-matrix solver written for this audit (all layer amplitudes in one linear system, 50-digit arithmetic) agrees with `dram.Forward` to 3.6e-7 and with `impedance` to 2.8e-8 over 60 random models, half of them with slow sediment, under 150 to 5000 m of water. A gradient split into n sublayers converges as 1/n^2 (16 sublayers: 0.3 %). [test: `e01_forward_checks.py`]

### 3.3 What compliance is sensitive to

![Figure 1](figures/fig01_physics.png)

*Figure 1. Left: compliance of an oceanic crust under 3000 m of water without sediment (black) and with 300 m of sediment (blue); the dashed line adds only the uniaxial term k h/(rho Vp^2) to the bare rock; the band is the range when the sediment Vs is halved or doubled. Below 10 mHz the sediment acts through its P modulus alone. Right: change of compliance for a 1 % change of Vs, Vp or density in all layers. Read it as: compliance is about as sensitive to Vp as to Vs.*

- Half-space: `eta = (1 - nu)/mu`. Between Vp/Vs 1.73 and 2.4 the same compliance gives Vs 11 % apart. [source: Crawford 2004 eq. 5; arithmetic]
- Summed over a normal crust, `d ln eta/d ln Vs` is -0.80 to -0.96, `d ln eta/d ln Vp` -1.02 to -1.17, `d ln eta/d ln rho` -0.99. The published "2 to 9 times more sensitive to Vs" is per m/s. [test]
- Thin soft layer: 300 m of sediment adds 58 % (4 mHz) to 133 % (10 mHz) to bare-rock compliance; the added part equals `k h/(rho Vp^2)` to 3 % up to 10 mHz; halving the sediment Vs changes compliance by 1 % at 4 mHz and 3 % at 6 mHz. [test, both forward codes and the 50-digit solver]
- Attenuation: seismic velocities at 5 to 10 Hz exceed those at 10 mHz by `ln(700)/(pi Q)`: 0.4 % for Q = 500, 4 % for Q = 50, 10 % for Q = 20. Compliance Vs of sediment and of cracked upper crust should be read a few per cent below refraction Vs for this reason alone. [inferred: standard constant-Q dispersion, not re-read]

## 4. Errors of a compliance curve

![Figure 2](figures/fig02_errors.png)

*Figure 2. 172 stations of the archive. Left: share of stations whose stacked curve has at most the given relative error, by three measures; the dotted line is the 1 % floor of ComPy 2.x. Middle: correlation of the errors between frequency bins. Right: change of the level and of the tilt of the curve between the two halves of the record.*

| Quantity (172 stations) | Median [10 % to 90 %] |
|---|---|
| Bendat & Piersol error | 0.15 % [0.09 to 0.38] |
| Bootstrap over windows | 0.23 % [0.10 to 0.55]; 1.36 times Bendat & Piersol |
| Correlation of neighbouring bins; two apart | 0.48; 0.11 |
| Bootstrap over 30-day blocks | 0.38 % [0.14 to 2.2] |
| Share of the 30-day covariance that is a level / a tilt | 59 % / 11 % |
| First half against second half: level; tilt | 0.67 %; 0.83 % |
| H2/H1 of the stack | 1.09 |
| Compliance of the loudest fifth of windows over the quietest fifth | 0.985 [0.88 to 1.003] |

[test: `e02_error_model.py`; the restack equals the archive to 1.5e-14]

- **Estimator.** Noise on the pressure channel would make H1 rise with the pressure level. It falls slightly. So coherence is limited by noise on the vertical, H1 is the right estimator, and H2 (9 % higher) is not a useful bracket. The small fall is in the direction of the selection bias of the coherence gate found in the 2.9 audit. [test]
- **Gain.** The Rayleigh-wave gain and the microseism gain of 73 gauges differ by 2.7 % (robust scatter; 2.0 % for the 36 microseism pairs), with a 16 to 84 % range of -5.4 to +2.6 % and outliers beyond 25 %. A Gaussian prior of 5 % is supported; a hard +-5 % is not. [test]
- **Slope.** The tide-based correction of the pressure-gauge leak corner (2.9 audit) changes compliance by -3 % (median) at the bottom of the band, -6 % at YV.RR52. The allowance for a slope error is 2 % where that correction is applied, more where the tide and the phase estimate of it disagree, 4 % where no correction exists. This number is a judgement; section 5 shows what happens when it is exceeded. [inferred]

![Figure 3](figures/fig03_covariance_gain.png)

*Figure 3. Left and middle: the errors given to the inversion at YV.RR52 by ComPy 2.8 and by ComPy 3: their size, and their correlation between bins. Right: difference between two independent pressure calibrations; the curve is what two gains with 5 % error each would give.*

## 5. Calibration of the new inversion

Synthetic data sets, 144 per case: a true model is drawn, its compliance computed, noise added with the station's covariance, and the result inverted (`e03_sbc.py`). For 144 cases a correct 50 % interval contains the truth in 42 to 58 % of them and a 90 % interval in 85 to 94 %.

| Case | Parameters within range | Vs at 8 depths, 50 % | Vs, 90 % | Median bias of Vs |
|---|---|---|---|---|
| Truth from the prior, YV.RR52 (deep, 9 data) | 15 of 16 | 0.39 to 0.60 | 0.83 to 0.92 | within 1.1 % |
| Truth from the prior, 7B.NN04 (abyssal, 8 data) | 15 of 16 | 0.45 to 0.54 | 0.85 to 0.91 | within 2.4 % |
| Truth from the prior, 7K.AN3 (slope, thick sediment, 13 data) | 15 of 17 | 0.46 to 0.55 | 0.88 to 0.94 | within 1.4 % |
| Truth from the prior, 7D.J09B (shelf, 23 data) | 17 of 18 | 0.42 to 0.55 | 0.90 to 0.94 | within 2.9 % |
| Smooth-gradient Earth the layers cannot represent, RR52 | n/a | 0.47 to 0.78 below 300 m | 0.94 to 1.00 below 300 m; 0.66 to 0.78 above 300 m | -3 to +2 % below 700 m; +17 to +25 % above 300 m |
| The same, NN04 | n/a | 0.50 to 0.79 below 300 m | 0.91 to 1.00 below 300 m; 0.72 to 0.83 above 300 m | -5 to +8 % below 700 m; +16 to +26 % above 300 m |
| Unmodelled slope error of 4 %, RR52 | n/a | 0.37 to 0.59 | 0.82 to 0.92 | +2 % at 3 km, +4 % at 5 km, +7 % at 7 km |
| Gain error of 10 % (twice the prior), RR52 | n/a | 0.33 to 0.50 | 0.78 to 0.88 | within 1.4 % |

[test; full numbers in `e03_sbc/summary.json`]

Reading:

- With the truth drawn from the prior, coverage is nominal for every parameter type on all four station types. The few rows outside the range (five of 67) are what chance gives.
- With an Earth of smooth gradients, a depth-dependent Vp/Vs between 1.76 and 2.5 and another density rule, intervals below 700 m are conservative and the median error of Vs is 2 to 9 %. In the top 300 m the truth has a slow layer 2A, the reference does not, and the data cannot tell: Vs there is the reference. These are the depths the result marks as unresolved.
- A slope error twice the allowance biases Vs below 3 km by 2 to 7 %. This is the main systematic risk of deep Vs.

## 6. The earlier sampler still repeats

Variant `final` of ComPy 2.8 at YV.RR52, run with the ComPy 3 code, gives the archived draws exactly (maximum difference 0.0 over 8 x 5000 x 10 values). The additions to the sampler do not change a 2.x run. [test: `c01_legacy_identity.py`]

ARCHIVE_RESULTS_PLACEHOLDER
