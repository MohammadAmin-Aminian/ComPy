# ComPy checked against the papers (November 2026)

This note compares each stage of the ComPy pipeline with the papers it rests on. It was
written before ComPy 3.0.0 was merged. Several of these papers could not be obtained for
the inversion overhaul ([report](../inversion_overhaul_2026-11/REPORT.md), section 2); they
were read here in the original.

**Result.** No formula in the pipeline had to change. Three small things were added, none
of which changes a default result (section 3). Section 4 lists where ComPy differs from a
paper and why, and section 5 what is still not verified.

## 1. What was read

| Paper | Read | Used for |
|---|---|---|
| Bell, Forsyth & Ruan (2015), BSSA 105, 300-313 | in full | transfer functions and their errors, tilt removal, band limit |
| An, Cai, Zheng, Meng & Liu (2017), GRL 44, 10,272-10,281 | in full | pressure calibration |
| Doran & Laske (2021), GJI 224, 843-857 | in full | deep-water band, calibration, inversion choices |
| Mosher, Audet & Gosselin (2021), G-cubed 22, e2021GC009720 | in full | compliance measurement, tests of assumptions |
| Crawford, Webb & Hildebrand (1991), JGR 96, 16,151-16,160 | the page with eqs 3, 4 and 7, and the data-selection statement | estimator, error, half-space |
| Crawford, Webb & Hildebrand (1998), JGR 103, 9895-9916 | section 2 (eqs 1-10) | gravity terms |
| Zha & Webb (2016), JGR 121, 3220-3237 | sections 2-3 and appendix A | selection, calibration, parameterization |
| Deng, An, Cai & Ren (2022), GRL 49, e2021GL096952 | sections 1, 2 and 4 | calibration with compressible water |
| Ruan, Forsyth & Bell (2014), JGR 119, 6357-6371 | sediment model (eq. 5, Table 2) | sediment Vs |
| Webb et al. (2026), calibration of pressure sensors | the abstract and the passages on the capillary leak | pressure-gauge response |
| Crawford & Webb (2000), BSSA 90, 952-963 | the statements on horizontal compliance | tilt removal |
| Mosher, Eilon, Janiszewski & Audet (2021), GJI 227, 1879-1892 | the model assumptions only | inversion |

## 2. Stage by stage

### Tilt removal (stage B)

| Paper | ComPy | |
|---|---|---|
| Remove from the vertical what is coherent with the horizontals, by transfer functions (Crawford & Webb 2000; Bell et al. 2015) | `compy.preprocessing.remove_coherent_noise`, both horizontals, conditioned spectra | same |
| Tilt changes through a deployment and must be estimated at least daily (Bell et al., Fig. 8) | hourly tilt, median over 24 h, applied per sample; transfer functions per window of at most 24 h | same |
| Transfer functions from windows without earthquakes (Bell et al.) | catalogue spans and STA/LTA spans are left out of the estimate | same |
| The horizontals carry no usable compliance signal: it is 2-4 times smaller than on the vertical and below the horizontal noise (Crawford & Webb 2000) | assumed; pressure is never a noise channel | same assumption |

### Compliance (stage C)

| Paper | ComPy | |
|---|---|---|
| Transfer function `G_xy/G_xx` with pressure as input (Bell et al., eq. 4); compliance is `k` times `\|Z/P\|` (Crawford et al. 1991, eq. 3; Mosher et al., eq. 4) | `k \|S_pz\|/S_pp` from stacked spectra | same |
| Error `sqrt(1 - g^2) / (\|g\| sqrt(2 n_d))` (Crawford et al. 1991, eq. 4; Bell et al., eq. 8) | `compy.numerics.compliance_uncertainty`, with `n_d` the effective number of Welch averages | same |
| `omega^2 = g k tanh(kH)` (Mosher et al., eq. 2); Mosher et al. solve it with a rational approximation of tanh (eq. 5, error up to 2.5 %) | solved to rounding (`compy.numerics.wavenumber`) | exact; at 4.4-4.8 km depth the approximation would be 1.5-2.6 % off in the band |
| Gravity terms (Crawford et al. 1998, eqs 4-6; Zha & Webb 2016, section 2.4) | `gravity_corrected_compliance` | same; see the note below |
| Keep data with pressure-vertical coherence above 0.8 (Zha & Webb 2016; Mosher et al.) | median coherence 0.9, or 0.8 if fewer than 500 windows pass | same or stricter; see the note on coherence |
| Compliance exists below about `sqrt(g / 2 pi H)` (Bell et al., eq. 11; their measured limit is about 10 % higher, Fig. 10b). At 4.4 km and more it is usable at periods above about 60 s (Doran & Laske, section 2.1); their curves span 70-140 s | band searched up to `f_c = sqrt(g / 2 pi H)`; the inversion keeps bins up to 0.8 `f_c` (66-70 s at the Young ORCA depths) | same |
| Pressure response limits reliable compliance to periods shorter than about 140 s (Doran & Laske, section 2.1) | the band may extend to 4 mHz (250 s) where coherence holds; the inversion carries a slope term for the response | **differs**, section 4 |

*Gravity terms.* The code was written from a quotation of Crawford et al. (1998). In the
original, eqs 4 to 6 give the deformation term, the free-air term and the wave attraction
as upward accelerations. A sensor reads the acceleration of its frame minus the change of
gravity, so the free-air term adds to `omega^2` and the attraction is added back to the
measured acceleration, as in the code. Eq. 6 of the paper prints the decay of the
attraction as `exp(-2 pi k H)` although `k` is already `2 pi / lambda`; the exponent is
`k H` (Zha & Webb 2016 write it so).

*Coherence.* Crawford et al. (1991, 1998) and ComPy use the amplitude coherence
`\|S_pz\| / sqrt(S_zz S_pp)`. Bell et al. (eq. 6) call its square the coherence. ComPy's
0.8 is therefore 0.64 in the convention of Bell et al. Zha & Webb's remark that a
coherence of 0.8 can bias compliance by up to 40 % if the noise is on the pressure shows
that their 0.8 is the amplitude coherence too.

### Pressure calibration (stage D)

| Paper | ComPy | |
|---|---|---|
| `p = rho h a` under a Rayleigh wave, for incompressible water and wavelengths long against the depth; verified at 0.02-0.2 Hz in 350-3300 m of water (An et al., eq. 6) | `water_ratio`: `rho h tan(x)/x` with `x = omega h sqrt(1/c_0^2 - 1/c^2)` (Zha & Webb 2016, A5; Deng et al. 2022, eq. 5), which reduces to `rho h` | the general form; needed at depth |
| One constant per gauge, by least squares against the prediction (An et al.) | one gain per event over bins with squared coherence of 0.95 or more, median of at least three events | same idea, more events |
| Events of magnitude 7 and more (Zha & Webb 2016) | magnitude 6 and more, shallower than 150 km, all eligible events | **differs**, section 4 |
| Scripps differential gauges read with reversed polarity (An et al., section 4.1) | pressure sign -1 by default, or from the Z/P phase | same |
| The gain may depend on frequency (Deng et al., Fig. 4g; Zha & Webb, appendix A; 3 % per decade in Webb et al. 2026) | one gain, with a 5 % prior and a slope term in the inversion | allowed for, not corrected |

At 4700 m `rho h` alone is 5 % too small at 0.02 Hz, 23 % at 0.04 Hz and 77 % at 0.06 Hz,
so the relation of An et al. cannot be used unchanged at the Young ORCA depths. The upper
end of the calibration band is then close to the first resonance of the water column
(`x` = 1.1 against `pi/2`), where the prediction depends on the Rayleigh phase velocity
(1 % for a 3 % change, at the top of the band) and on the sound speed (2 % for a 1 %
change).

### Inversion (stage E)

| Paper | ComPy 3 (`compy.inversion.auto`) | |
|---|---|---|
| Half-space compliance `Vp^2 / (2 rho Vs^2 (Vp^2 - Vs^2))` (Crawford et al. 1991, eq. 7; Doran & Laske, eq. 4; Mosher et al., eq. 3) | the propagator gives it to 0.1-0.3 % for a uniform model (the size of the inertial term the formula leaves out) | same |
| A linearised inversion holds for changes up to 10 % (Doran & Laske, Fig. 6) | the full forward model at every sample | not limited |
| Crust in layers of 100 m doubling with depth; sediment 0-2000 m; Hamilton (1979) Vp in sediment; DRAM sampler (Zha & Webb 2016) | the same | same |
| Vp/Vs fixed: 1.7 in the crust (Doran & Laske), 1.85 (Zha & Webb) | Vp from Vs by the mafic line (Vp/Vs 1.82-1.85) times a sampled factor with sigma 7 % | both inside the prior (1.7 is 1.2 sigma away) |
| Density `1.85 + 0.165 Vp` (Zha & Webb; Mosher et al., eq. 7) or not scaled (Doran & Laske) | Carlson & Raskin (1984) times a sampled factor, sigma 3 % | the two relations agree to 0.1 % at 6.5 km/s and 2.5 % at 4.5 km/s |
| Sediment Vs fixed to Ruan et al. (2014) (Zha & Webb) | Hamilton profile times a sampled scale (sigma 0.3) | the Ruan profile is 0.70-0.83 of Hamilton in the top 100 m, inside the prior; now available as `dram.ruan_vs` |
| Sediment density 2.0 g/cc (Zha & Webb; Ruan et al.) | Nafe-Drake curve, 1.64 g/cc at the seafloor | **differs**, section 4 |
| Uncertainty by Markov chains (Zha & Webb; Doran & Laske, section 3.2); report only what the data constrain (Mosher et al., section 3.1) | posterior intervals, prior shown, depths flagged as resolved or not | same intent |
| Test the result against the fixed assumptions (Mosher et al., section 4.1) | sensitivity runs in the overhaul report, section 8 | done there |

## 3. What was added

1. **Stage C records a check of the pressure response** (`response_check` in
   `<STA>_summary.json`). Bell et al. (Fig. 14) find that differential gauges are a few
   hundredths of a cycle out of phase with the vertical at 0.01 Hz after the nominal
   response is removed. Doran & Laske (section 2.2) say the corner of the gauge affects
   periods of 60 s and longer. Webb et al. (2026) measure leak time constants from 100 s
   in the laboratory to about 300 s at 5000 m and 2 °C. The Young ORCA metadata give 80 s
   (one pole at 2 mHz). Stage C now fits the factor on the metadata corner that explains
   the Z/P phase and writes it with the compliance bias it implies. Nothing is applied.
2. **Stage E has a variant `auto_phase`**, the default inversion on the curve corrected
   for that corner. It shows how much the response matters at a station. It is not the
   default because the fitted corner is a model: the whole-archive audit found that it
   predicts the tide only roughly.
3. **`dram.ruan_vs`**, the sediment profile of Ruan et al. (2014), as an alternative to
   the Hamilton profile.

The docstring of the gravity correction now states what the original paper says.

## 4. Differences left in place

- **Low-frequency limit of the band.** Doran & Laske stop at 140 s because of the pressure
  response. ComPy keeps lower frequencies where they are coherent and gives the inversion
  a slope term (4 % across the band). The response check of section 3 shows, station by
  station, whether that allowance is enough.
- **Calibration events.** Magnitude 6 and more instead of 7, with a strict coherence gate,
  which gives more events per station.
- **Calibration band 0.02-0.06 Hz at all depths.** Zha & Webb used 0.04-0.1 Hz at 2 km,
  which stays further from the water resonance (`x` up to 0.8) than 0.06 Hz does at 4.7 km
  (`x` = 1.1). A depth-dependent upper limit would reduce the dependence on the assumed
  phase velocity and sound speed. Not changed: the effect is about 1 % on the gain, inside
  its 5 % prior.
- **Sediment density.** Nafe-Drake instead of 2.0 g/cc, chosen in the overhaul because thin
  sediment acts through `rho Vp^2`; the P-modulus factor (sigma 15 %, bounds 0.6-1.6)
  spans the difference.
- **Vp/Vs and density are sampled, not fixed**, so their uncertainty is in the Vs
  intervals. The papers fix them.
- **The estimator assumes the incoherent noise is on the vertical.** If part of it is on
  the pressure, compliance is low by up to one minus the squared coherence (Zha & Webb
  2016, citing Crawford & Singh 2008). The `final_h2` variant brackets it.

## 5. Not verified

- Doran & Laske (2019), Brocher (2005), Hamilton (1976, 1979), Carlson & Raskin (1984)
  and Christensen & Shaw (1970) were not available in the original. The relations taken
  from them through other sources stay marked in the code.
- Whether the Z/P phase residual is the leak corner or something else (seismometer
  response, timing) cannot be told from one station; `compy.anomaly.phase_models` compares
  the hypotheses.
- Tilmann et al. (2020) on data uncertainty in Markov-chain inversions was not used.
