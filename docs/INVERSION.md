# Shear-velocity inversion in ComPy 2.3

This note covers the inversion half of ComPy: what was wrong with the earlier
code, what `compy_inversion` does instead, how it was checked, and its limits.
The compliance half is covered in `docs/AUDIT_2026-10.md`.

## 1. Audit of the earlier inversion code

Three code generations were compared with the methods sections of Aminian et
al. (2025, GJI, §4.3), the thesis (Aminian 2024, §3.3) and the cited
literature. The generations are:

- v1, the code in the repository before 2.0 (commit `7400c4e`);
- v2, ComPy 2.0/2.1 (`inv_compy._inversion`);
- the `inversion-improvements` branch.

### v1 (the code behind the published Vs models)

| # | Finding | Effect |
|---|---|---|
| V1 | `Roughness(vs, 2)` is `sqrt(sum(np.gradient(vs, 2)))`. The `2` is the grid *spacing*, not a derivative order, and the sum of a first derivative telescopes. The term equals √((V_bottom − V_top)/2) exactly and does not change when structure inside the profile changes. | The "roughness regularization" applied no smoothing. It only penalized a faster bottom relative to the top. |
| V2 | `liklihood` uses `exp(-0.5 ‖r‖₂)`, the norm and not its square. | A Laplace-like, much flatter likelihood. Posterior widths are not Gaussian-error widths. |
| V3 | Each Vs move also sets Vp = √3 Vs (Poisson 0.25) and density from Gardner. | The text says Vp and density were held at CRUST1.0 values; the code did otherwise. |
| V4 | Out-of-bounds Vs proposals are redrawn in a loop (shallow layers) or reflected (deepest layer). | Asymmetric proposals near the bounds with no Hastings correction, so detailed balance fails there. |
| V5 | A thickness move adds the same step to the chosen layer and to the deepest layer, and can pick the half-space row. A negative thickness is redrawn with \|r\|. | The total thickness is not conserved, and the proposal is asymmetric. |
| V6 | Vs bounds are 0.15–1.2 × the starting model (0.9–1.1 for the deepest layer). The paper states 0 to +25 % (shallow) and −100 % to +25 % (deep). | The prior depends on the starting model and does not match the text. |
| V7 | The data uncertainty is the standard deviation of the per-window estimates (thesis §3.3.4: 1.25e-11 at RR52, 30–60 % of η). The thesis also describes annealing it during the run, but that code is commented out. | Errors were roughly √n too large for a stacked mean. The fit hardly constrained the model, so deep Vs came from the bounds. |
| V8 | Compliance computed before ComPy 2.0 multiplied the wave-attraction term by ω²/(ω²+c), i.e. `k·γ·ad·(pa_ratio + |Z/P|)` in `Calculate_Compliance`, so the term was effectively dropped. The same line writes the free-air constant as `3.07 * 10e-6`, which is 3.07e-5, ten times too large. | The published η is low: by up to 25–30 % at 7 mHz at RR52 from the dropped term, plus 1.4 % at 7 mHz from the constant. Both bias deep Vs high. 2.0 fixed both. |

Figure code (`plot_inversion_density_all`, `final_plot`, `start_model_plot*` in `inv_compy`):

| # | Finding | Effect |
|---|---|---|
| F1 | The heat maps average the 1-m Vs profile over 100 samples, so each of the 100 rows is 100 m. The y ticks put row 12 at −2000 m (167-m rows), and the CRUST1.0 reference lines are built on 167-m rows (interfaces at rows 14 and 42 = 2340 and 6900 m). | Posterior heat maps and medians in Figs. 9–10 are drawn 5/3 too deep relative to their labels, while the CRUST1.0 lines are at the right depths. A feature labelled at 3.5 km is at about 2.1 km. Check against the saved `Inversion_container`, which is not on this machine. |
| F2 | Serpentinization is computed as (9000 − 2 Vs)/30, i.e. 0 % at 4.5 km/s and 100 % at 3.0 km/s. The Fig. 10 caption states 0 % at 4.2 km/s and 100 % at 2.2 km/s (Miller & Christensen 1997). | Fig. 10c overstates serpentinization: Vs 3.5 km/s gives 67 % in the code and 35 % with the caption's relation. |

Documentation errata in the paper:
- Eq. 13 accepts unconditionally when L(m_i) ≤ L(m_{i−1}). With L a likelihood, this is reversed; the code was correct.
- The misfit (1/N)·√Σ(·)² is neither the χ² used in the code nor a standard one.
- The 70–90 % acceptance target is not Roberts et al. (1997), whose optimum is 23 % (44 % in 1-D).

### v2 and the `inversion-improvements` branch

These fixed V1–V6:
- a standard Gaussian log-likelihood;
- a second-derivative roughness term;
- symmetric proposals with rejection outside the bounds;
- conserved thickness;
- explicit bounds;
- an optional data covariance, hierarchical smoothing, tempering, Occam and kernels.

The forward model was checked against an independent Aki & Richards propagator.

Remaining issues:

- **A1.** The default `alpha=0.25` acts on Vs in m/s on a 1-m grid. A unit step has roughness 0.25, so the default is a Gaussian prior on Vs jumps with σ = 4 m/s. Unless `alpha_bounds` is used, the smoothing prior dominates any compliance data.
- **A2.** Vp and density are fixed. As section 4 shows, this makes the posterior multimodal for real data (R-hat up to 2.3 with plain chains) and shifts Vs. It is a modelling choice, not a neutral one.
- **A3.** No guidance on data errors. The 0.2 % Bendat & Piersol error of a stacked estimate is far smaller than the month-to-month changes (2–3 %, highly correlated across frequency) and the 5 % gain uncertainty.
- **A4.** Single-site Metropolis moves mix slowly, and every forward call costs about 0.3 ms.

The legacy functions are kept unchanged for reproducibility. The new workflow
is a separate module.

## 2. The 2.3 workflow (`compy_inversion`)

```python
import compy_inversion as ci

column = ci.crust1_column(lat, lon, "path/to/crust1")      # CRUST1.0 cell
data = ci.ComplianceData.from_npz("RR52_compliance_compy.npz",
                                  subsets=monthly_estimates,  # (n_months, n_freqs)
                                  gain_rel=0.05, floor_rel=0.01)
model = ci.CrustalModel(moho_depth=column["moho_depth"], sediment_max=1000.0)
post = ci.sample(data, model, n_chains=8, n_iter=80000, n_burn=30000, thin=10)
prior = ci.sample(data, model, prior_only=True, n_chains=8, n_iter=80000, n_burn=30000)
post.save("RR52_main")            # .npz draws + .json summary with R-hat/ESS
q = post.vs_quantiles(z)          # Vs-depth credible intervals
```

### Data and covariance

`ComplianceData.build` sums independent sources:

- **random**: the Bendat & Piersol error of the stacked estimate (diagonal).
- **subsets**: the covariance of independent sub-record estimates (e.g. monthly stacks) divided by their number, and shrunk halfway towards its diagonal. This captures non-stationarity and the correlation between frequency bins.
- **gain**: σ_g² v vᵀ, with σ_g = 5 % (Aminian et al. 2025, Table 1) and v = ∂η/∂ln(gain) = η − k·(a_w/p)/(ω² + 2g/R). The wave-attraction part does not scale with the pressure gain, and this is exact for `gravity_corrected_compliance`. A gain error shifts the whole curve together and is not counted once per frequency.
- **floor**: 1 % per frequency for what a 1-D model and the estimator neglect (2-D structure, Crawford et al. 1998; Hulme et al. 2005).

The likelihood is Gaussian with this covariance (Mahalanobis χ²).

### Parameterization (Zha & Webb 2016, §3.1)

| Part | Parameters | Prior |
|---|---|---|
| Sediment | thickness h | uniform on [0, 1000] m. Vs(z) and Vp(z) follow Hamilton (1976, 1979) terrigenous silt-clay profiles in 8 equal sublayers, with ρ = 2000 kg/m³. Optionally a Vs multiplier in [0.5, 2] (variant `sedscale`). |
| Crust | one Vs per layer: 100, 200, 400, 800 and 1600 m below the sediment, then one layer to the Moho | uniform on [500, 4700] m/s, plus a Gaussian prior on second differences of Vs whose scale σ ∈ [20, 3000] m/s is itself sampled (Jeffreys prior; Malinverno & Briggs 2004). |
| Mantle | Vs of the half-space | uniform on [3800, 4700] m/s |
| Moho | depth from CRUST1.0 at the station; optional ±1 km shift (variant `moho`, as in Aminian et al. 2025) | uniform |

Vp = 1.85 Vs and ρ = 1.85 + 0.165 Vp (g/cc, km/s; Christensen & Shaw 1970) in
the crust and mantle, as in Zha & Webb (2016). Alternatives:
- `vp_vs=None`: Brocher (2005) Vp(Vs).
- `reference=crust1_model(column)`: Vp and density fixed at CRUST1.0, as in Crawford et al. (1991) and Aminian et al. (2025).

Sediment Vp is the larger of Hamilton's Vp and Brocher's. Models with
Vp/Vs < √2 (negative Poisson ratio) are outside the prior.

Hamilton's Vs is used instead of the Ruan et al. (2014) profile of Zha & Webb.
The Ruan paper was not accessible, and its relation is not quoted in the
sources read. Any profile can be passed as `sediment_vs` / `sediment_vp`.

### Sampler

`dram_chain` is delayed-rejection adaptive Metropolis (Haario et al. 2001, 2006),
the DRAM used by Zha & Webb (2016):
- **Proposal.** A Gaussian whose covariance is that of the second half of the chain so far. A global scale is tuned by Robbins-Monro to 23 % acceptance (Andrieu & Thoms 2008; Roberts et al. 1997).
- **Adaptation.** It stops at the end of burn-in, so retained draws come from a fixed Metropolis kernel.
- **Delayed rejection.** A second stage, 10 times smaller, uses the Tierney & Mira (1999) acceptance.
- **Smoothing weight.** Gibbs-sampled after every step.
- **Tempering.** `temperatures` turns each chain into a parallel-tempering ensemble (Geyer 1991; Sambridge 2014): likelihood^(1/T), prior untempered, adjacent swaps, and only T = 1 kept.

Chains start from independent prior draws. `Posterior.diagnostics` gives the
rank-normalized split R-hat and bulk ESS (Vehtari et al. 2021) of every
parameter and of log α.

### Forward model and sensitivity

`Forward` is the minor-vector propagator of `inv_compy.raydep` (Crawford et al.
1991, from W. C. Crawford's code), evaluated at the seafloor only. With the
surface vector (0, −1/y₃, 0, 0), the downward pass gives v = 1 and
σ_zz = y₅/y₃, so η = −k y₃/(ω y₅). Numba compiles it when available: 2 µs
instead of 290 µs per call.

`vs_kernels` gives d ln η / d ln Vs per km, with Vp/Vs fixed.
`sensitivity_depth` gives the depth above which a fraction of the integrated
|kernel| lies. This is a physical measure of how deep the data see, which a
posterior/prior width ratio does not provide when the layer geometry is fixed
(see the limits below).

### Gain as an unknown and several epochs

- **Sampled gain.** `sample(..., gain_bounds=(0.95, 1.05))` adds a pressure-gain factor g, the true gain over the gain used, with a uniform prior: the ±5 % of Aminian et al. (2025, §4.1.5). The predicted data are η_a + g (η − η_a), with η_a the gain-independent wave attraction. The covariance must then exclude the gain term (`gain_rel=0`). Without it, the gain enters as the covariance term above; with a Gaussian prior the two are equivalent to first order.
- **Several epochs.** `ComplianceData.join(parts, names)` stacks estimates of one station, e.g. before and after a change in pressure sensitivity, each computed with its own Rayleigh gain. The covariance is block diagonal, and each part gets its own gain factor, while the Earth model is shared. The ratio of the factors tests whether the two calibrations agree.

### Estimator bias (H1 and H2)

Compliance uses H1 = |S_pz|/S_pp, unbiased when the noise is on the vertical. Noise on the pressure gauge biases H1 low by the squared coherence; H2 = S_zz/|S_pz| is unbiased in that case. The truth lies between the two (Bendat & Piersol 2010; Crawford & Singh 2008). With stacked coherence of 0.86–0.97, H2/H1 is 1.07–1.34 for the YV stations, largest at the lowest frequencies. Inverting H2 (variant `final_h2`) brackets this systematic error.

## 3. Validation (tests/test_compy_inversion.py)

- **Forward model:**
  - `Forward` equals `raydep` for 120 random models, including slow sediments with propagating shear waves, to 2e-6 relative. The remaining difference is round-off amplified by cancellation in the propagator (condition number up to 1e7 when c ≪ Vs); the median difference is 3e-9.
  - A uniform half-space tends to the static (λ+2μ)/(2μ(λ+μ)) as (c/Vs)².
- **Gain sensitivity:** the gain vector equals a finite-difference derivative of `gravity_corrected_compliance`.
- **Covariance:** the terms sum to the total; a pure gain shift costs χ² < 1.5, while the same amplitude alternating across frequency costs > 20.
- **Sampler:**
  - DRAM reproduces the mean and variance of a correlated, truncated 2-D Gaussian, with and without delayed rejection.
  - The hierarchical smoothing weight reproduces the analytic marginal.
  - The truncated Gamma draws match their analytic mean.
  - Parallel tempering samples both modes of a target with a log-density barrier of 200 equally, where a single adaptive chain stays in one mode.
- **Kernels:** they integrate to the response to a uniform change (3 %); lower frequencies sense deeper; sliver slices at the Moho are removed.
- **Synthetic recovery:** an RR52-like model with realistic covariance is recovered (truth inside every 90 % interval, R-hat < 1.1).
- **Joint data and gain:** a joint data set has a block-diagonal covariance and additive χ². The gain factor maps compliance to measured data exactly. A 3 % gain offset between two epochs is recovered to ±1 %, and the misfit trace includes burn-in.

## 4. Results for the YV broadband stations

See `docs/yv_inversion_2026-10-29/` (README, figures, tables, scripts).

## 5. Limits

- **The data see the top few kilometres.** Half of the lowest-frequency sensitivity lies above 0.2–2.1 km (station dependent) and 90 % above 5–8 km. Below that, Vs reflects the prior and the layer geometry.
- **Fixed layer geometries make deep posteriors look narrower than the data warrant.** A thick bottom layer averages Vs over several kilometres. Compare with the prior-only run and the kernels before interpreting deep layers.
- **Vs depends on the Vp assumption.** Vp/Vs = 1.85, Brocher, and Vp fixed at CRUST1.0 give different absolute Vs for the same data. With fixed Vp the posterior becomes multimodal and needs tempering.
- **Sediment and upper-crust Vs trade off.** A free sediment Vs scale (`sedscale`) changes the top kilometre at RR28 and RR34.
- **Low-frequency bias at RR38.** None of the parameterizations with Vp/Vs coupling fits RR38 (χ² 22 for 6 data): its low-frequency compliance is too low. Pressure-channel noise biases η low (Crawford & Singh 2008), and RR38 also has the largest month-to-month scatter. Treat RR38 as unreliable.
- The roughness-prior normalization ignores the truncation by the Vs bounds, which slightly favours smoother models when σ approaches the bound width.
- The 1-D model ignores the 2-D/3-D structure of ridge-axis sites (Crawford et al. 1998; Hulme et al. 2005), partly covered by the 1 % floor.

## References

- Aminian, M. A., et al. (2025). GJI 242, ggaf253.
- Andrieu, C., & Thoms, J. (2008). Stat. Comput. 18, 343–373.
- Brocher, T. M. (2005). BSSA 95, 2081–2092.
- Christensen, N. I., & Shaw, G. H. (1970). Geophys. J. R. Astron. Soc. 20, 271–284.
- Constable, S. C., Parker, R. L., & Constable, C. G. (1987). Geophysics 52, 289–300.
- Crawford, W. C., Webb, S. C., & Hildebrand, J. A. (1991, 1998). JGR 96, 16151; 103, 9895.
- Crawford, W. C. (2004). GJI 157, 1130–1145.
- Crawford, W. C., & Singh, S. C. (2008). Geophys. Prospect. 56, 313–325.
- Geyer, C. J. (1991). Computing Science and Statistics, 156–163.
- Haario, H., Saksman, E., & Tamminen, J. (2001). Bernoulli 7, 223–242.
- Haario, H., Laine, M., Mira, A., & Saksman, E. (2006). Stat. Comput. 16, 339–354.
- Hamilton, E. L. (1976). Geophysics 41, 985–996.
- Hamilton, E. L. (1979). JASA 65, 909–922; 66, 1093–1101.
- Hulme, T., Crawford, W. C., & Singh, S. C. (2005). GJI 163, 547–558.
- Laske, G., et al. (2013). CRUST1.0, EGU 15, 2658.
- Malinverno, A., & Briggs, V. A. (2004). Geophysics 69, 1005–1016.
- Mosher, S. G., et al. (2021). GJI 227, 1879–1892.
- Roberts, G. O., Gelman, A., & Gilks, W. R. (1997). Ann. Appl. Probab. 7, 110–120.
- Sambridge, M. (2014). GJI 196, 357–374.
- Tierney, L., & Mira, A. (1999). Stat. Med. 18, 2507–2515.
- Vehtari, A., et al. (2021). Bayesian Anal. 16, 667–718.
- Zha, Y., & Webb, S. C. (2016). JGR 121, 3220–3237.
