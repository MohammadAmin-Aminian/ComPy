# YV microseism pressure–vertical check: all stations

Independent Claude implementation; branch `claude/microseism-yv`.
[Method and RR52 pilot](microseism_pressure_vertical_check.md).

> **Superseded for comparison purposes** by the all-network run in
> [microseism_all_networks.md](microseism_all_networks.md).
>
> * That run used each station's default centred week rather than the common
>   2013-06-10 week used here.
> * It read the pressure and vertical poles from each channel's response
>   snapshot, not from `yv_rayleigh_independent_validation.json`.
> * Some labels therefore differ, for example RR01 and RR26.
> * The claims marked *(revised)* below were corrected in the 15 October
>   self-review.
Status: `diagnostic_only`. No gain was fitted, applied or stored. No master,
score, saved product or raw file was changed; master hashes were checked
before and after the run.

## Run

```bash
bash scripts/run_yv_microseism_check.sh            # ~15 min; new run folder
PYTHONPATH=src python scripts/summarize_yv_microseism.py <run folder> --docs docs
```

Completed run: `results/EXPERIMENTS/MICROSEISM_YV/run_20261014T175553Z_30002f23`.

* **Coverage.** All 43 YV pressure epochs with completed PPSDs.
* **Interval.** The common week from 2013-06-10. Six stations whose coverage
  ends earlier (RR31, RR43–RR46, RR48) used the paired-spectra default week
  centred on their coverage.
* **Settings.** Depth from the master `Profondeur`, ρw = 1028 kg/m³,
  α = 1500 m/s, historical RR52 layer model; no reference gain going in.
* **Outcome.** 41 checked. RR12 (flat pressure) and RR42 (one day of pressure)
  are unresolved for the pairing reasons already documented.

Curated outputs:

* [yv_microseism_check.csv](yv_microseism_check.csv) — one row per station;
* [yv_microseism_check.json](yv_microseism_check.json) — phase-screen
  details;
* [figures/yv_microseism_check.png](figures/yv_microseism_check.png).

## Phase screen

Phase is the strongest network-wide result. For any layered elastic medium and
any mixture of plane waves, the P/Z transfer function is real. With sign −1 and
first-branch waves the phase is 0°, and it is never ±90°. The 0.1–0.3 Hz
microseism is coherent at almost every station, even where the long-wave band is
not.

**Method.** For bins with coh² ≥ 0.5 in 0.05–0.30 Hz (phase standard error ≈ 3°),
four hypotheses are compared:

1. nominal response;
2. the recorded positive-real Laplace pressure pole `+p` truly at `−p`
   (phase −2·atan(p/ω), amplitude unchanged);
3. hypothesis 1 with reversed relative polarity;
4. hypothesis 2 with reversed relative polarity.

The pole values come from
[yv_rayleigh_independent_validation.json](yv_rayleigh_independent_validation.json).
A label is assigned when the best circular RMS residual is ≤ 10°; above 20° the
station is "unexplained". These are screening labels, not scores.

| Label | Stations | Residual |
| --- | --- | --- |
| nominal | RR11, RR28, RR29, RR34, RR36, RR38, RR40, RR50, RR52 | 3–9° |
| reflected pressure pole | RR03, 05, 06, 07, 08, 09, 10, 13, 14, 16, 19, 20, 21, 25, 30, 41, 43, 44, 45, 46, 48, 55 | 1–10° |
| reversed polarity + reflected pole | RR47 | 3.0° |
| unexplained | RR17, RR18, RR22, RR26, RR56 | 23–60° |
| no coherent bins | RR01, RR31, RR53, RR57 | — |

* **Reflected pressure pole (22 of the 34 stations with a positive pole).**
  The observed phase follows −2·atan(p/ω) closely. At RR20, for example, it is
  −29°, −18°, −13°, −11° at 0.08, 0.12, 0.16, 0.25 Hz, against −26°, −17°,
  −13°, −8° predicted. A clock offset would instead grow linearly with
  frequency. This is independent, model-free support for treating those poles
  as sign errors in the metadata, as tested by the optional pole-reflection
  trial. At 0.03–0.07 Hz the same error is 20–40°, exceeding the Rayleigh
  calibration's 20° phase limit, which plausibly explains many of the original
  "phase convention" rejections.
* **RR11.** It has a recorded pole (0.0518 rad/s), but the nominal hypothesis
  fits better (3.9° against 7.3°). This needs a response check.
* **RR47.** The phase is about 169°, and 147.5° in the 17 coherent bins below
  its water resonance, where 152.6° is predicted. Below resonance the inference
  needs no structure model. Relative P/Z polarity is reversed compared with the
  other YV stations. Its vertical metadata is normal (dip −90°), so the cause
  is a pressure or vertical channel polarity that P–Z data alone cannot assign.
  RR47 has no Rayleigh gain; its sign assumption would fail that calibration's
  phase QC.
* **RR22 and RR56.** About +110° and −112°, which no wave mixture can produce
  and polarity plus the recorded pole cannot explain (best residual 60–68°).
  This points to a response or signal problem. RR22 was already unresolved in
  the Rayleigh runs.
* **RR17 and RR18.** These stations have much larger recorded poles (0.747 and
  0.766 rad/s, ≈ 0.12 Hz). Their phase follows about −atan(p/ω) (single), not
  −2·atan(p/ω).
* **RR26.** Its coherent low-frequency signal (0.029–0.037 Hz, coh² ≈ 0.75) has
  a phase of +75°. The coherent component there is not the seismic P/Z relation,
  so its amplitude is not interpreted.
* **No coherent bins.** RR01, RR31, RR53 and RR57 have no coherent P/Z bins even
  in the secondary microseism, which suggests a channel problem.

## Gain consistency

A station counts as gain-informative when two conditions hold:

* its best fixed sub-band has conditioning ratio (structure-free H2 bound ÷
  model-conditional H1 value) ≤ 1.5;
* its phase is explained.

Exactly the eight Trillium-240 stations analysed by
[Aminian et al. (2025)](https://doi.org/10.1093/gji/ggaf253) qualify. Elsewhere
the 0.03–0.1 Hz band is incoherent in ambient noise (Ruan et al. 2014, Fig. 3)
and the amplitude is not interpreted.

| Station | Depth (m) | Sub-band (Hz) | Microseism g (H1–H2) | Baseline g | vs baseline | Aminian 2025 | vs published |
| --- | --- | --- | --- | --- | --- | --- | --- |
| RR28 | 4550 | 0.045–0.058 | 0.734–0.811 | 0.820 | −10.5…−1.1% | 0.78 | −5.9…+3.9% |
| RR29 | 4829 | 0.034–0.043 | 0.684–0.730 | 0.730 | −6.3…+0.1% | 0.70 | −2.3…+4.3% |
| RR34 | 4265 | 0.048–0.061 | 0.564–0.582 | 0.585 | −3.6…−0.5% | 0.58 | −2.7…+0.4% |
| RR36 | 3560 | 0.055–0.072 | 0.623–0.651 | 0.654 | −4.7…−0.5% | 0.64 | −2.7…+1.7% |
| RR38 | 4560 | 0.045–0.058 | 0.612–0.615 | 0.636 | −3.8…−3.3% | 0.60 | +2.0…+2.6% |
| RR40 | 4780 | 0.034–0.044 | 0.713–0.751 | 0.743 | −3.9…+1.1% | 0.74 | −3.6…+1.4% |
| RR50 | 4118 | 0.049–0.063 | 0.762–0.859 | 0.717 | +6.4…+19.9% | 0.68 | +12.1…+26.3% |
| RR52 | 2918 | 0.047–0.063 | 0.671–0.681 | 0.679 | −1.2…+0.4% | 0.67 | +0.1…+1.7% |

* **Agreement.** Six of the eight brackets lie within about 4% of both the
  baseline and the published gains. The H1–H2 bracket contains the published
  value at RR28, RR29, RR34, RR36 and RR40, and the baseline value at RR29,
  RR40 and RR52.
* **RR50.** Microseism is higher than both references. RR50 is at the Rodrigues
  Triple Junction, where Aminian et al. (2025, Fig. 9g) find upper-crustal Vs of
  only 2–2.4 km/s. The historical RR52 model is faster, which would overestimate
  the admittance height and therefore g. *(revised)* This is a hypothesis, not a
  demonstrated model limitation: the Rayleigh reference uses the same model (in
  a different band), and no local-model run has been made.
* **RR38.** The structure-free bound flags the baseline 0.636 in 29% of
  long-wave bins (RR36: 11%). The published 0.60 satisfies it. The bound
  scales with depth (±2% depth moves it about ±3.5%), so this is a lead, not a
  conclusion.
* **Caveat.** The values are model-conditional and share the reference's water
  and Earth-model assumptions. They test agreement between noise and
  earthquakes, not absolute calibration. *(revised)* The Aminian et al. (2025)
  gains come from the same teleseismic Rayleigh P/a method family as the
  baseline, so agreement with both is not two independent confirmations.

## Leads for review

1. Cross-check the reflected-pole phase result against the pole-reflection
   trial and the RESIF/instrument documentation, without silently changing
   metadata.
2. RR47: determine whether the pressure or the vertical channel is reversed,
   for example from the first-motion polarity of teleseismic P waves on the
   vertical against neighbouring stations.
3. Review the responses of RR22, RR56, RR17, RR18 and RR11, and the RR26
   low-frequency signal.
4. Check the channels of RR01, RR31, RR53 and RR57, which show no coherent P/Z
   bins.
5. RR50/RR36/RR38: repeat the model-conditional values with local crustal models
   (for example the Aminian et al. 2025 Vs profiles) and verified depths.
6. For a structure-free per-band gain, use near-vertical teleseismic body waves
   (Deng et al. 2022). This also tests for a frequency-dependent DPG response.

## Limitations

* One week per station. Hourly windows are correlated, so there are no
  confidence intervals; H1–H2 is a noise bracket.
* No tilt correction and no horizontal channels.
* Phase screening uses coherent bins only. Above the first water resonance it
  assumes first-branch (fundamental-mode) dominance; below it, the inference
  needs no structure model.
* All results constrain g/g_z and relative P/Z polarity. The vertical response
  is assumed correct.
