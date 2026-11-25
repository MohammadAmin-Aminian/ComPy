# Upper end of the compliance band (ComPy 1.3.2)

**Decision.** Stage C searches the compliance band up to **0.8 f_c**, where f_c = √(g / 2πH) is the frequency at
which the water wavelength equals the depth (kH = 2π). In deep water 0.8 f_c is **kH = 4.0** (wavelength 1.57 H).
The inversion has dropped bins above 0.8 f_c since ComPy 1.3.0 (`auto.compliance_data(max_fraction_fc=0.8)`); stage C
now uses the same limit, so the windows are selected on the band that is inverted.
Setting: `BAND_MAX_FRACTION_FC` in `pipeline/c_compliance.py` (1.0 restores ComPy 1.2.2-1.3.1).

## What the papers do

| Paper | Water depth | Upper limit used or stated | In units of f_c |
|---|---|---|---|
| Crawford, Webb & Hildebrand (1991), JGR 96 | — | waves of wavelength H/2 to 2H reach the seafloor: f_max = √(g / 2π n H), 0.5 < n < 2, depending on the surface amplitude | 0.71-1.41 |
| Crawford et al. (1998), JGR 103 | about 2.2 km | 0.0035-0.03 Hz (wavelengths 2-40 km) | up to 1.13 |
| Willoughby & Edwards (2000), GRL 27 | 1.3 km | the rule of Crawford et al. (1991) with n = 0.5: up to 0.049 Hz | 1.41 |
| Crawford & Singh (2008) | — | "approximately the frequency at which the ocean surface gravity wavelength equals the water depth" | 1.0 |
| Ruan et al. (2014), JGR | 3 km | "the high-frequency cutoff is 0.02 to 0.03 Hz" | 0.88-1.32 |
| Bell, Forsyth & Ruan (2015), BSSA 105 | — (regression over their stations) | f_c (their eq. 11) has the slope of the measured coherence limit; the measured limit is about 10 % higher | about 1.1 |
| Zha & Webb (2016), JGR 121 | about 2 km | 0.006-0.02 Hz | 0.72 |
| Doran & Laske (2021), GJI 224 | 4.4 km and more | usable at periods above about 60 s; curves span 70-140 s | 0.76-0.88 |
| Mosher et al. (2021), G-cubed 22 | per station | chosen from the coherence and the sensitivity kernels of each station | per station |
| Aminian et al. (2025), GJI 242 | — | widest run with P-Z coherence above the threshold, below f_c | below 1.0 |

Depths are given where the paper states them for the quoted band. The papers agree that the limit is set by the decay of the surface-wave pressure with depth, as 1/cosh(kH), and
place it between 0.7 and 1.4 f_c. The upper values come from shallow sites or from the coherence alone.

## What the archive shows

- **Pressure.** Of the 379 screened gauges of the pressure atlas, the 164 that have an infragravity plateau and fall
  20 dB below it do so at a fixed kH = 4.27, interquartile 4.07-4.33, from 500 to 5500 m and for every sensor type
  (124 gauges deeper than 500 m); the period of that edge grows as
  H^0.44 (network-bootstrap interval 0.39-0.52; theory 0.5). At f_c the pressure is 33 dB below the plateau
  (interquartile 30-36 dB) at gauges deeper than 1 km. kH = 4.27 is 0.82 f_c.
- **Coherence.** In stage C of the 287 stations (wide gates, the input of stage E), the median Z-P coherence of the
  bins is 0.95 at 0.5-0.8 f_c, **0.76** at 0.8-0.9 f_c and **0.45** at 0.9-1.0 f_c (173 stations).
- **Shape of the curve.** Extrapolating each station's curve from its lower bins: one step lower in the band
  (0.55-0.7 f_c from 0.35-0.55 f_c) the compliance keeps rising, +4.6 % (interquartile +2.3 to +8.0 %); at 0.8-0.9 f_c
  the growth stops (+0.7 %, interquartile -5.0 to +3.4 %); at 0.9-1.0 f_c it collapses for a quarter of the stations
  (-2.7 %, interquartile -27 to +6 %). Where the pressure signal fades, noise on the pressure biases S_pp up and
  the compliance down.
- **Band ends of ComPy 1.3.1.** Median 0.75 f_c (kH 3.5); 77 of 287 stations ended above 0.8 f_c, none above f_c.

## Effect of the change

Restacked from the cached spectra for the 77 stations whose band ended above 0.8 f_c (wide gates):

| Quantity | Median | Interquartile | Largest |
|---|---|---|---|
| Band end, before (f_c) | 0.87 | 0.84-0.91 | 0.99 |
| Band end, after (f_c) | 0.78 | 0.76-0.79 | 0.80 |
| Windows kept | +0.8 % | +0.4 to +1.0 % | +2.2 % |
| Compliance in the bins both keep (median per station) | -0.01 % | -0.03 to +0.01 % | 0.69 % |
| Relative error of those bins | 0.13 % | 0.10-0.19 % | |

No station changed its coherence threshold. The other 210 stations do not change. Vs does not change, since the
inversion already dropped the bins above 0.8 f_c.

Scripts and tables: run folder `runs/pressure_2026-11-18` of the ComPy workspace (`s02_band_limit.py`,
`s03_band_restack.py`, `s05_band_effect.py`; `tables/band_*`). Test: `tests/test_paper_check.py::test_stage_c_band_search_stops_at_kh_4`.
