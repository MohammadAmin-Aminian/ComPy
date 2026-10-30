# Conditional microseism estimates: frequency and week consistency

Rule `claude-frequency-v1-20261015` was designed on YL; all values assume the vertical
response and historical RR52 model. No calibration is accepted/applied.
At least two eligible bands <=0.10 Hz must agree within 10%, and at least
three independent passing weeks must agree within 10% (stable <=5%).
Legacy best-band values remain in microseism_weeks_table.md.
H1-H2 envelopes and temporal spread are descriptive, not confidence limits.

| Station | Passing / checked weeks | Status | Conditional g | Week spread | Rayleigh baseline | Rayleigh trial | Published scalar g |
|---|---:|---|---:|---:|---:|---:|---:|
| 7D.J28A | 0/0 | no_passing_weeks | — | — | — | — | — |
| 7D.J29A | 2/5 | insufficient_passing_weeks | — | 0.8 | — | — | — |
| 7D.J30A | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.981 | 0.981 | — |
| 7D.J35A | 0/5 | no_passing_weeks | — | — | — | — | — |
| 7D.J36A | 0/5 | no_passing_weeks | — | — | — | — | — |
| 7D.J38A | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.780 | 0.780 | — |
| 7D.J43A | 0/5 | no_passing_weeks | — | — | — | — | — |
| 7D.J44A | 0/5 | no_passing_weeks | — | — | — | — | — |
| 7D.J45A | 3/4 | stable | 0.794 | 4.8 | 0.829 | 0.829 | — |
| 7D.J52A | 0/4 | no_passing_weeks | — | — | — | — | — |
| 7D.J53A | 0/0 | no_passing_weeks | — | — | — | — | — |
| 7D.J54A | 3/4 | stable | 1.020 | 2.7 | 1.080 | 1.080 | — |
| 7D.J55A | 1/5 | insufficient_passing_weeks | — | 0.0 | — | — | — |
| 7D.J63A | 0/0 | no_passing_weeks | — | — | — | — | — |
| 7D.J67A | 3/4 | moderate | 0.973 | 5.1 | — | — | — |
| XO.LA33 | 0/5 | no_passing_weeks | — | — | 0.934 | 0.934 | — |
| XO.LA34 | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.963 | 0.963 | — |
| XO.LT02 | 0/5 | no_passing_weeks | — | — | — | — | — |
| XO.WD48 | 0/1 | no_passing_weeks | — | — | — | — | — |
| XO.WD66 | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.909 | 0.909 | — |
| XO.WD70 | 2/5 | insufficient_passing_weeks | — | 1.7 | 0.930 | 0.930 | — |
| YL.A02W | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.696 | 0.696 | 0.667 |
| YL.A06W | 1/4 | insufficient_passing_weeks | — | 0.0 | 0.692 | 0.692 | 0.654 |
| YL.A08W | 0/5 | no_passing_weeks | — | — | 0.730 | 0.730 | 0.699 |
| YL.A10W | 0/5 | no_passing_weeks | — | — | 0.585 | 0.585 | 0.549 |
| YL.A12W | 2/5 | insufficient_passing_weeks | — | 1.6 | 0.639 | 0.639 | 0.599 |
| YL.A14W | 0/5 | no_passing_weeks | — | — | — | — | 0.508 |
| YL.B01W | 0/5 | no_passing_weeks | — | — | — | — | 5.263 |
| YL.B04W | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.739 | 0.739 | 0.709 |
| YL.B09W | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.793 | 0.793 | 0.758 |
| YL.B11W | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.559 | 0.559 | 0.526 |
| YL.B14W | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.650 | 0.650 | 0.617 |
| YL.C01W | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.576 | 0.576 | 0.541 |
| YL.C03W | 0/5 | no_passing_weeks | — | — | 0.688 | 0.688 | 0.645 |
| YL.C05W | 3/5 | stable | 0.575 | 4.4 | 0.585 | 0.585 | 0.565 |
| YL.C07W | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.269 | 0.526 | 0.503 |
| YL.C08W | 0/5 | no_passing_weeks | — | — | 0.593 | 0.593 | 0.562 |
| YL.C09W | 3/5 | stable | 0.584 | 2.0 | 0.575 | 0.575 | 0.543 |
| YL.C10W | 0/5 | no_passing_weeks | — | — | — | — | 0.662 |
| YL.C11W | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.578 | 0.578 | 0.546 |
| YL.C16W | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.591 | 0.591 | 0.559 |
| YL.C17W | 2/5 | insufficient_passing_weeks | — | 0.4 | 0.731 | 0.731 | 0.699 |
| YL.F02W | 2/5 | insufficient_passing_weeks | — | 1.2 | 0.732 | 0.732 | 0.680 |
| YL.N02W | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.610 | 0.610 | 0.578 |
| YL.N03W | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.581 | 0.581 | 0.556 |
| YL.S01W | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.594 | 0.594 | 0.552 |
| YV.RR01 | 0/5 | no_passing_weeks | — | — | — | 1.041 | — |
| YV.RR03 | 0/5 | no_passing_weeks | — | — | 1.124 | 1.142 | — |
| YV.RR05 | 0/5 | no_passing_weeks | — | — | — | 1.120 | — |
| YV.RR06 | 0/5 | no_passing_weeks | — | — | — | 1.140 | — |
| YV.RR07 | 0/5 | no_passing_weeks | — | — | 1.079 | 1.067 | — |
| YV.RR08 | 1/5 | insufficient_passing_weeks | — | 0.0 | 1.160 | 1.160 | — |
| YV.RR09 | 1/5 | insufficient_passing_weeks | — | 0.0 | 1.194 | 1.175 | — |
| YV.RR10 | 1/5 | insufficient_passing_weeks | — | 0.0 | — | 1.084 | — |
| YV.RR11 | 1/5 | insufficient_passing_weeks | — | 0.0 | 1.083 | 1.054 | — |
| YV.RR12 | 0/0 | no_passing_weeks | — | — | — | — | — |
| YV.RR13 | 0/5 | no_passing_weeks | — | — | 1.224 | 1.201 | — |
| YV.RR14 | 1/5 | insufficient_passing_weeks | — | 0.0 | 1.096 | 1.065 | — |
| YV.RR16 | 0/5 | no_passing_weeks | — | — | — | 1.078 | — |
| YV.RR17 | 0/5 | no_passing_weeks | — | — | — | — | — |
| YV.RR18 | 0/5 | no_passing_weeks | — | — | — | — | — |
| YV.RR19 | 0/5 | no_passing_weeks | — | — | — | 1.130 | — |
| YV.RR20 | 1/5 | insufficient_passing_weeks | — | 0.0 | — | 1.342 | — |
| YV.RR21 | 0/5 | no_passing_weeks | — | — | — | 1.114 | — |
| YV.RR22 | 0/5 | no_passing_weeks | — | — | — | — | — |
| YV.RR25 | 0/5 | no_passing_weeks | — | — | 1.163 | 1.157 | — |
| YV.RR26 | 0/5 | no_passing_weeks | — | — | — | 1.089 | — |
| YV.RR28 | 2/5 | insufficient_passing_weeks | — | 0.7 | 0.820 | 0.820 | 0.780 |
| YV.RR29 | 2/5 | insufficient_passing_weeks | — | 2.5 | 0.730 | 0.730 | 0.700 |
| YV.RR30 | 0/5 | no_passing_weeks | — | — | — | 1.186 | — |
| YV.RR31 | 0/3 | no_passing_weeks | — | — | — | — | — |
| YV.RR34 | 1/3 | insufficient_passing_weeks | — | 0.0 | 0.585 | 0.585 | 0.580 |
| YV.RR36 | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.654 | 0.654 | 0.640 |
| YV.RR38 | 2/5 | insufficient_passing_weeks | — | 0.3 | 0.636 | 0.636 | 0.600 |
| YV.RR40 | 4/5 | stable | 0.749 | 3.2 | 0.743 | 0.743 | 0.740 |
| YV.RR41 | 0/3 | no_passing_weeks | — | — | — | — | — |
| YV.RR42 | 0/0 | no_passing_weeks | — | — | — | — | — |
| YV.RR43 | 0/3 | no_passing_weeks | — | — | — | 1.140 | — |
| YV.RR44 | 0/3 | no_passing_weeks | — | — | — | 1.307 | — |
| YV.RR45 | 0/1 | no_passing_weeks | — | — | — | — | — |
| YV.RR46 | 0/3 | no_passing_weeks | — | — | — | 1.015 | — |
| YV.RR47 | 0/3 | no_passing_weeks | — | — | — | — | — |
| YV.RR48 | 0/3 | no_passing_weeks | — | — | — | 1.387 | — |
| YV.RR50 | 1/5 | insufficient_passing_weeks | — | 0.0 | 0.717 | 0.717 | 0.680 |
| YV.RR52 | 3/5 | stable | 0.685 | 1.1 | 0.679 | 0.679 | 0.670 |
| YV.RR53 | 0/5 | no_passing_weeks | — | — | — | — | — |
| YV.RR55 | 0/5 | no_passing_weeks | — | — | — | 1.592 | — |
| YV.RR56 | 0/5 | no_passing_weeks | — | — | — | — | — |
| YV.RR57 | 0/5 | no_passing_weeks | — | — | — | — | — |
