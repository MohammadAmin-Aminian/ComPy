# Data of the November 2026 pressure audit

Copied from the BRUIT-FM repository (`MohammadAmin-Aminian/bruit-fm`, commit
`8595452`), where the dated reports stay:

| File | Source in BRUIT-FM | Change |
|---|---|---|
| `pressure_qc_calibration_table.csv` | `docs/reports/2026-10-28_seismic_calibration/table/pressure_qc_calibration_table.csv` | Local file-path columns dropped (`inventory_file`, `ppsd_npz_file`, `*_result_file`, `*_plot_file`, `*_summary_file`, `tilt_vertical`); values unchanged |
| `metadata_offsets_channels.csv` | `docs/reports/2026-10-28_provisional_decisions/metadata_offsets_channels.csv` | None |

One row per pressure channel epoch (846), BRUIT-FM v1.2.0 with provisional
thresholds; no gain applied.
