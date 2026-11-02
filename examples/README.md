# Examples

Run every example from the repository root. The offline examples need no data;
CI runs `synthetic_inversion.py`, `reproducible_inversion.py` and
`5_Pressure_QC.py` on every push.

## Offline (synthetic data)

| Script | What it shows |
|---|---|
| [`synthetic_inversion.py`](synthetic_inversion.py) | Forward model and a short inversion on a synthetic compliance curve |
| [`reproducible_inversion.py`](reproducible_inversion.py) | Four-chain benchmark with saved arrays, checksums and diagnostics |
| [`converged_inversion.py`](converged_inversion.py) | Converged four-chain run with adaptive steps; plot with `plot_posterior_heatmap.py` |
| [`deep_crust_inversion.py`](deep_crust_inversion.py) | Deep-crust test (4200 m water, 3.5–18.5 mHz; about 25 min on two cores); plot with `plot_deep_posterior.py` |
| [`4_Vs_inversion.py`](4_Vs_inversion.py) | Bayesian Vs inversion with `compy.inversion.dram` (synthetic without arguments) |
| [`5_Pressure_QC.py`](5_Pressure_QC.py) | Pressure-gauge QC on a generated SDS archive (synthetic without arguments, a few seconds) |

```bash
MPLBACKEND=Agg python examples/synthetic_inversion.py
MPLBACKEND=Agg python examples/reproducible_inversion.py --output synthetic-results
python examples/5_Pressure_QC.py --output qc-example
```

## Real data, in workflow order

| Step | Script | Input |
|---|---|---|
| 0 | [`0_Download_data.py`](0_Download_data.py) | FDSN download, decimation and response removal for one station and interval |
| 1 | [`1_Transients_EQ.py`](1_Transients_EQ.py) | Earthquake masking and periodic-transient removal on the vertical |
| 2 | [`2_Rotation.py`](2_Rotation.py) | Tilt correction of response-corrected displacement |
| 3 | [`3_DPG_Calibration.py`](3_DPG_Calibration.py) | Rayleigh-wave DPG gain from an SDS archive, StationXML and a ComCat catalogue |
| 4 | [`4_Vs_inversion.py`](4_Vs_inversion.py) | Vs inversion of a `stacked_compliance` result (`.npz`) |
| 5 | [`5_Pressure_QC.py`](5_Pressure_QC.py) | QC of every pressure gauge in an SDS archive |

```bash
python examples/0_Download_data.py YV RR52 2012-12-01 2012-12-02 output.mseed
python examples/3_DPG_Calibration.py SDS/ stationxml/ YV RR52 events.geojson
python examples/4_Vs_inversion.py RR52_compliance_compy.npz \
    --lat -20.4723 --lon 68.1094 --crust1 path/to/crust1.0
```

For a whole network, use the [pipeline](../pipeline/README.md) instead.
