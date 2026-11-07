"""Pressure-gauge audit (November 2026): anomaly scan, figures and numbers.

ObsPy-free. Runs on the BRUIT-FM master table copied into
``pressure/docs/audit_2026-11/data`` and writes the figures, the ranked
anomaly table and ``numbers.json`` quoted in that folder's REPORT.md::

    MPLBACKEND=Agg python pressure/scripts/pressure_audit_2026_11.py
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from compy import anomaly as an  # noqa: E402
from compy.plotting.diagnostics import plot_anomaly_ranking  # noqa: E402
from compy_pressure import audit_plots, table_audit, water_layer  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "pressure/docs/audit_2026-11"
DISPERSION = (
    ROOT / "compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv"
)


def save(fig, path):
    fig.savefig(path, dpi=110)
    plt.close(fig)
    print("wrote", path.relative_to(ROOT))


def water_layer_table(dispersion):
    """Bias of a 30-70 mHz Rayleigh multiplier for a perfect gauge, by depth."""
    f = np.arange(0.03, 0.0701, 0.0025)
    c = np.interp(f, dispersion[:, 0], dispersion[:, 1])
    rows = []
    for h in (1000, 2000, 3000, 4000, 4500, 5000, 5400):
        keep = water_layer.usable(f, c, h)
        sens = water_layer.log_sensitivity(f[keep], c[keep], h)
        t = water_layer.ratio(f[keep], c[keep], h)
        weight = t**2 / np.sum(t**2)
        rows.append(dict(
            depth_m=h, bins_kept=int(keep.sum()), bins=len(f), max_x=round(float(sens["x"].max()), 2),
            top3_weight_pct=round(100 * float(np.sort(weight)[-3:].sum()), 0),
            max_sens_sound_speed=round(float(np.abs(sens["sound_speed"]).max()), 2),
            sound_speed_plus1_pct=round(100 * (water_layer.gain_bias(f, c, h, true_sound_speed=1515.0) - 1), 2),
            sound_speed_plus1_pct_log=round(100 * (water_layer.gain_bias(f, c, h, true_sound_speed=1515.0,
                                                                          weights="log") - 1), 2),
            phase_velocity_minus5_pct=round(100 * (water_layer.gain_bias(f, c, h, true_phase_velocity=0.95 * c) - 1), 2),
            phase_velocity_plus5_pct=round(100 * (water_layer.gain_bias(f, c, h, true_phase_velocity=1.05 * c) - 1), 2),
            column_mean_density=round(float(water_layer.column_mean_density(h)), 1),
            density_pct=round(100 * (water_layer.gain_bias(
                f, c, h, true_density=float(water_layer.column_mean_density(h))) - 1), 2),
        ))  # fmt: skip
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--table", type=Path, default=AUDIT / "data/pressure_qc_calibration_table.csv"
    )
    parser.add_argument(
        "--offsets", type=Path, default=AUDIT / "data/metadata_offsets_channels.csv"
    )
    parser.add_argument("--output", type=Path, default=AUDIT)
    args = parser.parse_args()
    figures = args.output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    table = table_audit.load_table(args.table, args.offsets)
    found, info = table_audit.scan_table(table)
    ranked = an.rank(found, 2.0)
    (args.output / "anomalies.md").write_text(an.to_markdown(ranked) + "\n")
    an.to_csv(ranked, args.output / "anomalies.csv")
    dispersion = np.loadtxt(DISPERSION, delimiter=",", skiprows=1)
    save(audit_plots.plot_gain_consistency(table), figures / "gain_consistency.png")
    save(audit_plots.plot_gain_depth(table), figures / "gain_vs_depth.png")
    save(audit_plots.plot_water_layer_sensitivity(dispersion[:, 0], dispersion[:, 1]),
         figures / "water_layer_sensitivity.png")  # fmt: skip
    save(
        audit_plots.plot_infragravity_residuals(table),
        figures / "infragravity_residuals.png",
    )
    save(audit_plots.plot_network_scores(table), figures / "network_scores.png")
    save(plot_anomaly_ranking(ranked, top=30, title="Ranked pressure-channel anomalies"),
         figures / "anomaly_ranking.png")  # fmt: skip
    g = table_audit.gains(table)
    agreement = {}
    for a, b in (("rayleigh", "microseism"), ("rayleigh", "local_p"), ("microseism", "local_p"),
                 ("rayleigh", "peer_offset"), ("microseism", "peer_offset")):  # fmt: skip
        both = np.isfinite(g[a]) & np.isfinite(g[b])
        d = 20 * np.log10(g[a][both] / g[b][both])
        agreement[f"{a}/{b}"] = dict(n=int(both.sum()), median_db=round(float(np.median(d)), 2),
                                     within_3db=int(np.sum(np.abs(d) <= 3)), within_6db=int(np.sum(np.abs(d) <= 6)))  # fmt: skip
    by_type = {}
    for method in ("rayleigh", "microseism", "local_p"):
        for kind in audit_plots.TYPES:
            sel = np.isfinite(g[method]) & (table["type"] == kind)
            if sel.sum() >= 3:
                by_type[f"{method}/{kind}"] = dict(n=int(sel.sum()), median_g=round(float(np.median(g[method][sel])), 3),
                                                   p10=round(float(np.percentile(g[method][sel], 10)), 3),
                                                   p90=round(float(np.percentile(g[method][sel], 90)), 3))  # fmt: skip
    res = table["infragravity_median_residual_db"]
    ig = {}
    for cls in ("good", "all"):
        for basin in ("Pacific", "Atlantic", "Indian", "Mediterranean"):
            for kind in (*audit_plots.TYPES, "all"):
                sel = np.isfinite(res) & (table["basin"] == basin)
                sel &= (table["final_class"] == cls) if cls != "all" else True
                sel &= (table["type"] == kind) if kind != "all" else True
                if sel.sum() >= 5:
                    ig[f"{cls}/{basin}/{kind}"] = dict(n=int(sel.sum()), median_db=round(float(np.median(res[sel])), 1),
                                                       q25=round(float(np.percentile(res[sel], 25)), 1),
                                                       q75=round(float(np.percentile(res[sel], 75)), 1))  # fmt: skip
    numbers = dict(n_channels=info["n_channels"], n_with_gain=info["n_with_gain"], n_anomalies=len(ranked),
                   kinds={k: sum(a.kind == k for a in ranked) for k in sorted({a.kind for a in ranked})},
                   agreement=agreement, gains_by_type=by_type, depth=info["depth"], rolloff=info["rolloff"],
                   polarity={k: {n: list(v) for n, v in d.items()} for k, d in info["polarity"].items()},
                   infragravity=ig, water_layer=water_layer_table(dispersion))  # fmt: skip
    (args.output / "numbers.json").write_text(
        json.dumps(numbers, indent=1, default=str) + "\n"
    )
    print("wrote", (args.output / "numbers.json").relative_to(ROOT))


if __name__ == "__main__":
    main()
