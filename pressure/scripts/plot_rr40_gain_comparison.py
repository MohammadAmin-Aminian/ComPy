"""Illustrative gain shifts of saved YV PSDs; never modify source products."""

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy.imaging.cm import pqlx
from obspy.signal import PPSD

from compy_pressure.plot_network_ppsd import add_references, sha256
from compy_pressure.pressure_test import load_references

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10


def pressure_gain_db(gain):
    """An amplitude multiplier changes PSD by 10 log10(g²), not 10 log10(g)."""
    gain = float(gain)
    if not np.isfinite(gain) or gain <= 0:
        raise ValueError("Pressure gain must be finite and positive")
    return float(20 * np.log10(gain))


def band_centre_shifts(frequencies, bands):
    """Illustration only; cannot unsmooth saved one-octave PSDs at band edges."""
    shift = np.full(len(frequencies), np.nan)
    occupied = np.zeros(len(frequencies), dtype=bool)
    for band, item in bands.items():
        lo, hi = map(float, band.split("-"))
        if not 0 < lo < hi:
            raise ValueError("Invalid local-P frequency band")
        keep = (frequencies >= lo) & (frequencies < hi)
        if np.any(occupied & keep):
            raise ValueError("Overlapping local-P display bands")
        shift[keep] = pressure_gain_db(item["gain_median"])
        occupied |= keep
    return shift


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--station", choices=["RR40", "RR52"], default="RR40")
    parser.add_argument("--show", action="store_true", help="Open the saved comparison")
    args = parser.parse_args()
    station = args.station
    root = Path(__file__).resolve().parents[1]
    review_path = root / "results/REVIEWS/SIX_CHECKS/final_20261019/review.json"
    review = json.loads(review_path.read_text())
    rows = [r for r in review["rows"] if r["station"] == f"YV.{station}"]
    if len(rows) != 1:
        raise ValueError(f"Expected one exact {station} source epoch")
    row = rows[0]
    epoch = {
        "RR40": "20121016T060000",
        "RR52": "20121019T090000",
    }[station]
    source = (
        root
        / f"results/NETWORKS/TEST/YV/{station}"
        / f"YV.{station}.00.BDH__{epoch}__58afc2d6f005/ppsd.npz"
    )
    protected = {str(p): sha256(p) for p in (source, review_path)}
    reference_dir = root / "../compy_pressure/reference_data"
    reference_hashes = {
        str(path): sha256(path)
        for path in (
            reference_dir / "pressure_references.json",
            reference_dir / "brown_hydroacoustic_fullband.csv",
        )
    }
    p = PPSD.load_npz(str(source), allow_pickle=False)
    if p.id != row["Seed_id"] or p.special_handling != "hydrophone":
        raise ValueError("Pressure source identity mismatch")
    order = np.argsort(1 / p.period_bin_centers)
    f = (1 / p.period_bin_centers)[order]
    values = np.asarray(p.psd_values, dtype=float)[:, order]
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite PSD")
    edges = np.r_[1 / p.period_xedges[::-1]]
    checks = row["checks"]
    ray = float(checks["rayleigh"]["metrics"]["pressure_multiplier"])
    micro = float(checks["microseism"]["metrics"]["microseism_conditional_g"])
    bands = json.loads(checks["local_p"]["metrics"]["local_p_gain_by_band"])
    local_shift = band_centre_shifts(f, bands)
    variants = [
        ("01_nominal", "Nominal response; no additional gain", np.zeros(len(f))),
        (
            "02_rayleigh",
            f"Rayleigh g={ray:.4f} ({pressure_gain_db(ray):+.2f} dB); "
            "hypothetical full band",
            np.full(len(f), pressure_gain_db(ray)),
        ),
        (
            "03_microseism",
            f"Microseism g={micro:.4f} ({pressure_gain_db(micro):+.2f} dB); "
            "hypothetical full band",
            np.full(len(f), pressure_gain_db(micro)),
        ),
        (
            "04_local_p",
            "Local-P: approximate band-centre shifts; estimated bands only",
            local_shift,
        ),
    ]
    out = (
        root
        / f"results/REVIEWS/{station}_GAIN_PPSD"
        / datetime.now(UTC).strftime("run_%Y%m%dT%H%M%S%fZ")
    )
    out.mkdir(parents=True, exist_ok=False)
    db = np.arange(np.floor(values.min()) - 15, np.ceil(values.max()) + 5, 1)
    fig, axes = plt.subplots(
        2, 2, figsize=(16, 10), sharex=True, sharey=True, layout="constrained"
    )
    for (name, title, shift), ax in zip(variants, axes.flat):
        shifted = values + shift
        probability = np.full((len(db) - 1, len(f)), np.nan)
        for j in np.flatnonzero(np.isfinite(shift)):
            counts = np.histogram(shifted[:, j], db)[0]
            assert counts.sum() == len(values)
            probability[:, j] = counts * 100 / len(values)
        probability = np.ma.masked_invalid(np.ma.masked_equal(probability, 0))

        def draw(axis):
            mesh = axis.pcolormesh(
                edges, db, probability, cmap=pqlx, vmin=0, vmax=30, shading="flat"
            )
            add_references(axis, True, f.min(), f.max())
            # Approximate environmental regions from Webb, not sensor passbands.
            for lo, hi, label in (
                (f.min(), 0.03, "IG / depth-dependent"),
                (0.03, 0.1, "notch*"),
                (0.1, min(5, f.max()), "pressure microseisms*"),
            ):
                if hi > lo:
                    axis.text(
                        np.sqrt(lo * hi),
                        0.98,
                        label,
                        transform=axis.get_xaxis_transform(),
                        ha="center",
                        va="top",
                        fontsize=8,
                        color="0.3",
                    )
            axis.set(
                xscale="log",
                xlim=(edges[0], edges[-1]),
                ylim=(db[0], db[-1]),
                xlabel="Frequency [Hz]",
                ylabel="Pressure PSD [dB re 1 Pa²/Hz]",
                title=title,
            )
            axis.grid(alpha=0.2)
            return mesh

        mesh = draw(ax)
        if name == "01_nominal":
            ax.legend(fontsize=8, loc="lower left")
        single, sax = plt.subplots(figsize=(12, 7), layout="constrained")
        sm = draw(sax)
        sax.legend(fontsize=8, loc="lower left")
        single.colorbar(sm, ax=sax, label="Occurrence [% per 1-dB bin]")
        single.suptitle(f"YV.{station}.00.BDH — {len(values):,} saved PSD windows")
        single.savefig(out / f"{name}.png", dpi=150)
        plt.close(single)
    fig.colorbar(mesh, ax=axes.ravel().tolist(), label="Occurrence [% per 1-dB bin]")
    fig.suptitle(
        f"{station} pressure PPSD — identical windows and axes; "
        "illustrative, not accepted calibration\n"
        "*Approximate environmental regions, not acceptance or sensor limits. "
        "Brown: in-water reference, not seafloor bounds.\n"
        "Scalar gains extrapolated beyond estimation bands; local-P is legacy "
        "and band-centre approximate (saved smoothing crosses edges)."
    )
    fig.savefig(out / "comparison.png", dpi=150)
    plt.close(fig)
    delta, dax = plt.subplots(figsize=(11, 4), layout="constrained")
    dax.axhline(0, color="0.6", label="Nominal")
    dax.axhline(pressure_gain_db(ray), label="Rayleigh scalar (hypothetical full band)")
    dax.axhline(
        pressure_gain_db(micro),
        ls="--",
        label="Microseism scalar (hypothetical full band)",
    )
    # Draw each local-P band's gain separately: never join across boundaries.
    for band, item in bands.items():
        lo, hi = map(float, band.split("-"))
        dax.plot(
            [lo, hi],
            [pressure_gain_db(item["gain_median"])] * 2,
            color="tab:red",
            lw=2,
            label="Local-P band estimates (legacy diagnostic)"
            if band == next(iter(bands))
            else "_nolegend_",
        )
    dax.set(
        xscale="log",
        xlim=(edges[0], edges[-1]),
        xlabel="Frequency [Hz]",
        ylabel="Illustrative PSD shift [dB]",
        title=f"{station}: gains only, not average PSDs or proof of calibration",
    )
    dax.legend(fontsize=8)
    dax.grid(alpha=0.2)
    delta.savefig(out / "gain_shifts.png", dpi=150)
    plt.close(delta)
    assert all(sha256(path) == digest for path, digest in protected.items())
    manifest = {
        "protected_sha256": protected,
        "source_unchanged": True,
        "windows": len(values),
        "rayleigh_g": ray,
        "microseism_g": micro,
        "local_p_bands": bands,
        "psd_shift_formula": "20*log10(g)",
        "scalar_full_band_extrapolation": True,
        "local_p_legacy_unvalidated": True,
        "local_p_display_bins": int(np.isfinite(local_shift).sum()),
        "local_p_band_center_approximation": True,
        "local_p_exact_frequency_dependent_response_correction": False,
        "raw_data_read_or_modified": False,
        "saved_ppsd_modified": False,
        "reference_sha256": reference_hashes,
        "references": load_references()[0],
        "environmental_regions_are_sensor_or_acceptance_limits": False,
        "illustrated_time_invariance_not_validated": True,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(out)
    if args.show:
        subprocess.run(
            ["open", str(out / "comparison.png"), str(out / "gain_shifts.png")],
            check=True,
        )


if __name__ == "__main__":
    main()
