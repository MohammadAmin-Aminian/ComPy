"""Time-frequency diagnostics of one station from cleaned (or raw physical) MiniSEED.

The figures of ``compy.plotting.diagnostics`` that need waveforms:
spectrograms of Z, P, H1 and H2, coherograms (Z-P, Z-H1, Z-H2, P-H1, P-H2)
and the tilt azimuth/admittance per window, over the whole deployment::

    python examples/compy3_waveform_diagnostics.py ~/runs/yv/work/RR52/clean \\
        --station RR52 --band 0.0073 0.0171 --output figures/RR52

The input folder holds the stage B segments (``<STA>_seg*.mseed``) or the
stage A ``phys/`` day files. Each segment is cut into ``--window`` second
windows; windows with gaps are skipped, so nothing is interpolated.
Needs ObsPy. Writes PNGs and ``<STA>_window_spectra.npz`` (band-averaged
coherence, admittance and tilt per window) for ``compy.anomaly``.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from compy import anomaly, timefreq  # noqa: E402
from compy.plotting import diagnostics as dg  # noqa: E402

NAMES = {"Z": "*Z", "P": "*H", "1": "*1", "2": "*2"}


def segments(folder, station):
    """Yield (t0, fs, {name: array}) for each gap-free stretch common to all channels."""
    from obspy import Stream, read

    files = sorted(Path(folder).glob(f"{station}*.mseed"))
    if not files:
        raise SystemExit(f"no MiniSEED files for {station} in {folder}")
    st = Stream()
    for path in files:
        st += read(str(path))
    st.merge(method=0)
    st = st.split()
    picked = {k: st.select(channel=v) for k, v in NAMES.items()}
    for tr_z in picked["Z"]:
        t0, t1 = tr_z.stats.starttime, tr_z.stats.endtime
        parts = {}
        for name, traces in picked.items():
            for tr in traces:
                a, b = max(t0, tr.stats.starttime), min(t1, tr.stats.endtime)
                if b - a > 0.5 * (t1 - t0):
                    t0, t1 = a, b
                    parts[name] = tr
                    break
        if len(parts) < 4:
            continue
        fs = tr_z.stats.sampling_rate
        data = {
            k: tr.slice(t0, t1, nearest_sample=False).data.astype(float)
            for k, tr in parts.items()
        }
        n = min(len(v) for v in data.values())
        yield float(t0.timestamp), fs, {k: v[:n] for k, v in data.items()}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("folder", type=Path)
    parser.add_argument("--station", required=True)
    parser.add_argument(
        "--band",
        type=float,
        nargs=2,
        default=(0.007, 0.017),
        help="compliance band [Hz]",
    )
    parser.add_argument("--tilt-band", type=float, nargs=2, default=(0.005, 0.035))
    parser.add_argument(
        "--window", type=float, default=7200.0, help="window length [s]"
    )
    parser.add_argument("--nseg", type=int, default=1024)
    parser.add_argument("--output", type=Path, default=Path("."))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    parts = []
    for t0, fs, data in segments(args.folder, args.station):
        if min(len(v) for v in data.values()) < args.window * fs:
            continue
        parts.append(timefreq.window_spectra(data, fs, args.window, args.nseg, t0=t0))
    if not parts:
        raise SystemExit("no complete windows")
    spec = {k: np.concatenate([p[k] for p in parts]) for k in parts[0] if k != "f"}
    spec["f"] = parts[0]["f"]
    t, f = spec["times"], spec["f"]
    sta = args.station
    psd = {"Z": spec["ZZ"], "P": spec["PP"], "H1": spec["11"], "H2": spec["22"]}
    fig = dg.plot_spectrogram(t, f, psd, band=args.band, title=f"{sta}: spectrograms")
    fig.savefig(args.output / f"{sta}_spectrogram.png", dpi=130)
    pairs = {
        "Z–P": ("P", "Z"),
        "Z–H1": ("1", "Z"),
        "Z–H2": ("2", "Z"),
        "P–H1": ("P", "1"),
        "P–H2": ("P", "2"),
    }
    coh = {k: timefreq.coherence(spec, *v) for k, v in pairs.items()}
    fig = dg.plot_coherogram(
        t, f, coh, band=args.band, threshold=0.9, title=f"{sta}: coherograms"
    )
    fig.savefig(args.output / f"{sta}_coherogram.png", dpi=130)
    tilt = timefreq.tilt_direction(spec, tuple(args.tilt_band))
    fig = dg.plot_tilt(t, tilt["azimuth"], tilt["admittance"], strength=tilt["coherence"],
                       angle_label="Admittance |S_HZ|/S_HH", strength_label="Z–H coherence",
                       title=f"{sta}: tilt direction and strength")  # fmt: skip
    fig.savefig(args.output / f"{sta}_tilt.png", dpi=130)
    plt.close("all")
    band = {k: timefreq.band_average(f, v, args.band) for k, v in coh.items()}
    adm = timefreq.band_average(f, timefreq.admittance(spec, "P", "Z"), args.band)
    np.savez(args.output / f"{sta}_window_spectra.npz", times=t, admittance_zp=adm,
             tilt_azimuth=tilt["azimuth"], tilt_admittance=tilt["admittance"],
             tilt_coherence=tilt["coherence"], **{k.replace("–", "_"): v for k, v in band.items()})  # fmt: skip
    ok = np.isfinite(adm)
    for name, series in (
        ("Z-P admittance", np.log10(adm[ok])),
        ("tilt admittance", tilt["admittance"][ok]),
    ):
        for c in anomaly.changepoints(t[ok], series, min_size=50):
            print(
                f"{sta} {name}: step {anomaly._iso(c.time)} shift {c.shift:+.3g} (σ {c.significance:.1f})"
            )
        hits = {
            k: v
            for k, v in anomaly.periodicity(t[ok], series).items()
            if v["p_tail"] < 0.01
        }
        if hits:
            print(f"{sta} {name}: periodic at {', '.join(hits)}")


if __name__ == "__main__":
    main()
