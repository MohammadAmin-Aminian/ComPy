"""Stage G: compliance after the pressure/vertical sensitivity jump.

For RR36, RR38, RR40 and RR50 the P/Z sensitivity changes during the
deployment (compliance run, Stage H). The 2026-10-28 inversion used only
the stable epoch. Here the windows after the jump date are stacked with the
median Rayleigh-wave gain of the events after the jump (ComPy calibration,
Stage D), using the stable epoch's band and the 0.9 coherence threshold so
the two epochs can be compared and inverted together. Output goes to
compliance_after/<STA>/ in the same format as the 2.2 compliance files.
"""

import json
import sys
from pathlib import Path

import numpy as np
from obspy import UTCDateTime

ROOT = Path(__file__).resolve().parents[1]
COMPLIANCE = ROOT.parent / "compliance_yv_2026-10-28"
sys.path.insert(0, str(COMPLIANCE / "scripts"))
sys.path.insert(0, str(ROOT.parent / "ComPy"))

import c_compliance as C  # noqa: E402
from compy_processing import SelectionCriteria, stacked_compliance  # noqa: E402

OUT = ROOT / "compliance_after"


def main():
    station = sys.argv[1]
    step = json.loads(
        (COMPLIANCE / "results" / "tables" / "stability_steps.json").read_text()
    )[station]
    jump = UTCDateTime(step["step_time"])
    cal = json.loads((C.WORK / station / "calibration" / "compy_gain.json").read_text())
    after = [
        (t, g)
        for t, g in zip(cal["events_used"], cal["event_gains"])
        if UTCDateTime(t) >= jump
    ]
    if not after:
        raise SystemExit(f"{station}: no Rayleigh calibration event after {jump}")
    gain = float(np.median([g for _, g in after]))
    stable = np.load(
        COMPLIANCE / "results" / station / "epoch" / f"{station}_compliance_compy.npz"
    )
    band = tuple(stable["band"])
    depth = float(stable["depth"])  # station depth from the stable-epoch file
    windows, dropped = C.load_windows(station)
    windows = [w for w in windows if w[0].stats.starttime >= jump]
    criteria = SelectionCriteria(0.9, (10.0, 50.0), (-170.0, -150.0), "gji2025_0.9")
    r = stacked_compliance(
        windows, depth, gain, band, criteria, C.NSEG, output_band=C.OUTPUT_BAND
    )
    out = OUT / station
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / f"{station}_compliance_after.npz", freqs=r.freqs, compliance=r.compliance,
             uncertainty=r.uncertainty, coherence=r.coherence, phase_deg=r.phase_deg,
             window_mean=r.window_mean, window_std=r.window_std, band=np.array(r.band),
             all_freqs=r.all_freqs, mean_coherence_all=r.mean_coherence_all,
             S_zz=r.stack["S_zz"], S_pp=r.stack["S_pp"], S_pz=r.stack["S_pz"], gain=gain,
             depth=depth, window_starts=np.array([float(t) for t in r.window_starts]))  # fmt: skip
    in_band = (r.freqs >= band[0] - 1e-12) & (r.freqs <= band[1] + 1e-12)
    ratio = r.compliance[in_band] / np.interp(
        r.freqs[in_band], stable["freqs"], stable["compliance"]
    )
    summary = dict(station=station, jump=str(jump), gain_after=gain, events_after=after,
                   gain_stable=float(stable["gain"]), band_hz=band, n_windows=r.n_windows,
                   windows_available=len(windows), rejected=r.rejected,
                   monthly=C.monthly(windows, r, depth, gain),
                   ratio_after_over_stable=ratio.tolist(),
                   median_ratio=float(np.median(ratio)))  # fmt: skip
    (out / f"{station}_summary.json").write_text(
        json.dumps(summary, indent=1, default=str)
    )
    print(f"{station}: jump {jump.date}, gain {summary['gain_stable']:.3f} -> {gain:.3f} "
          f"({len(after)} events), {r.n_windows} windows, after/stable compliance "
          f"{np.round(ratio, 3).tolist()}", flush=True)  # fmt: skip


if __name__ == "__main__":
    main()
