"""Stage G: independent TiSKitPy estimate on the windows ComPy accepted.

TiSKitPy's SpectralDensity (prolate 4π windows) and ResponseFunctions
(noise on the output channel) give Z/P and its uncertainty; TiSKitPy's own
Compliance._zp_to_ncompl converts to normalized compliance. The same
gravity correction as ComPy is then applied, because TiSKitPy 2.3.1 does
not implement ``correct_gravity_terms``. Pressure is multiplied by the
'bruitfm' gain. Three 819.2-s segments are taken in each accepted 35-min
window (starts 0, 640, 1280.8 s).
"""

import json
import sys
from pathlib import Path

import numpy as np
from obspy import UTCDateTime, read

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "ComPy"))

import tiskitpy as tiskit  # noqa: E402
from compy_processing import gravity_corrected_compliance  # noqa: E402

SEGMENT_S = 819.2
OFFSETS = (0.0, 640.0, 1280.8)


def main():
    import logging

    logging.disable(logging.INFO)
    station = sys.argv[1]
    res = dict(
        np.load(ROOT / "results" / station / f"{station}_compliance_bruitfm.npz")
    )
    gain, depth = float(res["gain"]), float(res["depth"])
    st = read(str(ROOT / "work" / station / "clean" / f"{station}_seg*.mseed"))
    st = st.select(channel="BHZ") + st.select(channel="BDH")
    st.merge(method=0)
    st = st.split()
    for tr in st:
        tr.data = tr.data.astype(np.float64)
        if tr.stats.channel == "BDH":
            tr.data *= gain
    starts = [UTCDateTime(t) + o for t in res["window_starts"] for o in OFFSETS]
    # from_stream needs one contiguous trace per channel; use the longest segment set.
    out = {"freqs": [], "zp": [], "unc": [], "coh": []}
    sds = []
    bounds = sorted({(float(tr.stats.starttime), float(tr.stats.endtime)) for tr in st})
    for a, b in bounds:
        seg_id = (UTCDateTime(a), UTCDateTime(b))
        part = st.slice(*seg_id)
        if len(part) != 2:
            continue
        inside = [t for t in starts if seg_id[0] <= t and t + SEGMENT_S <= seg_id[1]]
        if len(inside) < 10:
            continue
        sds.append((len(inside), tiskit.SpectralDensity.from_stream(
            part, window_s=SEGMENT_S, starttimes=inside)))  # fmt: skip
    n, sd = max(sds, key=lambda x: x[0])
    ids = sd.ids
    p_id = next(i for i in ids if i.endswith("BDH"))
    z_id = next(i for i in ids if i.endswith("BHZ"))
    rfs = tiskit.ResponseFunctions(sd, p_id, [z_id], noise_channel="output")
    f = rfs.freqs
    zp = np.abs(rfs.value(z_id))
    unc_mult = rfs.uncert_mult(z_id)
    k_ratio = tiskit.Compliance._zp_to_ncompl(f, zp, depth, "M")
    eta = gravity_corrected_compliance(
        f, depth, zp, np.sqrt(np.abs(sd.autospect(p_id)))
    )
    band = (f >= res["band"][0]) & (f <= res["band"][1])
    ours = np.interp(f[band], res["freqs"], res["compliance"])
    report = {
        "segments_used": int(n), "band": res["band"].tolist(),
        "median_ratio_tiskit_over_compy": float(np.median(eta[band] / ours)),
        "max_abs_diff_percent": float(np.max(np.abs(eta[band] / ours - 1)) * 100),
        "freqs": f[band].tolist(), "tiskit_eta": eta[band].tolist(),
        "tiskit_eta_no_gravity": k_ratio[band].tolist(),
        "tiskit_rel_uncertainty": unc_mult[band].tolist(),
        "compy_eta": ours.tolist(),
    }  # fmt: skip
    (ROOT / "results" / station / f"{station}_tiskit_check.json").write_text(
        json.dumps(report, indent=1)
    )
    print(station, {k: report[k] for k in ("segments_used", "median_ratio_tiskit_over_compy",
                                           "max_abs_diff_percent")})  # fmt: skip


if __name__ == "__main__":
    main()
