"""Stage H: H1 versus H2 transfer-function estimates (noise on Z or on P).

The 2.2 compliance uses H1 = |S_pz| / S_pp, unbiased when all noise is on
the vertical. If the pressure channel also carries noise, H1 is biased low
by the coherence squared (Bendat & Piersol 2010; Crawford & Singh 2008),
while H2 = S_zz / |S_pz| is unbiased for noise on P alone. The truth lies
between the two, and H2/H1 = 1/coherence**2. Both are gravity corrected the
same way. Uses the stored stacked spectra of every compliance file.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COMPLIANCE = ROOT.parent / "compliance_yv_2026-10-28"
sys.path.insert(0, str(ROOT.parent / "ComPy"))

from compy_processing import gravity_corrected_compliance  # noqa: E402

STATIONS = ("RR28", "RR29", "RR34", "RR36", "RR38", "RR40", "RR50", "RR52")
EPOCH = {"RR36", "RR38", "RR40", "RR50"}


def estimators(path):
    d = np.load(path)
    f, lo, hi = d["all_freqs"], *d["band"]
    keep = (f >= lo - 1e-12) & (f <= hi + 1e-12)
    szz, spp, spz = d["S_zz"][keep], d["S_pp"][keep], np.abs(d["S_pz"][keep])
    depth = float(d["depth"])
    h1 = gravity_corrected_compliance(f[keep], depth, spz / spp, np.sqrt(spp))
    h2 = gravity_corrected_compliance(f[keep], depth, szz / spz, np.sqrt(spp))
    coh = spz / np.sqrt(szz * spp)
    return dict(freqs=f[keep].tolist(), h1=h1.tolist(), h2=h2.tolist(), coherence=coh.tolist(),
                h2_over_h1=(h2 / h1).tolist(), stored=np.interp(f[keep], d["freqs"], d["compliance"]).tolist())  # fmt: skip


def main():
    out = {}
    for s in STATIONS:
        sub = "epoch/" if s in EPOCH else ""
        out[s] = {
            "stable": estimators(
                COMPLIANCE / "results" / s / f"{sub}{s}_compliance_compy.npz"
            )
        }
        after = ROOT / "compliance_after" / s / f"{s}_compliance_after.npz"
        if after.exists():
            out[s]["after"] = estimators(after)
        for k, v in out[s].items():
            r = np.array(v["h2_over_h1"])
            check = np.max(np.abs(np.array(v["h1"]) / np.array(v["stored"]) - 1))
            print(
                f"{s} {k:6s} H2/H1 {r.min():.3f}-{r.max():.3f}  (lowest f {r[0]:.3f})  "
                f"min coherence {min(v['coherence']):.3f}  H1 vs stored {check:.1e}"
            )
    (ROOT / "results" / "tables" / "h1_h2.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
