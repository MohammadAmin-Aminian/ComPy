"""Stage F: ComPy 2.1 vs 2.2 legacy estimator on identical cleaned RR52 data.

Both versions run ``estimate(..., beta=True)`` (the historical beta gates,
2-h windows stepped by 5 min, nseg=4096) on the same Stage B output and the
same gain, so differences come only from the code fixes D1 (gravity term not
scaled by coherence) and D2 (uncertainty over all windows). Each version is
run in its own subprocess so the two module sets never mix.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable
CODE = r"""
import sys, json, glob, numpy as np
sys.path.insert(0, sys.argv[1])
from obspy import read
from compy_processing import estimate
st = read(sorted(glob.glob(sys.argv[2]))[0])
st = st.select(channel="BHZ") + st.select(channel="BDH")
for tr in st:
    tr.data = tr.data.astype(float)
curves, coh, acc, fc, f, scatter, theory, _, _ = estimate(
    st, float(sys.argv[3]), float(sys.argv[4]), 2, 0.007, 0.019, True, 4096)
json.dump({"fc": fc.tolist(), "mean": np.mean(curves, 0).tolist(),
           "theory": np.asarray(theory).tolist(), "std": scatter.tolist(),
           "n": len(acc)}, open(sys.argv[5], "w"))
"""


def main():
    station, depth = "RR52", 2918.0
    summary = json.loads(
        (ROOT / "results" / station / f"{station}_summary.json").read_text()
    )
    gain = summary["gains"]["bruitfm"]
    seg = str(ROOT / "work" / station / "clean" / f"{station}_seg*.mseed")
    out = {}
    for label, code_dir in (("v2.1", ROOT / "work" / "compy_v2.1"),
                            ("v2.2", ROOT.parent / "ComPy")):  # fmt: skip
        path = ROOT / "results" / station / f"legacy_beta_{label}.json"
        subprocess.run([PY, "-c", CODE, str(code_dir), seg, str(depth), str(gain), str(path)],
                       check=True)  # fmt: skip
        out[label] = json.loads(path.read_text())
    fc = np.array(out["v2.1"]["fc"])
    old, new = np.array(out["v2.1"]["mean"]), np.array(out["v2.2"]["mean"])
    band = (fc >= 0.007) & (fc <= 0.019)
    report = {
        "windows": out["v2.2"]["n"],
        "median_change_percent_in_band": float(np.median((new / old - 1)[band] * 100)),
        "max_change_percent_in_band": float(np.max(np.abs(new / old - 1)[band] * 100)),
        "change_percent_by_freq": {
            f"{f:.4f}": float(c) for f, c in zip(fc[band], (new / old - 1)[band] * 100)
        },  # fmt: skip
        "theory_v2.1_median_rel": float(
            np.median((np.array(out["v2.1"]["theory"]) / old)[band])
        ),
        "theory_v2.2_median_rel": float(
            np.median((np.array(out["v2.2"]["theory"]) / new)[band])
        ),
    }
    (ROOT / "results" / station / "legacy_v21_vs_v22.json").write_text(
        json.dumps(report, indent=1)
    )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    axes[0].plot(
        fc[band], old[band], "o-", color="#2a78d6", label="ComPy 2.1 (γ × gravity term)"
    )
    axes[0].plot(fc[band], new[band], "s-", color="#eb6834", label="ComPy 2.2")
    axes[0].set(xlabel="Frequency [Hz]", ylabel="Normalized compliance η [1/Pa]",
                title=f"RR52 legacy beta estimator, {report['windows']} windows")  # fmt: skip
    axes[0].legend(frameon=False)
    axes[1].plot(fc[band], (new / old - 1)[band] * 100, "o-", color="#0b0b0b")
    axes[1].axhline(0, color="#9a9994", lw=1)
    axes[1].set(
        xlabel="Frequency [Hz]", ylabel="2.2 / 2.1 − 1 [%]", title="Effect of fix D1"
    )
    fig.savefig(
        ROOT / "figures" / "RR52_legacy_v21_vs_v22.png", dpi=110, bbox_inches="tight"
    )
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
