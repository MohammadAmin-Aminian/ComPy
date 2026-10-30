"""Stage B: Occam's inversion (Constable et al. 1987) as an independent check.

The minimum-structure method of Crawford et al. (1991): Vs of 20 layers
growing by 1.2 from about 50 m down to the CRUST1.0 Moho is inverted for the
smoothest model whose chi-square on the data is at most N. The starting Vs,
Vp and density are those of CRUST1.0 at the station, as Crawford started from
refraction models; Vp and density stay fixed. The errors are the marginal
standard deviations of the Stage A covariance, so the gain error is treated
as independent per frequency.
"""

import json
import sys

import numpy as np

from a_invert import ROOT, build_model, load_data

import compy_inversion as ci  # noqa: E402  (path set by a_invert)
import inv_compy as inv  # noqa: E402


def main():
    station = sys.argv[1]
    data, _ = load_data(station)
    model, column = build_model(station)
    ref = ci.crust1_model(column)
    thick = 50.0 * 1.2 ** np.arange(20)
    thick *= column["moho_depth"] / thick.sum()
    mid = np.cumsum(thick) - thick / 2
    tops = np.concatenate([[0.0], np.cumsum(ref[:-1, 0])])
    rows = ref[np.searchsorted(tops, mid, side="right") - 1]
    start = np.vstack(
        [np.column_stack([thick, rows[:, 1:3], np.maximum(rows[:, 3], 300.0)]), ref[-1]]
    )
    result = inv.invert_occam(data.eta, data.freqs, data.depth, start, data.sigma)
    m = result["model"]
    out = dict(chi2=result["chi2"], chi2_start=float(np.sum(((data.eta - inv.calc_norm_compliance(
        data.depth, data.freqs, start)) / data.sigma) ** 2)), iterations=result["iterations"],
        thickness=m[:, 0].tolist(), vs=m[:, 3].tolist(), vs_start=start[:, 3].tolist())  # fmt: skip
    (ROOT / "results" / station / f"{station}_occam.json").write_text(
        json.dumps(out, indent=1)
    )
    print(station, {k: out[k] for k in ("chi2_start", "chi2", "iterations")})


if __name__ == "__main__":
    main()
