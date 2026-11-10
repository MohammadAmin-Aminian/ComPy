"""Command line of the automatic inversion: ``compy-invert``.

Example::

    compy-invert RR52_compliance_compy.npz --lat -20.47 --lon 68.11 --crust1 ~/crust1

reads a stage C result (or a CSV with the columns ``freq_hz``, ``compliance``
and ``uncertainty`` in Hz and 1/Pa, with ``--depth``), builds the model of the
station from CRUST1.0, inverts, and writes ``<out>.npz``, ``<out>.json``,
``<out>_profile.csv`` and ``<out>.png`` next to the input. ``--crust1`` is the
folder with the CRUST1.0 files; ``--fetch-crust1`` downloads them there first
(Laske et al. 2013, about 1 MB).
"""

from __future__ import annotations

import argparse
import sys
import tarfile
import urllib.request
from pathlib import Path

import numpy as np

from compy.inversion import auto, dram
from compy.inversion.errors import default_covariance

CRUST1_URL = "https://igppweb.ucsd.edu/~gabi/crust1/crust1.0.tar.gz"
CRUST1_FILES = ("crust1.bnds", "crust1.vp", "crust1.vs", "crust1.rho")


def fetch_crust1(directory):
    """Download the four CRUST1.0 tables into ``directory`` if they are missing."""
    directory = Path(directory).expanduser()
    if all((directory / name).exists() for name in CRUST1_FILES):
        return directory
    directory.mkdir(parents=True, exist_ok=True)
    archive = directory / "crust1.0.tar.gz"
    urllib.request.urlretrieve(CRUST1_URL, archive)  # noqa: S310 (fixed https URL)
    with tarfile.open(archive) as tar:
        for name in CRUST1_FILES:  # only the named data tables, written as plain files
            member = tar.extractfile(name)
            (directory / name).write_bytes(member.read())
    return directory


def read_compliance(path, depth=None, band=None):
    """Frequencies, compliance, relative covariance, depth and a record from a file.

    A stage C ``.npz`` of ComPy 3 carries the measured covariance; an earlier
    one, or a CSV, gets ``errors.default_covariance`` from its uncertainties.
    """
    path = Path(path)
    if path.suffix == ".npz":
        d = dict(np.load(path))
        depth = float(d["depth"]) if depth is None else depth
        lo, hi = d["band"] if band is None else band
        if "cov_rel" in d and band is None:
            keep = d["cov_usable"].astype(bool)
            f = d["cov_freqs"][keep]
            eta = np.interp(f, d["freqs"], d["compliance"])
            return (
                f,
                eta,
                d["cov_rel"][np.ix_(keep, keep)],
                depth,
                "bootstrap of the windows",
            )
        sel = (d["freqs"] >= lo - 1e-12) & (d["freqs"] <= hi + 1e-12)
        f, eta, sigma = d["freqs"][sel], d["compliance"][sel], d["uncertainty"][sel]
    else:
        table = np.genfromtxt(path, delimiter=",", names=True)
        names = table.dtype.names
        col = {key: next(n for n in names if n.startswith(key)) for key in
               ("freq", "compliance", "uncertainty")}  # fmt: skip
        f, eta, sigma = (table[col[k]] for k in ("freq", "compliance", "uncertainty"))
        if "in_band" in names and band is None:
            sel = table["in_band"] > 0
            f, eta, sigma = f[sel], eta[sel], sigma[sel]
        elif band is not None:
            sel = (f >= band[0]) & (f <= band[1])
            f, eta, sigma = f[sel], eta[sel], sigma[sel]
        if depth is None:
            raise SystemExit("--depth is required for a CSV file")
    return (
        f,
        eta,
        default_covariance(f, eta, sigma),
        depth,
        "default structure from uncertainties",
    )


def main(argv=None):
    """Entry point of ``compy-invert``."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    parser.add_argument(
        "compliance", help="stage C .npz or a CSV of the compliance curve"
    )
    parser.add_argument(
        "--lat", type=float, required=True, help="station latitude [deg]"
    )
    parser.add_argument(
        "--lon", type=float, required=True, help="station longitude [deg]"
    )
    parser.add_argument(
        "--crust1", required=True, help="folder with the CRUST1.0 tables"
    )
    parser.add_argument(
        "--fetch-crust1", action="store_true", help="download CRUST1.0 if missing"
    )
    parser.add_argument(
        "--depth", type=float, help="water depth [m] (read from an .npz)"
    )
    parser.add_argument("--band", type=float, nargs=2, metavar=("F_LO", "F_HI"),
                        help="frequencies used [Hz]; default the stored band")  # fmt: skip
    parser.add_argument("--gain", default="calibrated", choices=sorted(auto.GAIN_SIGMA),
                        help="how the pressure gain was obtained (sets its prior sigma)")  # fmt: skip
    parser.add_argument(
        "--gain-sigma", type=float, help="prior sigma of the gain factor"
    )
    parser.add_argument("--response", default="unchecked", choices=sorted(auto.RESPONSE_TILT),
                        help="low-frequency pressure response: checked against the tide or not")  # fmt: skip
    parser.add_argument("--sediment-prior", type=float, nargs=2, metavar=("MEAN", "SIGMA"),
                        help="Gaussian prior on sediment thickness [m], e.g. from a sediment map")  # fmt: skip
    parser.add_argument("--out", help="output stem; default <input>_auto")
    parser.add_argument("--processes", type=int, default=None)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-figure", action="store_true")
    args = parser.parse_args(argv)

    crust1 = (
        fetch_crust1(args.crust1)
        if args.fetch_crust1
        else Path(args.crust1).expanduser()
    )
    f, eta, cov_rel, depth, source = read_compliance(
        args.compliance, args.depth, args.band
    )
    station = Path(args.compliance).stem.split("_compliance")[0]
    data = auto.compliance_data(f, eta, depth, cov_rel=cov_rel, response=args.response,
                                station=station)  # fmt: skip
    column = dram.crust1_column(args.lat, args.lon, crust1)
    overrides = {}
    if args.sediment_prior:
        overrides["sediment_prior"] = tuple(args.sediment_prior)
    model = auto.model_from_crust1(column, **overrides)
    sigma = args.gain_sigma or auto.GAIN_SIGMA[args.gain]
    info = dict(file=str(args.compliance), latitude=args.lat, longitude=args.lon,
                covariance=source, response=args.response, crust1=column)  # fmt: skip
    result = auto.invert(data, model, gain_sigma=sigma, processes=args.processes, seed=args.seed,
                         info=info)  # fmt: skip
    stem = (
        Path(args.out)
        if args.out
        else Path(args.compliance).with_name(station + "_auto")
    )
    result.save(stem)
    if not args.no_figure:
        from compy.plotting.auto import station_figure

        station_figure(result, stem.with_suffix(".png"))
    q = result.quality()
    print(f"{station}: {q['label']}; chi-square {q['chi2']:.1f} for {q['n_data']} data; "
          f"max R-hat {q['max_rhat']:.3f}; gain factor {q['gain']:.3f}")  # fmt: skip
    for reason in q["reasons"]:
        print(f"  note: {reason}")
    prof = result.profile()
    res = prof["z"][prof["resolved"]]
    if len(res):
        print(
            f"  Vs resolved between {res.min():.0f} and {res.max():.0f} m below the seafloor"
        )
    print(
        f"  written: {stem}.npz, .json, _profile.csv"
        + ("" if args.no_figure else ", .png")
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
