# ComPy 2.0 validation record

## Scope

Reviewed the four original Python modules, example scripts, existing regression
tests and README at base commit `57182e5383b2dc5cc8d83f16b2ceb570b32744dd`.
Core numerical and processing paths were refactored and tested; the original
station-oriented publication layouts remain compatibility utilities.

## Scientific sources checked

- Aminian et al. (2025), GJI, DOI: https://doi.org/10.1093/gji/ggaf253.
  Full article reviewed from the coauthor's IPGP copy:
  https://www.ipgp.fr/~stutz/2025-aminian-crawford-stutzmann-al-GJI.pdf.
  Equations 6–8 identify displacement/pressure compliance, amplitude coherence
  and its error estimate. Section 4.1.6 describes gravity acceleration terms;
  Section 4.3 describes squared-derivative roughness and fixed Vp/density.
- OBStools compliance implementation:
  https://github.com/nfsi-canada/OBStools/blob/master/obstools/comply/classes.py.
  Reviewed the original approximate gravity-wave solver and spectral workflow.
  ComPy v2 instead checks its bounded solver directly against the dispersion
  equation. This review is not an independent elastic-propagator benchmark.
- TiSKitPy 2.3.1 documentation and installed source, including `CleanRotator`,
  `DataCleaner`, `TimeSpans.from_eqs`, and rotation application:
  https://tiskitpy.readthedocs.io/latest/.
- NPL's Mackenzie equation reference:
  https://resource.npl.co.uk/acoustics/techguides/soundseawater/underlying-phys.html.
  Used to verify polynomial coefficients and signs.

These are targeted methodological references, not an exhaustive literature survey.

## Numerical checks

The automated tests include:

- Dispersion-equation residuals for DC and positive frequencies across water depths.
- An analytical homogeneous half-space static compliance limit, `(1-nu)/mu`.
- Invariance under splitting a homogeneous finite layer into identical sublayers
  (relative tolerance 1e-8 accommodates propagator cancellation).
- Known linear Z/P transfer, coherence, and gravity-corrected compliance.
- Single-event synthetic pressure calibration, arbitrary gain scaling and disba smoke test.
- Saved inversion predictions compared with the corresponding accepted model;
  fixed half-space, conserved finite depth, evaluated initial state, RNG isolation.
- Exact sample coverage, overlap, copies and MiniSEED read/write round trips.
- Invalid arguments, zero coherence, likelihood underflow and contextual failures.
- Real TiSKitPy rotation and horizontal-noise cleaning integration with a synthetic tilted signal.

The synthetic event test bypasses instrument-response removal because its input
is already in physical units. It does not validate a particular StationXML response.

## Local environment

Python 3.12; NumPy 2.3.5; ObsPy 1.5.1; TiSKitPy 2.3.1; disba 0.7.0.
Local verification: 42 regression tests passed; Ruff checks and formatting passed;
source distribution and wheel built successfully; the offline synthetic example ran. CI is configured
for Python 3.10–3.13; unobserved CI runs are not represented as locally verified.

## Scientific limits and behavior differences

- No RHUM-RUM waveform dataset/response inventory was supplied for a complete
  version-1/version-2 station regression or paper reproduction.
- The elastic propagator is restricted to finite real-valued results; fluid
  sediment layers with Vs=0 are not supported by this solid-layer interface.
- Calibration depends on usable response metadata and genuinely suitable
  teleseismic windows. The strongest-amplitude window heuristic needs inspection.
- Mean-Welch spectra, exact wave dispersion and corrected gravity conversion
  alter derived compliance relative to the original code.
- The v2 sampler uses a standard summed Gaussian data term; the paper describes
  averaged chi-square. Retune regularization rather than reuse alpha blindly.
- This sampler does not sample gauge gain or enforce monotonic Vs. Observation
  covariance, effective independent counts and calibration errors need scientific
  treatment appropriate to the dataset. Acceptance alone does not prove convergence.
- Large depth-profile arrays remain optional for historical tuple compatibility.
- Deployment-specific publication layouts are not covered by the core scientific
  tests. They require their documented data containers and adequate accepted samples.

Tests provide evidence for checked properties, not a guarantee of zero defects.
