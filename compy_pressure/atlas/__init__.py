"""Global seafloor pressure-noise atlas (ComPy 3.2).

From an SDS archive and StationXML to station-balanced noise models of seafloor
pressure gauges (DPG, APG, hydrophone), stratified by sensor, depth, distance to
the coast and to the shelf break, basin, latitude and season.

Modules:
    native      one estimator for native-rate spectra, 0.02 s to about 1260 s
    geography   distance to the coastline and to the shelf break (200 m isobath)
    classify    sensor type and brand, basin, depth and distance classes, seasons
    cube        per-channel window screens, histograms and daily band levels
    screens     channel cascade, level gates, P-Z calibration, day and band screens
    models      station-balanced percentiles and Peterson-style line fits
    figures     standard figures
    band_rule   rule W: two limits per channel and window, from the instrument and the window
    pipeline    the whole build from a channel table
    cli         the ``compy-pressure-atlas`` command
    config      every threshold, with the values of the November 2026 atlas

Developed on the ComPy pressure archive (975 channels, 35 networks) in November
2026; see docs/PRESSURE_ATLAS.md.
"""
