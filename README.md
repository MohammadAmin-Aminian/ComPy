# ComPy: Seafloor Compliance Analysis Tool
<p align="center">
  <img src="_Images/ComPy.png" width="225">
</p>


# Overview

ComPy is a specialized software tool designed for the automated processing and analysis of seafloor compliance data. It enhances the precision of subsurface shear velocity models using Broadband Ocean-Bottom Stations data from the Indian Ocean. ComPy is tailored to manage extensive data processing steps, ensuring high resolution and accuracy in geological analysis.

# Features

Automation of data pre-processing steps including glitch removal, tilt effect minimization, and pressure gauge calibration.
Implementation of the Metropolis-Hastings algorithm for robust depth-velocity inversion.
Utilization of advanced signal processing techniques to handle seafloor compliance data.

# Installation

Before installing ComPy, ensure you have Python and the necessary packages installed. ComPy requires Python 3.x.

# Clone the repository
git clone https://github.com/your-repository/ComPy.git

# Navigate to the ComPy directory
cd ComPy

# Install required Python packages
pip install numpy matplotlib scipy obspy tiskitpy

# Usage

Here's how you can use ComPy to process your seafloor compliance data:

import compy

rotated_stream,azimuth,angle,variance = compy.Rotate(stream_decim,time_window = 1)

**Function Overview:**
The `compy.Rotate` function rotates seismic data to minimize tilt effects and removes coherence noise, enhancing data accuracy for compliance analysis. The default processing window is set to 1 hour but can be adjusted as needed.

- **rotated_stream**: The seismic data stream after rotation and noise removal.
- **azimuth**: The direction of the rotation applied to correct the tilt in degrees.
- **angle**: The angle of tilt correction applied to the seismic data.
- **variance**: The reducted variance ratio (After/Before), indicating the effectiveness of noise reduction.


<p align="center">
  <img src="_Images/RR52_Tilt.png" width="700">
</p>

# Plotting Functions

## `fp.coherogram_spectrogram_alpha(rotated_stream)`

This function plots spectrograms that feature average values within the compliance frequency band. The output includes several panels:

- **(a) Calibrated Pressure:** Shows pressure data adjusted for calibration errors.
- **(b) Median Values of Pressure:** Displays the median values of the calibrated pressure data.
- **(c) Corrected Vertical Acceleration:** Presents the vertical acceleration data after corrections.
- **(d) Median Values of Vertical Acceleration:** Shows the median values of the corrected vertical acceleration data.
- **(e) Coherogram:** Plots the coherogram between the calibrated pressure and the corrected vertical acceleration.
- **(f) Average Values:** Illustrates the average values across the dataset. 

Black dashed lines indicate the frequency limits of the compliance band. Green shaded areas in panel (f) highlight selected time windows with coherency exceeding the coherence threshold, marked by a red dashed line set at 0.8. This visualization aids in identifying significant patterns and anomalies in the data.
<p align="center">
  <img src="_Images/RR52_window_selection.png" width="700">
</p>

# Contributing

We welcome contributions from the community. Please review CONTRIBUTING.md for guidelines on how to submit improvements to ComPy.

# License

This project is licensed under the MIT License - see the LICENSE file for details.

# Citation
_Inference of Shallow Subsurface Structures of the Indian Ocean Derived from Compliance Function Analysis
# Acknowledgments

This tool was developed at the Institut de Physique du Globe de Paris and funded by the SPIN project, an Innovative Training Network (ITN) supported by the European Commission under the Horizon 2020 Marie Skłodowska-Curie Actions (MSCA). We extend our gratitude to all contributors and collaborators who have made this project possible.
<p align="center">
  <img src="_Images/H2020_acknowledgment.png" width="300" style="display: inline-block;">
  <img src="_Images/IPGP_UPC.png" width="300" style="display: inline-block;">
</p>



