#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:27:47 2024

@author: mohammadamin
"""
from obspy.clients.fdsn import Client
from obspy import UTCDateTime,read
import ffplot as fp
from tiskitpy import Decimator
import compy
compy.plt_params()

# Read the downloaded stream. Replace "Path" with the actual file path where your file is stored. 
# If the file was downloaded using the previous example, you can skip this line.

stream_decim = read("Path")


# Rotate the data and remove coherence noise to minimize tilt effects. Set the time window to 1 hour.

rotated_stream,azimuth,angle,variance = compy.Rotate(stream_decim,time_window = 24)

fp.psd(rotated_stream)

fp.coh(rotated_stream)

fp.coh_h(rotated_stream)