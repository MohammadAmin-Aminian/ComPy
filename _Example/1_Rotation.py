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

# Read The downloaded Stream , instead of "Path" Import the destination of your file. or if you have downloaded it using prevous example, skip this line

stream_decim = read("Path")


# Rotates the Data , and Remove coherence Noise to minimize Tilt effect on Hourly Basis time_window = 1 hour.
# it also remove instrument response

rotated_stream,azimuth,angle,variance = compy.Rotate(stream_decim,time_window = 6)

fp.psd(rotated_stream)

fp.coh(rotated_stream)

fp.coh_h(rotated_stream)
