#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:54:22 2024

@author: mohammadamin
"""
from obspy.clients.fdsn import Client
from obspy import UTCDateTime,read
import ffplot as fp
from tiskitpy import Decimator
import compy
compy.plt_params()
import Pressure_calibration as DPG


stream = read("Path")

gain_factor = DPG.calculate_spectral_ratio(stream,mag = 7 ,f_min=0.03,f_max=0.07,plot_condition = True)
