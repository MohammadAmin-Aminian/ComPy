#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 11:59:46 2024

@author: Mohammad-Amin AMINIAN


Institut de Physique du Globe de Paris

"""

# Download Data from Broadband Ocean Bottom Station 

from obspy.clients.fdsn import Client
from obspy import UTCDateTime,read
import ffplot as fp
from tiskitpy import Decimator
import compy
compy.plt_params()

client = Client("RESIF")
net = "4G"
sta = "AZBBA"
start_time = "2008-03-01T00:00:00"
end_time = "2008-03-10T00:00:00"


inv = client.get_stations(
        network=net,
        station=sta,
        channel="*",
        location="*",
        level="response")


starttime = UTCDateTime(start_time)
endtime = UTCDateTime(end_time)
stream = client.get_waveforms(network = net, 
                              station = sta , 
                              location = "*", 
                              channel = "*", 
                              starttime = starttime, 
                              endtime = endtime,
                              attach_response= True)

stream.merge(method=1)
# Downsample the data in 3 steps using FIR filter
decim = Decimator([5, 5, 2])
stream_decim = decim.decimate(stream)
stream_decim.remove_response()

# Plots Coherence Function Between Cahnnels
fp.coh(stream_decim)

# Plots Power Spectrum Density Function of various Cahnnels
fp.psd(stream_decim,nseg=2**12)

# Spectrogram and Coherogram (Coherence Over Time)
fp.coherogram_spectrogram_alpha(stream_decim,nesg=2**11)





































