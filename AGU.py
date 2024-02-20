#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 11:48:48 2023

@author: mohammadamin
"""


from obspy import read
from obspy.clients.fdsn import Client
import os
import numpy as np
import matplotlib.pyplot as plt
# import ffplots as fp
from obspy import UTCDateTime

client = Client("RESIF")
net = "YV"
sta = "RR52"

inv = client.get_stations(
    network=net,
    station=sta,
    channel="BHZ",
    location="*",
    level="response")
A = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/')
A.sort()
stream = read()
stream.clear()

# for i in range(0,len(A)):
#     stream = stream + read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/'+A[i])
# stream.merge(fill_value='interpolate')
#%% Loading rotated stream
A = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Rotated/')
A.sort()

rotated_stream = read()
rotated_stream.clear()

i = 5

rotated_stream = rotated_stream + read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Rotated/'+A[i])

rotated_stream.trim(starttime=UTCDateTime("2013-04-2T03:00:00"),endtime=UTCDateTime("2013-04-2T04:00:00"))

rotated_stream.filter("bandpass", freqmin=0.009, freqmax=0.011)


time = np.arange(len(rotated_stream.select(channel="BHZ")[0].data)) / rotated_stream.select(channel="BHZ")[0].stats.sampling_rate

plt.figure(dpi=300,figsize=(60,8))
plt.title("YV.RR52 1-Hour -- Filtered [0.009 Hz to 0.011 Hz]")
plt.plot(time,rotated_stream.select(channel="BDH")[0].data / np.max(rotated_stream.select(channel="BDH")[0].data),label="Pressure",linewidth=5)
plt.plot(time,rotated_stream.select(channel="BHZ")[0].data / np.max(rotated_stream.select(channel="BHZ")[0].data),label="Vertical Displacement",linewidth=5)
plt.axvspan(0, 800, color='lightgreen', alpha=0.5)
plt.axvspan(800, 2400, color='lightcoral', alpha=0.5)
plt.axvspan(2400, 3600, color='lightgreen', alpha=0.5)
plt.xlabel("Time [s]")
plt.ylabel("Normalized Amplitude")
plt.legend(loc='upper right',fontsize=25)
plt.xlim([0,3600])
plt.tight_layout()
#%%
