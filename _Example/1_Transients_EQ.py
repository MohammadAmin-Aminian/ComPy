#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 14:41:23 2023

@author: mohammadamin
"""

import matplotlib.pyplot as plt
import scipy
from obspy import UTCDateTime,read
from obspy.clients.fdsn import Client
import numpy as np
import obspy
import tiskit
from tiskitpy.rptransient import  Transients, PeriodicTransient as PT

sta = "RR52"

import os 
import compy
i = 1

# A = os.listdir("/Users/mohammadamin/Desktop/Data/Z3/decimated")
A = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+"/Decimated/")
A.sort()

stream = read()
stream.clear()

transients_stream = read()
transients_stream.clear()

# stream =  read("/Users/mohammadamin/Desktop/Data/Z3/decimated/"+A[i])
    
stream = stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])
stream.sort()

    
split_streams = compy.cut_stream_with_overlap(stream, 7*24*60*60, overlap = 0)
transients_stream.clear()

i = 1 
    
Rot = tiskit.CleanRotator(split_streams[i], remove_eq=False,filt_band=(0.001, 0.5))

split_streams[i] = Rot.apply(split_streams[i],horiz_too=True)
    
split_streams[i].sort(['channel', 'starttime'])
        
        
transients = {'RR28': [PT("1h", 3620.3, 0.05, [-820, 330], stream[0].stats.starttime)],
                      'RR29': [PT("1h", 3620.3, 0.05, [-290, 135], stream[0].stats.starttime)],
                      'RR31': [PT("1h", 3620.3, 0.05, [-290, 135], stream[0].stats.starttime)],
                      'RR34': [PT("1h", 3620.3, 0.05, [-290, 135], stream[0].stats.starttime)],
                      'RR36': [PT("1h", 3620.3, 0.05, [-290, 135], stream[0].stats.starttime)],
                      'RR38': [PT("1h", 3620.0, 0.05, [-235,180], stream[0].stats.starttime)],
                      'RR40': [PT("1h", 3620.3, 0.05, [-290, 135], stream[0].stats.starttime)],
                      'RR50': [PT("1h", 3620.3, 0.05, [-290, 300], stream[0].stats.starttime)],
                      'RR52': [PT("1h", 3619.7, 0.05, [-600,210], split_streams[i][0].stats.starttime)],
                      'A419A': [PT("1h", 3620.5, 0.05, [-260,140], split_streams[i][0].stats.starttime)]
                      }

rt = Transients(transients[sta])
    
zdata = split_streams[i].select(channel='*Z')[0]

# Generate timespans to avoid because of earthquakes

eq_spans = tiskit.TimeSpans.from_eqs(zdata.stats.starttime, 
                                     zdata.stats.endtime, 
                                     minmag=5,
                                     days_per_magnitude = 0.5,
                                     save_eq_file=False)
ss = eq_spans.zero(stream)
plt.figure(dpi=300,figsize=[30,10])
plt.plot(stream[3].data[1000000:2000000],color="gray",linewidth=5,label="Raw")
plt.plot(ss[3].data[1000000:2000000],color="black",linewidth=5,label="Raw - Global EQ")
plt.xlabel("Time Sample")
plt.ylabel("Amplitude")
plt.legend(loc="upper right")
plt.show()

# calcualates and stores a list of PeriodicTransents

rt.calc_timing(zdata, eq_spans)
# Calculate transient time parameters

rt.calc_transients(zdata, eq_spans, plot=False)
# Remove transient from data

cleaned = rt.remove_transients(zdata, plot=False, match=False, prep_filter=False)


plt.figure(dpi=300,figsize=[30,10])
plt.title("Residuals of One Day")
plt.plot((cleaned.data - zdata.data)[0:181440],linewidth=5)
plt.xlabel("Time Sample")
plt.ylabel("Amplitude")
plt.show()









































