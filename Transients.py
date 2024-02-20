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
# import tiskitpy as tiskit
import tiskit
from tiskitpy.rptransient import  Transients, PeriodicTransient as PT
# import compy

client = Client("RESIF")
start =UTCDateTime("2012-10-12")
net="YV"
sta="RR52"

invz = client.get_stations(
                           network=net,
                           station=sta,
                           channel="BHZ",
                           location="*",
                           level="response")
inv1 = client.get_stations(
                           network=net,
                           station=sta,
                           channel="BH1",
                           location="*", 
                           level="response")
inv2 = client.get_stations(
                           network=net,
                           station=sta,
                           channel="BH2",
                           location="*",
                           level="response")
invp = client.get_stations(
                           network=net,
                           station=sta,
                           channel="BDH",
                           location="*",
                           level="response")

#%%
import os 
import compy

j = 10

i = 0


client = Client("RESIF")
# start =UTCDateTime("2012-10-12")
net = "YV"
# sta = "RR28"


A = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/')
A.sort()

stream = read()
stream.clear()

transients_stream = read()
transients_stream.clear()


stream =  read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[j])
    
    # stream = stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])
stream.sort()
    # stream[0].stats.channel = "BHZ"
    # stream[3].stats.channel = "BDH"
    
split_streams = compy.cut_stream_with_overlap(stream, 1*24*60*60, overlap = 0)
transients_stream.clear()
        
Rot = tiskit.CleanRotator(split_streams[i], remove_eq=False,filt_band=(0.001, 1))
        # Rot = tiskit.CleanRotator(split_streams[i], remove_eq=False,filt_band=(0.002, 0.03))
    
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
                      'RR52': [PT("1h", 3619.7, 0.05, [-550,155], split_streams[i][0].stats.starttime)]}

rt = Transients(transients[sta])
    
zdata = split_streams[i].select(channel='*Z')[0]
    
eq_spans = tiskit.TimeSpans.from_eqs(zdata.stats.starttime, zdata.stats.endtime, minmag=5.5,days_per_magnitude = 0.5)
    
rt.calc_timing(zdata, eq_spans)

rt.calc_transients(zdata, eq_spans, plot=False)
    
cleaned = rt.remove_transients(zdata, plot=False, match=False, prep_filter=False)
