#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 11:18:01 2023

@author: mohammadamin
"""
from obspy import read
from obspy.clients.fdsn import Client
import ffplot as fp
from obspy import UTCDateTime
import obspy
# import tiskit
# import ffplot as fp
import os
import numpy as np
import matplotlib.pyplot as plt
# import com_forward as cm
# nhnm = obspy.signal.spectral_estimation.get_nhnm()
# nlnm = obspy.signal.spectral_estimation.get_nlnm()
import tiskit
import scipy
import scipy.stats
# import math as M
import obspy.signal
nhnm = obspy.signal.spectral_estimation.get_nhnm()
nlnm = obspy.signal.spectral_estimation.get_nlnm()
import compy
compy.plt_params()
import pickle

angle = []
azimuth = []

stations = ["RR28","RR29","RR34","RR36","RR38","RR40","RR50","RR52"]

for i in range(0,len(stations)):
    sta = stations[i]
    file_path = "/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Compliance/'



    with open(file_path+f'rotation_hourly_{sta}.pkl', 'rb') as file:
        loaded_rotation_container = pickle.load(file)
    

    stream_stats =loaded_rotation_container['stream stats']


    #RR52
    # endtime = UTCDateTime("2013-09-15T00")
    # #RR38
    # starttime = UTCDateTime("2012-12-15T00")
    
    endtime = stream_stats.endtime
    starttime = stream_stats.starttime
    

    aa = (endtime - starttime) // (3600)
    # aa = -1 
    
    angle.append(loaded_rotation_container["angle"][0:int(aa)])
    
    azimuth.append(loaded_rotation_container["azimuth"][0:int(aa)])



#Auto Corrolation
auto_angle= []
auto_Azimuth= []
stations = ["RR28","RR29","RR34","RR36","RR38","RR40","RR50","RR52"]
stations1 = ["RR28","RR29","RR34","RR36","RR38","RR40","RR50","RR52"]

jj = 60
stations1.pop(jj)
for i in range(0,len(angle)):
    if i ==jj:
        pass
    else:
        auto_angle.append(np.correlate(angle[i],angle[i] ,mode='full'))
        auto_Azimuth.append(np.correlate(azimuth[i],azimuth[i] ,mode='full'))


plt.figure(dpi=300, figsize=(25, 20))
plt.subplot(211)
for i in range(0,len(angle)):
    # plt.title("Correlation Angle " + str(stations[jj]))
    plt.title("Autoc orrelation [Incident Angle]")
    plt.plot(auto_angle[i][auto_angle[i].size//2:],linewidth=3,label=str(stations1[i]),linestyle="dashed")
    plt.grid(True)
    plt.legend(loc="upper right",fontsize=30)
    plt.xlabel("Time Lag [Hour]")
    plt.xlim([0,48])
    # plt.ylim([-50000,100000])
    plt.ylim([-1000000,8000000])
    plt.grid(True)

plt.subplot(212)
for i in range(0,len(angle)):
    # plt.title("Correlation Angle " + str(stations[jj]))
    plt.title("Auto correlation [Azimuth]")
    plt.plot(auto_Azimuth[i][auto_Azimuth[i].size//2:],linewidth=3,label=str(stations1[i]),linestyle="dashed")
    plt.grid(True)
    plt.legend(loc="upper right",fontsize=30)
    plt.xlabel("Time Lag [Hour]")
    plt.xlim([0,48])
    # plt.ylim([-50000,100000])
    plt.grid(True)
    # plt.yscale('log')
    plt.ylim([-1000000,25000000])
    
plt.tight_layout()






