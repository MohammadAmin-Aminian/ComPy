#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 13:44:13 2023

@author: mohammadamin
"""
# %% Inventory Download
# You should find the start time and end time from inventory and the put it in the Restriction function
# It should be automatic in future!
# the problem is that it save the channels seperatly no in one single stream....
# so i have to merge the channels before furthur process
# from obspy.signal.trigger import plot_trigger
# from obspy.signal.trigger import coincidence_trigger
# from obspy.signal.trigger import classic_sta_lta
# import matplotlib.pyplot as plt
# import scipy
from obspy import read
from obspy.clients.fdsn import Client
import ffplot as fp
# import numpy as np
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

client = Client("RESIF")
# start =UTCDateTime("2012-10-12")
net = "YV"
sta = "RR38"

inv = client.get_stations(
        network=net,
        station=sta,
        channel="BH*",
        location="*",
        level="response")
    
invp = client.get_stations(
        network=net,
        station=sta,
        channel="*H",
        location="*",
        level="response")
    
A = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/')

B = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/')

C = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Rotated/')

A.sort()
B.sort()
C.sort()

raw_stream = read()
raw_stream.clear()

transients_stream = read()
transients_stream.clear()

rotated_stream = read()
rotated_stream.clear()

i = 4
i = 2
x = 2 
for i in range(0,len(A)):
    raw_stream = raw_stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])
raw_stream.merge(fill_value='interpolate')

for i in range(0,len(B)):
    transients_stream = transients_stream + read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/'+B[i])
    # stream = stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])
transients_stream.merge(fill_value='interpolate')

for i in range(2,len(C)):
    rotated_stream = rotated_stream + read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Rotated/'+C[i])
    # stream = stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])
rotated_stream.merge(method = 1,fill_value='interpolate')

raw_stream.select(channel="BH*").remove_response(inventory=inv,
                                                              output="DISP", plot=False)

raw_stream.select(channel="*H").remove_response(inventory=invp,
                                                              output="DEF", plot=False)


transients_stream.select(channel="BH*").remove_response(inventory=inv,
                                                              output="DISP", plot=False)

transients_stream.select(channel="*H").remove_response(inventory=invp,
                                                              output="DEF", plot=False)






stream_rot_var = raw_stream.copy()

D = tiskit.CleanRotator(stream_rot_var, remove_eq=True)
    
stream_rot_var = D.apply(stream_rot_var)

    
eq_spans = tiskit.TimeSpans.from_eqs(raw_stream.select(channel='*Z')[0].stats.starttime,
                                     raw_stream.select(
                                         channel='*Z')[0].stats.endtime,
                                     minmag=6, days_per_magnitude=0.5)

raw_stream_zero = eq_spans.zero(raw_stream)
transients_stream_zero = eq_spans.zero(transients_stream)
rotated_stream_zero = eq_spans.zero(rotated_stream)
stream_rot_var_zero = eq_spans.zero(stream_rot_var)


fp.psd(raw_stream)
fp.psd(transients_stream)
fp.psd(rotated_stream)

fp.psd_all(rotated_stream_zero,raw_stream_zero,stream_rot_var_zero, transients_stream_zero)


fp.coherogram_spectrogram(raw_stream,nseg=2**12,tw =1,TP=5)

fp.coherogram_spectrogram(rotated_stream,nseg=2**12,tw =1,TP=5)

#%%

import pickle
with open('/Users/mohammadamin/Desktop/Data/YV/'+sta+'/Compliance/rotation_container.pkl', 'rb') as file:
    loaded_rotation_container = pickle.load(file)
    
ang = loaded_rotation_container["angle"]

azi = loaded_rotation_container["azimuth"]

stream_stats = loaded_rotation_container['stream stats']

dates = [str(stream_stats.starttime + ( stream_stats.endtime - stream_stats.starttime ) * 0)[0:10],
         str(stream_stats.starttime + ( stream_stats.endtime - stream_stats.starttime ) * 0.25)[0:10],
         str(stream_stats.starttime + ( stream_stats.endtime - stream_stats.starttime ) * 0.5)[0:10],
         str(stream_stats.starttime + ( stream_stats.endtime - stream_stats.starttime ) * 0.75)[0:10]]

plt.rcParams.update({'font.size': 25})
plt.figure(dpi=300,figsize=[25,15])
plt.subplot(211)
plt.title('Tile Correction Station ' + str(stream_stats.station))
plt.plot(ang,linewidth = 3,color="blue")
plt.ylabel('Angle')
plt.xticks(np.arange(0,len(azi),int(len(azi)/3.5)), dates, rotation=0)  # Set x-axis names with rotation

plt.subplot(212)
plt.plot(azi,linewidth = 3,color='blue')
plt.xticks(np.arange(0,len(azi),int(len(azi)/3.5)), dates, rotation=0)  # Set x-axis names with rotation
plt.ylabel('Azimuth')
plt.xlabel('Time')


at1 = 120
at2 = 150

dates = [str(stream_stats.starttime + ( 3600*24*at1 ) )[0:10],
         str(stream_stats.starttime +  3600*24*(at1 +(at2 - at1 ) * 0.25))[0:10],
         str(stream_stats.starttime +  3600*24*(at1 +(at2 - at1 ) * 0.5))[0:10],
         str(stream_stats.starttime +  3600*24*(at1 +(at2 - at1 ) * 0.75))[0:10]]

plt.rcParams.update({'font.size': 25})
plt.figure(dpi=300,figsize=[25,15])
plt.subplot(211)
plt.title('Tile Correction Station ' + str(stream_stats.station))
plt.plot(ang[4*at1:4*at2],linewidth = 3,color="blue")
plt.ylabel('Angle [°]')
plt.xticks(np.arange(0,len(azi[4*at1:4*at2]),int(len(azi[4*at1:4*at2])/3.5)), dates, rotation=0)  # Set x-axis names with rotation
plt.grid(True)

plt.subplot(212)
plt.plot(azi[4*at1:4*at2],linewidth = 3,color='blue')
plt.xticks(np.arange(0,len(azi[4*at1:4*at2]),int(len(azi[4*at1:4*at2])/3.5)), dates, rotation=0)  # Set x-axis names with rotation
plt.ylabel('Azimuth [°] ')
plt.xlabel('Time')
plt.grid(True)


plt.figure(dpi=300,figsize=(20,12))
plt.subplot(211)
plt.psd(ang[4*at1:4*at2],Fs=1/(6*3600),linewidth=3,color='blue')
plt.yscale('log')
plt.subplot(212)
plt.psd(azi[4*at1:4*at2],Fs=1/(6*3600),linewidth=3,color='blue')
plt.tight_layout()


#%% Compliance Container 
import pickle
import scipy
import compy

# Load compliance Container
print("Loading Compliance Container...")
with open('/Users/mohammadamin/Desktop/Data/YV/'+sta+'/Compliance/com_container.pkl', 'rb') as file:
    loaded = pickle.load(file)

compliance_high = loaded["Compliance"]
f = loaded["Frequency of Compliance"]
uncertainty = loaded['Uncertainty']
coh_high = loaded["Coherence"]
f_coh = loaded["Frequency of Coherence"]
stream = loaded["Stream"]

# compliance_high,coh_high,stream = compy.optimizer(compliance_high, coh_high, stream, f_coh,alpha= 0.85,f_min_com=0.008,f_max_com=0.016)

print("calculating Median of Coherence And Compliance funtion")
Coherence_all = np.median(coh_high,axis=0)
Data = np.median(compliance_high,axis=0)

a1= 5
a2= -7
plt.rcParams.update({'font.size': 25})
plt.figure(dpi=300,figsize=(12,12))
plt.subplot(211)
plt.title(sta)
plt.plot(f,Data)
plt.vlines(x = f[a1], ymin=0, ymax=1e-10,color='black',linestyles="dashed")
plt.vlines(x = f[a2], ymin=0, ymax=1e-10,color='black',linestyles="dashed")
plt.ylabel('Compliance')
plt.grid(True)

plt.subplot(212)
plt.semilogx(f_coh,Coherence_all)
plt.vlines(x = f[a1], ymin=0, ymax=1,color='black',linestyles="dashed")
plt.vlines(x = f[a2], ymin=0, ymax=1,color='black',linestyles="dashed")
plt.hlines(y = 0.8, xmin = f[a1] , xmax = f[a2],linestyles='dashed',linewidth=3,color = 'r')
plt.ylim([0,1])
plt.xlabel('Frequency [Hz]')
plt.ylabel('Coherence')
plt.grid(True)

Data1 = np.median(compliance_high[:,a1:a2],axis=0)
Data = scipy.signal.savgol_filter(Data1,5,2)
f = f[a1:a2]

plt.figure(dpi=300,figsize=(12,12))
plt.subplot(211)
plt.title(sta)
plt.plot(f,Data1,'blue',linewidth=3,label='Compliance')
plt.plot(f,Data,'r',linewidth=3,linestyle='dashed',label='Smoothed Compliance')
plt.ylabel('Compliance')
plt.legend(loc='upper left')
plt.grid(True)

plt.subplot(212)
plt.plot(f_coh,Coherence_all)
plt.hlines(y = 0.8, xmin = f[0] , xmax = f[-1],linestyles='dashed',linewidth=3,color = 'r')
plt.xlim(f[0],f[-1])
plt.xlabel('Frequency [Hz]')
plt.ylabel('Coherence')
plt.grid(True)
plt.tight_layout()
#%%
import obspy
nhnm = obspy.signal.spectral_estimation.get_nhnm()
nlnm = obspy.signal.spectral_estimation.get_nlnm()
import pickle
import scipy
# Load compliance Container
print("Loading Compliance Container...")
with open('/Users/mohammadamin/Desktop/Data/YV/'+sta+'/Compliance/com_container.pkl', 'rb') as file:
    loaded = pickle.load(file)

compliance_high = loaded["Compliance"]
f_com = loaded["Frequency of Compliance"]
uncertainty = loaded['Uncertainty']
coh_high = loaded["Coherence"]
f_coh = loaded["Frequency of Coherence"]
stream = loaded["Stream"]


compliance_high,coh_high,stream = compy.optimizer(compliance_high, coh_high, stream, f_coh,alpha= 0.90,f_min_com=0.008,f_max_com=0.016)

print("calculating Median of Coherence And Compliance funtion")

Coherence_all = np.median(coh_high,axis=0)
Data = np.median(compliance_high,axis=0)
nseg= 2**11
f,Dz = scipy.signal.welch(stream[0].select(component='Z')[0],fs=stream[0][0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
Czp = np.zeros([len(stream),len(Dz)])

Dp = np.zeros([len(stream),len(Dz)])
Dz = np.zeros([len(stream),len(Dz)])

for i in range(0,len(stream)):
    stream[i].detrend(type="simple")
    stream[i].detrend(type="demean")
    stream[i].detrend(type="simple")
    
    f,Dz[i] = scipy.signal.welch(stream[i].select(component='Z')[0],fs=stream[i][0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    f,Dp[i] = scipy.signal.welch(stream[i].select(component='H')[0],fs=stream[i][0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    
    f,Czp[i] = scipy.signal.coherence(stream[i].select(component='Z')[0].data,
                               stream[i].select(component='H')[0].data,
                               fs=stream[i][0].stats.sampling_rate,
                               nperseg =nseg,noverlap=(nseg/2),
                               window=scipy.signal.windows.tukey(nseg,
                               (5*60*stream[i][0].stats.sampling_rate)/nseg))


plt.rcParams.update({'font.size': 40})
plt.figure(dpi=300,figsize=(28,20))
plt.subplot(221)
for i in range(0,len(Dz)):
    
  plt.plot(f,10*np.log10(Dz[i]*(2*np.pi*f)**4),linewidth = 1,color='blue')    
plt.plot(f,10*np.log10(Dz[i]*(2*np.pi*f)**4),label="BHZ",linewidth = 1,color='blue')    
# plt.plot(f,10*np.log10(np.median(Dz,axis=0)*(2*np.pi*f)**4),label="Median",linewidth = 3,color='red')    

    # plt.plot(f,10*np.log10(Dp)-140,label="BDH - 140dB",linewidth = 3)
    
plt.xscale('log')
# plt.title(st[0].stats.network+"."+st[0].stats.station+"  "+str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
plt.xlabel("Frequency (Hz)")
plt.ylabel("Acceleration$((m/s^2)^2/Hz)$ [dB] ")
plt.grid(True)
plt.ylim([-200 ,-80])
# plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
plt.xlim([0.004,0.1])
plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
plt.plot(1/nlnm[0],nlnm[1],'--k')
plt.legend(loc = "upper right",fontsize=30)
plt.tight_layout()

plt.subplot(222)
for i in range(0,len(Dz)):
    
  plt.plot(f,10*np.log10(Dp[i]*(2*np.pi*f)**4),linewidth = 1,color='blue')    
plt.plot(f,10*np.log10(Dp[i]*(2*np.pi*f)**4),label="BDH",linewidth = 1,color='blue')    
# plt.plot(f,10*np.log10(np.median(Dp,axis=0)*(2*np.pi*f)**4),label="Median",linewidth = 3,color='red')    

    # plt.plot(f,10*np.log10(Dp)-140,label="BDH - 140dB",linewidth = 3)
    
plt.xscale('log')
# plt.title(st[0].stats.network+"."+st[0].stats.station+"  "+str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
plt.xlabel("Frequency (Hz)")
plt.ylabel("Pressure $(Pa^2/Hz)$[dB]\n ")
plt.grid(True)
# plt.ylim([-200 ,-80])
# plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
plt.xlim([0.004,0.1])
# plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
# plt.plot(1/nlnm[0],nlnm[1],'--k')
plt.legend(loc = "upper right",fontsize=30)
plt.tight_layout()

plt.subplot(223)
for i in range(0,len(Dz)):
    
  plt.plot(f,Czp[i],linewidth = 0.5,color='blue')    
plt.plot(f,Czp[i],linewidth = 0.5,color='blue',label = "BHZ / BDH")    
plt.plot(f,np.median(Czp,axis=0),linewidth = 3,color='red',label="Median")    
 

plt.vlines(x = 0.0075, ymin=0, ymax=1,color='black',linestyles="dashed",linewidth=3)
plt.vlines(x = 0.016, ymin=0, ymax=1,color='black',linestyles="dashed",linewidth=3,label="Frequency Limit")

    # plt.plot(f,10*np.log10(Dp)-140,label="BDH - 140dB",linewidth = 3)
    
plt.xscale('log')
# plt.title(st[0].stats.network+"."+st[0].stats.station+"  "+str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
plt.xlabel("Frequency (Hz)")
plt.ylabel('Coherence')
plt.grid(True)
# plt.ylim([-200 ,-80])
# plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
plt.xlim([0.004,0.1])
# plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
# plt.plot(1/nlnm[0],nlnm[1],'--k')
plt.legend(loc = "upper right",fontsize=30)
plt.tight_layout()


plt.subplot(224)
for i in range(0,len(Dz)):
    
  plt.plot(f_com,compliance_high[i],linewidth = 0.5,color='blue')    
plt.plot(f_com,compliance_high[i],linewidth = 0.5,color='blue',label = "Compliance")    
plt.plot(f_com,np.median(compliance_high,axis=0),linewidth = 3,color='red',label="Median")    


    # plt.plot(f,10*np.log10(Dp)-140,label="BDH - 140dB",linewidth = 3)
    
plt.yscale('log')
# plt.title(st[0].stats.network+"."+st[0].stats.station+"  "+str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
plt.xlabel("Frequency (Hz)")
plt.ylabel('Compliance')
plt.grid(True)
# plt.ylim([-200 ,-80])
# plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
plt.xlim([0.0066,0.016])
plt.ylim([1e-11,1e-9])
# plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
# plt.plot(1/nlnm[0],nlnm[1],'--k')
plt.legend(loc = "upper right",fontsize=30)
plt.tight_layout()



    











