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

client = Client("RESIF")
# start =UTCDateTime("2012-10-12")
net = "YV"
sta = "RR52"
inv = client.get_stations(
        network=net,
        station=sta,
        channel="B*",
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

i = 7
for i in range(0,len(A)):
    raw_stream = raw_stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])
raw_stream.merge(fill_value='interpolate')

for i in range(0,len(B)):
    transients_stream = transients_stream + read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/'+B[i])
    # stream = stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])
transients_stream.merge(fill_value='interpolate')

for i in range(0,len(C)):
    rotated_stream = rotated_stream + read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Rotated/'+C[i])
    # stream = stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])
rotated_stream.merge(method = 1,fill_value='interpolate')

#RR52
endtime = UTCDateTime("2013-09-15T00")
raw_stream.trim(starttime=rotated_stream[0].stats.starttime,endtime = endtime)
transients_stream.trim(starttime=rotated_stream[0].stats.starttime,endtime = endtime)
rotated_stream.trim(starttime=rotated_stream[0].stats.starttime,endtime = endtime)


#RR38
# starttime = UTCDateTime("2012-12-15T00")
# raw_stream.trim(starttime=rotated_stream[0].stats.starttime,endtime = endtime)
# transients_stream.trim(starttime=rotated_stream[0].stats.starttime,endtime = endtime)
# rotated_stream.trim(starttime=starttime,endtime = rotated_stream[0].stats.endtime)

raw_stream.select(channel="BH*").remove_response(inventory=inv,
                                                              output="DISP", plot=False)

raw_stream.select(channel="*H").remove_response(inventory=invp,
                                                              output="DEF", plot=False)


transients_stream.select(channel="BH*").remove_response(inventory=inv,
                                                              output="DISP", plot=False)

transients_stream.select(channel="*H").remove_response(inventory=invp,
                                                              output="DEF", plot=False)


# fp.coherogram_spectrogram(raw_stream,nseg=2**11,tw =1)

# fp.coherogram_spectrogram(transients_stream,nseg=2**11,tw =1)

# fp.coherogram_spectrogram(rotated_stream,nseg=2**11,tw =1)


# stream_rot_var = raw_stream.copy()

# D = tiskit.CleanRotator(stream_rot_var, remove_eq=True)
    
# stream_rot_var = D.apply(stream_rot_var)


eq_spans = tiskit.TimeSpans.from_eqs(rotated_stream.select(channel='*Z')[0].stats.starttime,
                                     rotated_stream.select(
                                         channel='*Z')[0].stats.endtime,
                                     minmag=6, days_per_magnitude=1.5)


raw_stream_zero = eq_spans.zero(raw_stream)
transients_stream_zero = eq_spans.zero(transients_stream)
rotated_stream_zero = eq_spans.zero(rotated_stream)
# stream_rot_var_zero = eq_spans.zero(stream_rot_var)

fp.psd_h_all(raw_stream,raw_stream_zero,transients_stream_zero,rotated_stream_zero,tw = 1,treshhold_high=1.5e-14, treshhold_low=1e-17)

# fp.psd_h_all_beta(raw_stream,
#                   transients_stream,
#                   transients_stream_zero,
#                   rotated_stream_zero)

fp.coherogram_spectrogram_all(raw_stream)

fp.coherogram_spectrogram(raw_stream,nseg=2**11,tw =1)
fp.coherogram_spectrogram(transients_stream,nseg=2**11,tw =1)
fp.coherogram_spectrogram(transients_stream_zero,nseg=2**11,tw =1)


fp.coherogram_spectrogram_alpha(raw_stream,nseg=2**12,tw = 1)
fp.coherogram_spectrogram_alpha(rotated_stream_zero,nseg=2**12,tw = 1)



fp.psd_h(rotated_stream,nseg=2**12,tw =1)

fp.plot_transfer_function(transients_stream)

#%%
#Glitch plots

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

client = Client("RESIF")
# start =UTCDateTime("2012-10-12")
net = "YV"
sta = "RR52"
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


A.sort()
B.sort()

raw_stream = read()
raw_stream.clear()

transients_stream = read()
transients_stream.clear()

i=  2
N = 24

raw_stream = raw_stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])

transients_stream = transients_stream + read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/'+B[i])
    # stream = stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])


raw_stream.trim(starttime= UTCDateTime(raw_stream[3].stats.starttime),endtime=UTCDateTime(raw_stream[3].stats.starttime+N*3600))
transients_stream.trim(starttime= UTCDateTime(raw_stream[3].stats.starttime),endtime=UTCDateTime(raw_stream[3].stats.starttime+N*3600))


raw_stream.select(channel="BHZ").remove_response(inventory=inv,
                                                              output="DISP", plot=False)


transients_stream.select(channel="BHZ").remove_response(inventory=inv,
                                                              output="DISP", plot=False)


compy.plt_params()
plt.figure(dpi=300,figsize=(20,10))

plt.plot(raw_stream.select(channel="BHZ")[0],linewidth=5,label="Raw")

plt.plot(transients_stream.select(channel="BHZ")[0],linewidth=5,label="Cleaned")

plt.plot(transients_stream.select(channel="BHZ")[0].data - raw_stream.select(channel="BHZ")[0].data ,linewidth=5,label="Cleaned")

plt.title(raw_stream[3].stats.network + "."+raw_stream[3].stats.station)


#%% Rotation Container 
file_path = "/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Compliance/'
file_path_save = "/Users/mohammadamin/Desktop/figures_1Feb"

import pickle
with open(file_path+f'rotation_hourly_{sta}.pkl', 'rb') as file:
    loaded_rotation_container = pickle.load(file)
    

stream_stats =loaded_rotation_container['stream stats']


#RR52
endtime = UTCDateTime("2013-09-15T00")
#RR38
starttime = UTCDateTime("2012-12-15T00")

endtime = stream_stats.endtime
starttime = stream_stats.starttime



eq_spans = tiskit.TimeSpans.from_eqs(starttime,
                                     endtime,
                                     minmag=6, days_per_magnitude=1.5)


aa = (endtime - starttime) // (3600)
# aa = -1 

angle = loaded_rotation_container["angle"][0:int(aa)]

azimuth = loaded_rotation_container["azimuth"][0:int(aa)]

t1 = np.arange(0,len(angle))
    
dates = []
for i in range(0, int(len(t1) // (24*30))):
    print(i)
    dates.append(str(starttime + (endtime - starttime) * i / int(len(t1) // (24*30)))[0:10])
    
    
    # Assuming t1 is a list that contains the time values
    # Calculate the positions of the ticks to align with the dates
tick_positions = [int(i * len(t1) / (len(dates) - 1)) for i in range(len(dates))]
    
# plt.rcParams.update({'font.size': 40})
plt.figure(dpi=300, figsize=(30, 20))

plt.subplot(211)
plt.title("YV." + str(stream_stats.station) +'  '+ str(starttime)[0:10]+'--'+str(endtime)[0:10]+" Tilt [Hourly]")
plt.plot(azimuth+180, '.', color='black')
# for i in range(0,len(eq_spans)):
    # plt.plot(int((eq_spans.start_times[i] - starttime) // (1*3600)),
              # azimuth[int(eq_spans.start_times[i] - starttime) // (1*3600)],'o', color='red', markersize=10)
# plt.plot(int((eq_spans.start_times[i] - starttime) // (1*3600)),
              # azimuth[int(eq_spans.start_times[i] - starttime) // (1*3600)],'o', color='red', markersize=10,label="Event")

plt.xticks([])
plt.ylabel("Azimuth ["u"\u00b0]")
plt.legend(loc="lower right")
plt.grid(True)

plt.subplot(212)
plt.plot(angle, '.', color='black')
    
# for i in range(0,len(eq_spans)):
#     plt.plot(int((eq_spans.start_times[i] - starttime) // (1*3600)),
#               angle[int(eq_spans.start_times[i] - starttime) // (1*3600)],'o', color='red', markersize=10)

# plt.plot(int((eq_spans.start_times[i] - starttime) // (1*3600)),
              # angle[int(eq_spans.start_times[i] - starttime) // (1*3600)],'o', color='red', markersize=10,label="Event")      
plt.xlabel("Time [Date]")
plt.ylabel("Incident Angle ["u"\u00b0]")
plt.grid(True)
plt.legend(loc="lower right")
    
    # Set x-axis names with rotation and align ticks with the dates generated
plt.xticks(tick_positions, dates, rotation=65)
plt.ylim([-5, 5])
plt.tight_layout()


Day_start = 90
Day_end = Day_start + 30
azimuth_c = azimuth[Day_start*24:Day_end*24]
angle_c = angle[Day_start*24:Day_end*24]
ii = 1

padded_azimuth = np.pad(azimuth_c, (0, 1), 'constant', constant_values=0)  # padding with one zero
azimuth_reshaped = azimuth_c.reshape(24*ii, -1)

padded_angle = np.pad(angle_c, (0, 1), 'constant', constant_values=0)  # padding with one zero
angle_reshaped = angle_c.reshape(24*ii, -1)

plt.tight_layout()
plt.savefig(file_path_save + f"YV_{stream_stats.station}_Tilt_Hourly_{starttime}_{endtime}.pdf")




plt.figure(dpi=300, figsize=(30, 20))
plt.subplot(211)
plt.title("YV." + str(stream_stats.station) +" Tilt Daily Stacked [Hourly]")
for i in range(0,int(len(azimuth_c)//24*ii)):
    plt.plot(azimuth_c[i*24*ii:(i+1)*24*ii], color='black',linewidth= 0.5)
plt.plot(np.mean(azimuth_reshaped,axis=1),linewidth=3,color='b')     
plt.xticks([])
plt.ylabel("Azimuth ["u"\u00b0]")
# plt.legend(loc="lower right")
plt.grid(True)

plt.subplot(212)

for i in range(0,int(len(angle_c)//24*ii)):
    plt.plot(angle_c[i*24*ii:(i+1)*24*ii], color='black',linewidth=0.5)

plt.plot(np.mean(angle_reshaped,axis=1),linewidth=5,color="b")
plt.xlabel("Time [Hour]")
plt.ylabel("Incident Angle ["u"\u00b0]")
plt.grid(True)
# plt.legend(loc="lower right")
    
    # Set x-axis names with rotation and align ticks with the dates generated
# plt.xticks(tick_positions, dates, rotation=65)
plt.ylim([-2, 2])
plt.tight_layout()


#Auto Corrolation
Day_start = 0
Day_end = Day_start + 389
azimuth_c = azimuth[Day_start*24:Day_end*24]
angle_c = angle[Day_start*24:Day_end*24]

auto_azimuth = np.correlate(azimuth_c, azimuth_c,mode='full')
auto_angle = np.correlate(angle_c, angle_c,mode='full')
auto_azi_ang = np.correlate(angle_c, azimuth_c,mode='full')


plt.figure(dpi=300, figsize=(20, 20))
plt.suptitle("correlation  "+str(sta),y=0.95)
plt.subplot(311)
plt.title("Azimuth")
plt.plot(auto_azimuth[auto_azimuth.size//2:],linewidth=3,color='b')
plt.grid(True)
plt.xlim([0,48])

plt.subplot(312)
plt.title("Angle")
plt.plot(auto_angle[auto_angle.size//2:],linewidth=3,color='b')
plt.grid(True)
plt.xlim([0,48])

plt.subplot(313)
plt.title("Azimuth-Angle")
plt.plot(auto_azi_ang,linewidth=3,color='b')
plt.grid(True)
plt.tight_layout()

#%%
# Loading Rotation Container 
import pickle
with open(file_path+f'rotation_hourly_{sta}.pkl', 'rb') as file:
    loaded_rotation_container = pickle.load(file)
    
angle = loaded_rotation_container["angle"]
azimuth = loaded_rotation_container["azimuth"]
stream_stats =loaded_rotation_container['stream stats']

for i in range(0,len(angle)):
    if angle[i] < 0 :
        angle[i] = -angle[i]
        azimuth[i] =azimuth[i]+180
        
        
for i in range(0,len(angle)):
    if azimuth[i] < 0 :
        azimuth[i] = -azimuth[i]
                
t1 = np.arange(0,len(angle))
    
dates = []
for i in range(0, int(len(t1) // (24*30))):
    print(i)
    dates.append(str(stream_stats.starttime + (stream_stats.endtime - stream_stats.starttime) * i / int(len(t1) // (24*30)))[0:10])
    
    
    # Assuming t1 is a list that contains the time values
    # Calculate the positions of the ticks to align with the dates
tick_positions = [int(i * len(t1) / (len(dates) - 1)) for i in range(len(dates))]
    
plt.rcParams.update({'font.size': 40})
plt.figure(dpi=300, figsize=(30, 20))

plt.subplot(211)
plt.title("YV." + str(stream_stats.station) +'  '+ str(stream_stats.starttime)[0:10]+'--'+str(stream_stats.endtime)[0:10]+" Tilt [Hourly]")
plt.plot(azimuth, '.', color='black')
# for i in range(0,len(eq_spans)):
#     plt.plot(int((eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)),
#               azimuth[int(eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)],'o', color='red', markersize=10)

plt.xticks([])
plt.ylabel("Azimuth ["u"\u00b0]")
# plt.legend(loc="lower right")
plt.grid(True)

plt.subplot(212)
plt.plot(angle, '.', color='black')
    
# for i in range(0,len(eq_spans)):
#     plt.plot(int((eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)),
#               angle[int(eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)],'o', color='red', markersize=10)

# plt.plot(int((eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)),
                  # angle[int(eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)],'o', color='red', markersize=10,label="Event")       
plt.xlabel("Time [Date]")
plt.ylabel("Incident Angle ["u"\u00b0]")
plt.grid(True)
# plt.legend(loc="lower right", facecolor='lightgreen')
    
    # Set x-axis names with rotation and align ticks with the dates generated
plt.xticks(tick_positions, dates, rotation=65)
plt.ylim([-0.1, 3])
plt.tight_layout()

plt.savefig(file_path_save + f"YV_{stream_stats.station}_Tilt_Hourly_{starttime}_{endtime}.pdf")





#%% Rotation Container Variance
from obspy import UTCDateTime
import tiskitpy as tiskit
import numpy as np
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


file_path = "/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Compliance/'
import pickle
with open(file_path+f'rotation_hourly_Variance_raw_{sta}.pkl', 'rb') as file:
# with open(file_path+f'rotation_daily_Variance_raw_{sta}.pkl', 'rb') as file:
    loaded_rotation_container = pickle.load(file)
    
angle = loaded_rotation_container["angle"]
azimuth = loaded_rotation_container["azimuth"]
stream_stats =loaded_rotation_container['stream stats']
variance = loaded_rotation_container['Variance']

#RR52
endtime = stream_stats.endtime - 45*24*3600


# endtime = stream_stats.endtime 

angle = angle[0:int((endtime -stream_stats.starttime )//3600)]
azimuth = azimuth[0:int((endtime -stream_stats.starttime )//3600)]
variance = variance[0:int((endtime -stream_stats.starttime )//3600)]



eq_spans = tiskit.TimeSpans.from_eqs(stream_stats.starttime,endtime,
                                     minmag=6, days_per_magnitude=1.5)


time_window = 1
time_steps = np.arange(len(angle))

t1 = np.arange(0,len(angle))

time_interval  =  24 # if you do for daily put 1 , if you doing it for hourly put 24 
dates = []
for i in range(0, int(len(t1) // (time_interval*30))):
    print(i)
    dates.append(str(stream_stats.starttime + (endtime - stream_stats.starttime) * i / int(len(t1) // (time_interval*30)))[0:10])
    
    
    # Assuming t1 is a list that contains the time values
    # Calculate the positions of the ticks to align with the dates
tick_positions = [int(i * len(t1) / (len(dates) - 1)) for i in range(len(dates))]
    



   
for i in range(0,len(angle)):
    if angle[i] < 0 :
        angle[i] = -angle[i]
        azimuth[i] =azimuth[i]+180
        
azimuth = np.remainder(azimuth, 360.)    

              
variance = np.log10(variance+10e-0)
# variance = ((variance - 1) // 10) * 10 + 1

norm = plt.Normalize(variance.min(), variance.max())
norm = plt.Normalize(0, variance.max())
cmap = plt.cm.Greys  # Using 'viridis', but you can choose any colormap

# Plot
plt.figure(dpi=300,figsize=(30,20))
plt.subplot(211)
plt.title("YV." + str(stream_stats.station) +'  '+ str(stream_stats.starttime)[0:10]+'--'+str(endtime)[0:10]+" Tilt [Daily]")

plt.scatter(time_steps, azimuth, c=variance, cmap=cmap, norm=norm,s=25)
plt.plot(time_steps,np.poly1d(np.polyfit(time_steps,azimuth,3,w=variance))(time_steps),color="red",linewidth=5)
# for i in range(0,len(eq_spans)):
#     plt.plot(int((eq_spans.start_times[i] - stream_stats.starttime) // (time_window*3600)),
#               angle[int(eq_spans.start_times[i] - stream_stats.starttime) // (time_window*3600)],'o', color='white', markersize=5)
plt.ylim([0, 120])

plt.colorbar(label='Log10 Variance Reduction')  # Shows the mapping of color to 'varia' values
plt.xticks([])
plt.ylabel("Azimuth ["u"\u00b0]")
plt.grid(True)


plt.subplot(212)
plt.scatter(time_steps, angle, c=variance, cmap=cmap, norm=norm,s=25)
plt.plot(time_steps,np.poly1d(np.polyfit(time_steps,angle,3,w=variance))(time_steps),color="red",linewidth=5)

# for i in range(0,len(eq_spans)):
#     plt.plot(int((eq_spans.start_times[i] - stream_stats.starttime) // (time_window*3600)),
#               angle[int(eq_spans.start_times[i] - stream_stats.starttime) // (time_window*3600)],'o', color='white', markersize=5)

plt.colorbar(label='Log10 Variance Reduction')  # Shows the mapping of color to 'varia' values
plt.xlabel("Time [Date]")
plt.ylabel("Incident Angle ["u"\u00b0]")
plt.grid(True)
plt.xticks(tick_positions, dates, rotation=65)
plt.ylim([-0.1, 5])
plt.tight_layout()

plt.savefig(f"YV_{stream_stats.station}_Tilt_Daily.pdf")

# Fs = 1/3600  # In Hertz
# nseg = 2**8
# # Compute the spectrogram
# frequencies, times, S_azimuth = scipy.signal.spectrogram(azimuth, fs=Fs, nperseg=nseg, noverlap=nseg//2, nfft=2**12, window='hann', scaling='spectrum')
# frequencies, times, S_Angle = scipy.signal.spectrogram(angle, fs=Fs, nperseg=nseg, noverlap=nseg//2, nfft=2**12, window='hann', scaling='spectrum')

# t1_spec = np.arange(0,len(times))
# tick_positions_spec = [int(i * len(t1_spec) / (len(dates) - 1)) for i in range(len(dates))]


# # Plotting the Spectrogram
# plt.figure(dpi=300,figsize=(30,20))
# plt.subplot(211)
# plt.title("Spectrogram of Azimuth of  YV."+str(sta))
# plt.pcolormesh(times, frequencies, 10 * np.log10(S_azimuth), shading='gouraud')
# plt.ylabel('Frequency [Hz]')
# plt.yscale('log')
# plt.ylim(0.000001,0.0001)
# plt.xticks([])
# plt.colorbar(label='Intensity [dB]')

# plt.subplot(212)
# plt.title("Spectrogram of Angle")
# plt.pcolormesh(times, frequencies, 10 * np.log10(S_Angle), shading='gouraud')
# plt.ylabel('Frequency [Hz]')
# plt.xlabel('Time [sec]')
# plt.yscale('log')

# plt.ylim(0.000001,0.0001)
# plt.colorbar(label='Intensity [dB]')
# # plt.xticks(tick_positions_spec, dates, rotation=65)

# plt.tight_layout()
# plt.show()


#%%
#Rotation Auto Correltations
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

nseg = 2**10
angle = []
azimuth = []

stations = ["RR28","RR29","RR34","RR36","RR38","RR40","RR50","RR52"]

for i in range(0,len(stations)):
    sta = stations[i]
    file_path = "/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Compliance/'



    with open(file_path+f'rotation_hourly_Variance_raw_{sta}.pkl', 'rb') as file:
        loaded_rotation_container = pickle.load(file)
        
    # with open(file_path+f'rotation_daily_Variance_raw_{sta}.pkl', 'rb') as file:
        # loaded_rotation_container = pickle.load(file)        
    

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

for j in range(0,len(azimuth)):
    for i in range(0,len(angle[j])):
        if angle[j][i] < 0 :
            angle[j][i] = -angle[j][i]
            azimuth[j][i] =azimuth[j][i]+180
        

for j in range(0,len(azimuth)):
    for i in range(0,len(angle)):
        if azimuth[j][i] < 0 :
            azimuth[j][i] = azimuth[j][i]+180
              
        

#Auto Corrolation
auto_angle = []
auto_Azimuth = []
psd_angle = []
psd_Azimuth= []
psd_freq = []
stations = ["RR28","RR29","RR34","RR36","RR38","RR40","RR50","RR52"]
stations1 = ["RR28","RR29","RR34","RR36","RR38","RR40","RR50","RR52"]


jj = 60
# stations1.pop(jj)
for i in range(0,len(angle)):
    if i ==jj:
        pass
    else:
        auto_angle.append(np.correlate(angle[i],angle[i] ,mode='full') - np.mean(np.correlate(angle[i],angle[i] ,mode='full')))
        auto_Azimuth.append(np.correlate(azimuth[i],azimuth[i] ,mode='full') - np.mean(np.correlate(azimuth[i],azimuth[i] ,mode='full')))
        psd_Azimuth.append(scipy.signal.welch(azimuth[i],fs = (1/(3600)),nperseg=nseg,noverlap=nseg//2)[1])
        psd_angle.append(scipy.signal.welch(angle[i],fs = (1/(3600)),nperseg=nseg,noverlap=nseg//2)[1])
        psd_freq.append(scipy.signal.welch(angle[i],fs = (1/(3600)),nperseg=nseg,noverlap=nseg//2)[0])

Hour_start = 0
Hour_End = 96

plt.figure(dpi=300, figsize=(25, 20))
plt.subplot(211)
for i in range(0,len(angle)):
    # plt.title("Correlation Angle " + str(stations[jj]))
    plt.title("Incident Angle [\u00B0]")
    plt.plot(auto_angle[i][auto_angle[i].size//2:],linewidth=3,label=str(stations1[i]),linestyle="dashed")
    plt.grid(True)
    plt.legend(loc="upper right")
    # plt.xlabel([])
    plt.ylabel('Autocorrelation')
    # plt.ylim([-2e6,1e7])
    plt.xlim([Hour_start,Hour_End])
    plt.grid(True)
    
plt.subplot(212)
for i in range(0,len(angle)):
    # plt.title("Correlation Angle " + str(stations[jj]))
    plt.title("Azimuth[\u00B0]")
    # plt.figure(dpi=300, figsize=(25, 20))

    plt.plot(auto_Azimuth[i][auto_Azimuth[i].size//2:],linewidth=3,label=str(stations1[i]),linestyle="dashed")
    plt.grid(True)
    plt.legend(loc="upper right")
    plt.xlabel("Lag [Hour]")
    plt.ylabel('Autocorrelation')
    plt.xlim([Hour_start,Hour_End])
    plt.yscale('log')
    plt.ylim([1e7,2e8])
    plt.grid(True)
    # plt.yscale('log')
plt.tight_layout()


# plt.savefig(file_path_save + f"YV_{stream_stats.station}_AutoCorrelation.pdf")

# Automatic peaking
peaks, _ = scipy.signal.find_peaks(scipy.signal.savgol_filter(np.mean(psd_Azimuth,axis=0),9,1))
top_four_peaks_idx = np.argsort(np.mean(psd_Azimuth,axis=0)[peaks])[-4:]
top_four_frequencies = psd_freq[0][peaks][top_four_peaks_idx]
top_four_psd_values = np.mean(psd_Azimuth,axis=0)[peaks][top_four_peaks_idx]

# Manual peaking
top_four_peaks_idx = [5,42 ,83 ,125 ,167]
x_points = 1/psd_freq[0][top_four_peaks_idx ]
y_points = scipy.signal.savgol_filter(np.mean(psd_Azimuth, axis=0), 9, 1)[top_four_peaks_idx]
labels = ['8d 12h 48m','24h 22m', '12h 20m', '8h 11m', '6h 8m']
# Outputs
print("Top Four Frequencies:", top_four_frequencies)
print("Top Four PSD Values:", top_four_psd_values)



plt.figure(dpi=300, figsize=(25, 30))
plt.suptitle("PSD of Incident Angles and Azimuths")
plt.subplot(211)
for i in range(0,len(angle)):
    # plt.title("Correlation Angle " + str(stations[jj]))
    # plt.loglog(psd_freq[i],psd_angle[i],linewidth=3,label=str(stations1[i]))
    plt.loglog(1/psd_freq[i],scipy.signal.savgol_filter(psd_angle[i],9,1),linewidth=3,label=str(stations1[i]))
    
plt.title("Incident Angle")
plt.loglog(1/psd_freq[i],scipy.signal.savgol_filter(np.mean(psd_angle,axis=0),9,1),linewidth=5,label="Mean",color='black',linestyle="dashed")

plt.grid(True)
plt.legend(loc="upper right")
plt.ylabel('Amplitude [°$^2$/Hz]')
plt.grid(True)
    
plt.subplot(212)
for i in range(0,len(angle)):
    plt.loglog(1/psd_freq[i],scipy.signal.savgol_filter(psd_Azimuth[i],9,1),linewidth=3,label=str(stations1[i]))

plt.title("Azimuth")    
plt.loglog(1/psd_freq[i],scipy.signal.savgol_filter(np.mean(psd_Azimuth,axis=0),9,1),linewidth=5,label="Mean",color='black',linestyle="dashed")
plt.scatter(x_points, y_points, color='red', s=500)  # Example: s=100 for larger circles
y_points[0] = y_points[0]/1.3

for x, y, label in zip(x_points, y_points, labels):
    plt.text(x-2000, y**1.04, label, fontsize=40,rotation=-60)


plt.grid(True)
plt.legend(loc="upper right")
# plt.xlabel("Frequency [Hz]")
plt.xlabel("Period [s]")
plt.ylabel('Amplitude [°$^2$/Hz]')
plt.grid(True)
# plt.xlim(1e5,1e7)
plt.tight_layout()

plt.savefig("Tilt_PSD.pdf")

#%% Compliance Container
import obspy
# nhnm = obspy.signal.spectral_estimation.get_nhnm()
# nlnm = obspy.signal.spectral_estimation.get_nlnm()
import pickle
import scipy
import compy
# Load compliance Container
print("Loading Compliance Container...")
with open('/Users/mohammadamin/Desktop/Data/YV/'+sta+'/Compliance/com_container_old.pkl', 'rb') as file:
    loaded = pickle.load(file)

compliance_high = loaded["Compliance"]
f_com = loaded["Frequency of Compliance"]
uncertainty = loaded['Uncertainty']
coh_high = loaded["Coherence"]
f_coh = loaded["Frequency of Coherence"]
stream = loaded["Stream"]


a1 = 10
a2 = -158
a1 = 5
a2 = -10

f_com = f_com[a1:a2]

# compliance_high,coh_high,stream = compy.optimizer(compliance_high, 
#                                                   coh_high, 
#                                                   stream, 
#                                                   f_coh,
#                                                   alpha=0.70,
#                                                   beta = 0.70,
#                                                   zeta = 0.7,
#                                                   f_min_com=0.007,
#                                                   f_max_com=0.018)

compliance_high,coh_high,stream = compy.optimizer_rms(compliance_high, 
                                                      coh_high, 
                                                      stream,
                                                      f_com,
                                                      a1,
                                                      a2,
                                                      percentage = 50,
                                                      alpha=0.1,
                                                      beta=0.10)


print("calculating Median of Coherence And Compliance funtion")

Coherence_all = np.mean(coh_high,axis=0)
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

# buggy_indices = np.where((f >= 0.048) & (f <= 0.052))[0]
# if buggy_indices.size > 0:
#         # Use interpolation to estimate the correct values for buggy data points
#     good_indices = np.where((f < 0.048) | (f > 0.052))[0]
    

Dz_filtered = [Dz]
Czp_filtered = [Czp]
# for ii in range(0,len(Dz)):
#     correct_values1 = np.interp(buggy_indices, good_indices, Dz[ii][good_indices])
#     correct_values2 = np.interp(buggy_indices, good_indices, Czp[ii][good_indices])
#     Dz_filtered[0][ii][buggy_indices] = correct_values1
#     Czp_filtered[0][ii][buggy_indices] = correct_values2

Dz = np.array(Dz_filtered[0])
Czp = np.array(Czp_filtered[0])


# plt.rcParams.update({'font.size': 40})
plt.figure(dpi=300,figsize=(30,25))
plt.suptitle(stream[0][0].stats.network+"."+stream[0][0].stats.station,y=0.95)    
plt.subplot(221)
for i in range(0,len(Dz)):
    
  plt.plot(f,10*np.log10(Dz[i]*(2*np.pi*f)**4),linewidth = 0.5,color='green')    
plt.plot(f,10*np.log10(Dz[i]*(2*np.pi*f)**4),label="BHZ",linewidth = 2,color='green')    
# plt.plot(f,10*np.log10(scipy.signal.savgol_filter(np.median(Dz,axis=0)*(2*np.pi*f)**4,5,1)),label="Median",linewidth = 3,color='red')    
# plt.plot(f,10*np.log10(np.median(Dz,axis=0)*(2*np.pi*f)**4),label="Median",linewidth = 3,color='red')    

    # plt.plot(f,10*np.log10(Dp)-140,label="BDH - 140dB",linewidth = 3)
plt.vlines(x = f_com[0], ymin=-200, ymax=-80,color='black',linestyles="dashed",linewidth=3)
plt.vlines(x = f_com[-1], ymin=-200, ymax=-80,color='black',linestyles="dashed",linewidth=3,label="Frequency Limit")
    
plt.xscale('log')
# plt.title(st[0].stats.network+"."+st[0].stats.station+"  "+str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
# plt.xlabel("Frequency (Hz)")
plt.ylabel("Acceleration$((m/s^2)^2/Hz)$ [dB] ")
plt.grid(True)
plt.ylim([-200 ,-80])
# plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
plt.xlim([0.005,0.1])
plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
plt.plot(1/nlnm[0],nlnm[1],'--k')
plt.legend(loc = "upper left")
plt.tight_layout()

plt.subplot(222)
for i in range(0,len(Dz)):
    
  plt.plot(f,10*np.log10(Dp[i]),linewidth = 1,color='green')    
plt.plot(f,10*np.log10(Dp[i]),label="BDH",linewidth = 2,color='green')    
# plt.plot(f,10*np.log10(scipy.signal.savgol_filter((np.median(Dp,axis=0)),5,1)),label="Median",linewidth = 5,color='red')
# plt.plot(f,10*np.log10((np.median(Dp,axis=0))),label="Median",linewidth = 5,color='red')


plt.vlines(x = f_com[0], ymin=-40, ymax= 60,color='black',linestyles="dashed",linewidth=3)
plt.vlines(x = f_com[-1], ymin=-40, ymax= 60,color='black',linestyles="dashed",linewidth=3,label="Frequency Limit")

    # plt.plot(f,10*np.log10(Dp)-140,label="BDH - 140dB",linewidth = 3)
    
plt.xscale('log')
# plt.title(st[0].stats.network+"."+st[0].stats.station+"  "+str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
# plt.xlabel("Frequency (Hz)")
plt.ylabel("Pressure $(Pa^2/Hz)$[dB]\n ")
plt.grid(True)
# plt.ylim([-200 ,-80])
# plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
plt.xlim([0.005,0.1])
# plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
# plt.plot(1/nlnm[0],nlnm[1],'--k')
plt.legend(loc = "upper left")
plt.tight_layout()


plt.subplot(223)

for i in range(0,len(Dz)):
  plt.plot(f,Czp[i],linewidth = 0.5,color='green')    
  
plt.plot(f,Czp[i],linewidth = 2,color='green',label = "BHZ / BDH")    
# plt.plot(f,scipy.signal.savgol_filter(np.median(Czp,axis=0),10,5),linewidth = 5,color='red',label="Median")    
# plt.plot(f,np.median(Czp,axis=0),linewidth = 5,color='red',label="Median")    
 
plt.hlines(y = 0.8, xmin = f_com[0] , xmax = f_com[-1],linestyles='dashed',linewidth=3,color = 'r',label='0.8 Treshhold')

plt.vlines(x = f_com[0], ymin=0, ymax=1,color='black',linestyles="dashed",linewidth=3)
plt.vlines(x = f_com[-1], ymin=0, ymax=1,color='black',linestyles="dashed",linewidth=3,label="Frequency Limit")

    # plt.plot(f,10*np.log10(Dp)-140,label="BDH - 140dB",linewidth = 3)
    
plt.xscale('log')
# plt.title(st[0].stats.network+"."+st[0].stats.station+"  "+str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
plt.xlabel("Frequency [Hz]")
plt.ylabel('Coherence')
plt.grid(True)
# plt.ylim([-200 ,-80])
# plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
plt.xlim([0.005,0.1])
# plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
# plt.plot(1/nlnm[0],nlnm[1],'--k')
plt.legend(loc = "upper right")
plt.tight_layout()

plt.subplot(224)
for i in range(0,len(compliance_high)):
  plt.plot(f_com,compliance_high[i][a1:a2],linewidth = 0.5,color='green')    
# plt.plot(f_com,compliance_high[i],linewidth = 0.5,color='blue',label = "Compliance")    
# plt.plot(f_com,np.median(compliance_high,axis=0)[a1:a2],linewidth = 3,color='red',label="Median")    
plt.plot(f_com,compliance_high[i][a1:a2],linewidth = 0.5,color='green',label="Compliance")

plt.errorbar(f_com,np.median(compliance_high,axis=0)[a1:a2], yerr=np.std(compliance_high,axis=0)[a1:a2], 
             fmt='o',markersize=10,color='black', ecolor='black', capsize=10,linewidth=5,label="Median")

# plt.yscale('log')
    # plt.plot(f,10*np.log10(Dp)-140,label="BDH - 140dB",linewidth = 3)
    
# plt.yscale('log')
# plt.title(st[0].stats.network+"."+st[0].stats.station+"  "+str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
plt.xlabel("Frequency [Hz]")
plt.ylabel('Compliance')
plt.grid(True)
# plt.ylim([-200 ,-80])
# plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
plt.xlim([f_com[0],f_com[-1]])
plt.ylim([10e-13,0.2*10e-10])
# plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
# plt.plot(1/nlnm[0],nlnm[1],'--k')

plt.legend(loc = "upper left")
plt.tight_layout()

plt.savefig(file_path_save + f"Compliance_RR52.pdf")


#%% Coherence frequency bands
import pickle
import scipy
import compy
import numpy as np
import matplotlib.pyplot as plt 

sta = "RR52"
# Load compliance Container
print("Loading Compliance Container...")
# with open('/Users/mohammadamin/Desktop/Data/YV/'+sta+'/Compliance/com_container.pkl', 'rb') as file:
with open('/Users/mohammadamin/Desktop/Data/YV/'+sta+'/Compliance/com_container_old.pkl', 'rb') as file:
    loaded = pickle.load(file)
    

    
compliance_high = loaded["Compliance"]
f = loaded["Frequency of Compliance"]
uncertainty = loaded['Uncertainty']
coh_high = loaded["Coherence"]
f_coh = loaded["Frequency of Coherence"]
stream = loaded["Stream"]

Gzz = np.zeros([int(len(stream)/4) ,len(f_coh)])
Gpp = np.zeros([int(len(stream)/4) ,len(f_coh)])
for i in range(0,int(len(stream)/4)):
    F_G , Gzz[i] = scipy.signal.welch(stream[i].select(channel="BHZ"),fs = stream[0][0].stats.sampling_rate,nfft=2*len(f_coh)-1)
    F_G , Gpp[i] = scipy.signal.welch(stream[i].select(channel="BDH"),fs = stream[0][0].stats.sampling_rate,nfft=2*len(f_coh)-1)

Gpp = Gpp*((2*np.pi*F_G)**4)
Gzz = Gzz*((2*np.pi*F_G)**4)
# compliance_high,coh_high,stream = compy.optimizer(compliance_high, coh_high, stream, f_coh,alpha= 0.85,f_min_com=0.008,f_max_com=0.016)

print("calculating Median of Coherence And Compliance funtion")
Coherence_all = np.median(coh_high,axis=0)
Coherence_all_smoothed = scipy.signal.savgol_filter(Coherence_all,5,4)

Data = np.median(compliance_high,axis=0)
 

treshhold = 0.8

a1 = 10
a2 = -158
a1 = 4
a2 = -11

plt.rcParams.update({'font.size': 35})
plt.figure(dpi=300,figsize=(25,15))

# plt.subplot(211)
# plt.semilogx(f_coh,Coherence_all,'b')
plt.title("Coherence (BHZ / BDH) YV." + str(sta))
plt.semilogx(f_coh,Coherence_all_smoothed,'b',linewidth = 3)
plt.vlines(x = f[a1], ymin=0, ymax=1,color='black',linestyles="dashed",linewidth=3)
plt.vlines(x = f[a2-1], ymin=0, ymax=1,color='black',linestyles="dashed",linewidth=3,label="Compliance Frequency Band")
plt.hlines(y =treshhold, xmin = f[a1] , xmax = f[a2-1],linestyles='dashed',linewidth=3,color = 'r',label="0.8 Treshhold")
# plt.plot(F_G,np.mean(Gzz,axis=0)/max(np.mean(Gzz,axis=0)),label = "Normalized PSD (Z'')",linewidth=3)
# plt.plot(F_G,np.mean(Gpp,axis=0)/max(np.mean(Gpp,axis=0)),label = "Normalized PSD (P)",linewidth=3)

plt.ylim([0,1])
plt.xlabel('Frequency [Hz]')
plt.ylabel('Coherence')
plt.axvspan(0.004, 0.021, color='Blue', alpha=0.3, lw=0,label='Infragravity Band (I))')
plt.axvspan(0.021,0.05, color='Orange', alpha=0.3, lw=0,label='low correlation Band (II)')
plt.axvspan(0.05, 0.08, color='Green', alpha=0.3, lw=0,label='Microseismic Band (III) ')
plt.axvspan(0.08, 0.2, color='Pink', alpha=0.6, lw=0,label='Microseismic Band (IV)')

tick_positions = np.logspace(-2.4, 0, num=20)  # Adjust the range and number of ticks as needed
plt.xticks(tick_positions, [f"{val:.3f}" for val in tick_positions], rotation=90)  # Set the tick positions and labels
plt.gca().xaxis.set_minor_locator(plt.NullLocator())  # Disable minor ticks on the x-axis
plt.xlim([0.004,0.2])

plt.grid(True)
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig("/Users/mohammadamin/Desktop/Coherence_bands_RR52.pdf")


# plt.subplot(212)
# plt.title(sta)
# plt.plot(f,Data)
# plt.vlines(x = f[a1], ymin=0, ymax=1e-10,color='black',linestyles="dashed")
# plt.vlines(x = f[a2], ymin=0, ymax=1e-10,color='black',linestyles="dashed")
# plt.ylabel('Compliance')
# plt.grid(True)



# Data1 = np.median(compliance_high[:,a1:a2],axis=0)
# Data = Data1
# Data = scipy.signal.savgol_filter(Data1,10,5)
# f = f[a1:a2]




# plt.figure(dpi=300,figsize=(12,12))
# plt.subplot(211)
# plt.title(sta)
# plt.plot(f,Data1,'blue',linewidth=3,label='Compliance')
# plt.plot(f,Data,'r',linewidth=3,linestyle='dashed',label='Smoothed Compliance')
# plt.ylabel('Compliance')
# plt.legend(loc='upper left')
# plt.grid(True)

# plt.subplot(212)
# plt.plot(f_coh,Coherence_all)
# plt.hlines(y = treshhold, xmin = f[0] , xmax = f[-1],linestyles='dashed',linewidth=3,color = 'r')
# plt.xlim(f[0],f[-1])
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Coherence')
# plt.grid(True)
# plt.tight_layout()
#%%
# Tectonic Difference
import numpy as np
sta = ["RR28","RR29","RR34","RR36","RR38","RR40","RR50","RR52"]
frequency_limits = np.array([[8,-165],[8,-170],[8,-165],[8,-165],[8,-165],[8,-165],[8,-165],[8,-160]])
frequency_limits = np.array([[8,-155],[8,-155],[8,-155],[8,-155],[8,-155],[8,-155],[8,-155],[8,-155]])

percentage = np.array([25,25,25,50,25,25,25,25,])
alpha = np.array([0.85,0.85,0.85,0.85,0.85,0.80,0.85,0.80,])
beta = np.array([0.85,0.845,0.85,0.845,0.84,0.80,0.85,0.80,])

treshhold = 0.8

frequency_limits_2 = [] 
    
import pickle
import scipy
import compy
import numpy as np
import matplotlib.pyplot as plt 
import compy
compy.plt_params()

Data_All = []
s = []

# for RR52 use  old container and a1=4 a2=29

for ii in range(0,len(sta)):
    # Load compliance Container
    print("Loading Compliance Container...")
    with open('/Users/mohammadamin/Desktop/Data/YV/'+sta[ii]+'/Compliance/com_container_old.pkl', 'rb') as file:
        loaded = pickle.load(file)
    
    compliance_high = loaded["Compliance"]
    f = loaded["Frequency of Compliance"]
    uncertainty = loaded['Uncertainty']
    coh_high = loaded["Coherence"]
    f_coh = loaded["Frequency of Coherence"]
    stream = loaded["Stream"]
    a1= frequency_limits[ii][0]
    a2= frequency_limits[ii][1]
    
    # compliance_high,coh_high,stream = compy.optimizer_rms(compliance_high, 
    #                                                       coh_high, 
    #                                                       stream,
    #                                                       f,
    #                                                       a1,
    #                                                       a2,
    #                                                       percentage = percentage[ii],
    #                                                       alpha=alpha[ii],
    #                                                       beta=beta[ii])

    # Gzz = np.zeros([int(len(stream)/4) ,len(f_coh)])
    # Gpp = np.zeros([int(len(stream)/4) ,len(f_coh)])
    # for i in range(0,int(len(stream)/4)):
    #     F_G , Gzz[i] = scipy.signal.welch(stream[i].select(channel="BHZ"),fs = stream[0][0].stats.sampling_rate,nfft=2*len(f_coh)-1)
    #     F_G , Gpp[i] = scipy.signal.welch(stream[i].select(channel="BDH"),fs = stream[0][0].stats.sampling_rate,nfft=2*len(f_coh)-1)
    
    # Gzz = Gzz*((2*np.pi*F_G)**4)
    # compliance_high,coh_high,stream = compy.optimizer(compliance_high, coh_high, stream, f_coh,alpha= 0.85,f_min_com=0.008,f_max_com=0.016)
    
    print("calculating Median of Coherence And Compliance funtion")
    Coherence_all = np.median(coh_high,axis=0)
    Coherence_all_smoothed = scipy.signal.savgol_filter(Coherence_all,5,1)
    
    Data = np.median(compliance_high,axis=0)
    


    print(sta[ii])
    print(f[a1])
    print(f[a2])
    plt.figure(dpi=300,figsize=(30,20))
    
    # plt.subplot(211)
    # plt.semilogx(f_coh,Coherence_all,'b')
    plt.title("Coherence P/Z'' , YV."+str(sta[ii]))
    plt.semilogx(f_coh,Coherence_all_smoothed,'b',linewidth = 7)
    plt.vlines(x = f[a1], ymin=0, ymax=1,color='black',linestyles="dashed",linewidth=5)
    plt.vlines(x = f[a2], ymin=0, ymax=1,color='black',linestyles="dashed",linewidth=5,label="Compliance Frequency Band")
    plt.hlines(y = treshhold, xmin = f[a1] , xmax = f[a2],linestyles='dashed',linewidth=5,color = 'r',label="0.8 Treshhold")
    # plt.plot(F_G,np.mean(Gzz,axis=0)/max(np.mean(Gzz,axis=0)),label = "Normalized PSD (Z'')",linewidth=3)
    # plt.plot(F_G,np.mean(Gpp,axis=0)/max(np.mean(Gpp,axis=0)),label = "Normalized PSD (P)",linewidth=3)
    
    
    plt.ylim([0,1])
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Coherence')
    plt.axvspan(0.004, 0.02, color='Blue', alpha=0.3, lw=0,label='Infragravity Band (I))')
    plt.axvspan(0.02,0.05, color='Orange', alpha=0.3, lw=0,label='low correlation Band (II)')
    plt.axvspan(0.05, 0.08, color='Green', alpha=0.3, lw=0,label='Microseismic Band (III) ')
    plt.axvspan(0.08, 0.2, color='Pink', alpha=0.6, lw=0,label='Microseismic Band (IV)')
    tick_positions = np.logspace(-2.4, 0, num=20)  # Adjust the range and number of ticks as needed
    
    # Set the tick positions and use ScalarFormatter for scientific notation labels
    plt.xticks(tick_positions, [f"{val:.3f}" for val in tick_positions], rotation=90) 
     # Display tick labels in decimal form
    plt.grid(True)
    plt.xlim([0.004,0.2])

    plt.legend(loc='lower right')
    plt.tight_layout()
    
    # plt.subplot(212)
    # plt.title(sta)
    # plt.plot(f,Data)
    # plt.vlines(x = f[a1], ymin=0, ymax=1e-10,color='black',linestyles="dashed")
    # plt.vlines(x = f[a2], ymin=0, ymax=1e-10,color='black',linestyles="dashed")
    # plt.ylabel('Compliance')
    # plt.grid(True)
    
    
    
    Data1 = np.median(compliance_high,axis=0)[a1:a2]
    Data = Data1
    Data = scipy.signal.savgol_filter(Data1,3,1)
    f = f[a1:a2]
    
    frequency_limits_2.append(f)
    
    # plt.figure(dpi=300,figsize=(12,12))
    # plt.subplot(211)
    # plt.title(sta[ii])
    # plt.plot(f,Data1,'blue',linewidth=3,label='Compliance')
    # plt.plot(f,Data,'r',linewidth=3,linestyle='dashed',label='Smoothed Compliance')
    # plt.ylabel('Compliance')
    # plt.legend(loc='upper left')
    # plt.grid(True)
    
    # plt.subplot(212)
    # plt.plot(f_coh,Coherence_all)
    # plt.hlines(y = treshhold, xmin = f[0] , xmax = f[-1],linestyles='dashed',linewidth=3,color = 'r')
    # plt.xlim(f[0],f[-1])
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Coherence')
    # plt.grid(True)
    # plt.tight_layout()
    
    Data_All.append(Data)
    s.append(np.std(compliance_high,axis=0)[a1:a2])
    
    plt.savefig("/Users/mohammadamin/Desktop/Coherence_bands.pdf")



sta = ["RR28","RR29","RR34","RR36","RR38","RR40","RR50","RR52"]

plt.rcParams.update({'font.size': 35})
plt.figure(dpi=300, figsize=[25, 20])
plt.title("Nomalized Compliance")

for jj in range(0, len(Data_All)):
    plt.plot(frequency_limits_2[jj], Data_All[jj], label=str(sta[jj]), linewidth=7)
    # plt.plot(frequency_limits_2[jj], np.mean(Data_All,axis=0) - Data_All[jj], label=str(sta[jj]), linewidth=7)

# Get current handles and labels
handles, labels = plt.gca().get_legend_handles_labels()

# Define the new order of your labels here
# For example, to reverse the order, you can do:
new_order = list((range(len(sta))))
new_order = [5,2,0,6,4,3,1,7]

# Reorder handles and labels
ordered_handles = [handles[i] for i in new_order]
ordered_labels = [labels[i] for i in new_order]

# Create the legend with the new order
plt.legend(ordered_handles, ordered_labels, loc='upper left')

plt.xlabel("Frequency")
plt.grid(True)

# Display the plot
# plt.show()

sta = ["RR28","RR29","RR34","RR36","RR38","RR40","RR50","RR52"]

# plt.rcParams.update({'font.size': 35})
plt.figure(dpi=300, figsize=[30, 16])

for jj in range(0, len(Data_All)):
    plt.subplot(121)
    plt.title("Normalized Compliance")
    plt.legend(ordered_handles, ordered_labels, loc='upper left')
    # plt.plot(frequency_limits_2[jj], np.mean(Data_All,axis=0), label=str("Mean"), linewidth=10,color= 'black',linestyle='dashed')

    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Normalized Compliance")

    plt.grid(True)
    plt.plot(frequency_limits_2[jj], Data_All[jj], label=str(sta[jj]), linewidth=5)
    plt.subplot(122)
    plt.title("Normalized Compliance Deviations from Mean")
    plt.legend(ordered_handles, ordered_labels, loc='upper left')

    plt.xlabel("Frequency [Hz]")

    plt.grid(True)
    plt.plot(frequency_limits_2[jj],- np.mean(Data_All,axis=0) + Data_All[jj], label=str(sta[jj]), linewidth=7)

# Get current handles and labels
handles, labels = plt.gca().get_legend_handles_labels()

# Define the new order of your labels here
# For example, to reverse the order, you can do:
new_order = list((range(len(sta))))
new_order = [5,2,0,6,4,3,1,7]
    
# Reorder handles and labels
ordered_handles = [handles[i] for i in new_order]
ordered_labels = [labels[i] for i in new_order]

# Create the legend with the new order


# Display the plot
# plt.show()

plt.tight_layout()
plt.savefig(file_path_save + f"Compliance_Functions.pdf")

# plt.subplot(212)
# plt.title("Deviation of Station Compliance from Overall Mean")
# for jj in range(0,len(Data_All)):
#     plt.plot(frequency_limits_2[jj],np.mean(Data_All,axis=0) -Data_All[jj],label=str(sta[jj]),linewidth=3)
# plt.legend(loc = 'lower left')
# plt.grid(True)


#%%
# Inversion
import pickle
import inv_compy
import compy

sta1 = ["RR28","RR29","RR34","RR36","RR38","RR40","RR50","RR52"]

print("Loading Compliance Container...")
Inversion_container = [] 
for ii in range(0,len(sta1)):
    with open('/Users/mohammadamin/Desktop/Data/YV/'+sta1[ii]+'/Compliance/Inversion_container.pkl', 'rb') as file:
        Inversion_container.append(pickle.load(file))

sta = Inversion_container[jj]["Station"]
Data = Inversion_container[jj]["compliance Measured"]
f = Inversion_container[jj]["compliance Frequency"]
s = Inversion_container[jj]["uncertainty"]
ncompl = Inversion_container[jj]["compliance Forward"]
starting_model = Inversion_container[jj]["Models"]
vs = Inversion_container[jj]["Shear Velocity"]
vs0 = Inversion_container[jj]["Shear Velocity Starting"]
mis_fit = Inversion_container[jj]["Misfit Fucntion"]
likeli_hood = Inversion_container[jj]["Liklihood"]
sigma_v = Inversion_container[jj]["sigma_v"]
sigma_h = Inversion_container[jj]["sigma_h"]
alpha = Inversion_container[jj]["alpha"]
alpha_pressure = Inversion_container[jj]["alpha_pressure"]
iteration = Inversion_container[jj]["iteration"]
n_layer = Inversion_container[jj]["n_layer"]
burnin = Inversion_container[jj]["burnin"]
mis_fit_trsh = Inversion_container[jj]["mis_fit_trsh"]
accept_rate = Inversion_container[jj]["accept_rate"]


compy.plt_params()
inv_compy.plot_inversion_density(vs,vs0,mis_fit,Data,s,freq=f,sta=sta,burnin=burnin,ncompl=(ncompl),iteration=iteration,mis_fit_trsh = mis_fit_trsh)

inv_compy.plot_inversion_density_all(Inversion_container)


Vs_final = []
for ii in range(0,len(vs)):
    if mis_fit[0][ii] < 9:
        if vs[ii][9000][0] > vs[ii][5000][0]:
            Vs_final.append(vs[ii][0:10000])

vs_median = np.median(Vs_final,axis=0)
vs_dif = np.sum(Vs_final - vs_median,axis=0)

argmean = np.argmin(np.abs(vs_dif),axis=0)[0]
argmin = np.argmin(vs_dif,axis=0)[0]
argmax = np.argmax(vs_dif,axis=0)[0]

















































