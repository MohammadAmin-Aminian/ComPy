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

client = Client("RESIF")
# start =UTCDateTime("2012-10-12")
net = "YV"
sta = "RR28"

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


for i in range(0,len(A)):
    stream = stream + read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/'+A[i])
    # stream = stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])
stream.merge(fill_value='interpolate')

#%%
# DPG Calibration
import Pressure_calibration as DPG
gain_factor = DPG.calculate_spectral_ratio(stream,f_min=0.03,f_max=0.05)

#%%
#Compliance Funtion, Rotation's angles 
import compy 

rotated_stream,azimuth,angle = compy.Rotate(stream,time_window = 6)

High_Com_c,High_Czp,High_Com_Stream,f_c,f,uncertainty = compy.Calculate_Compliance(rotated_stream,gain_factor=gain_factor)

High_Com1,High_Czp1,High_Com_Stream1 = compy.optimizer(High_Com_c,High_Czp,High_Com_Stream,f=f,alpha=0.90)


#%%
# Saving Rotated Stream

output_dir = "/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Rotated/'

os.makedirs(output_dir)

compy.split_and_save_stream(rotated_stream,30*24*60,output_dir)

#%% Saving Rotation Container and loading 
import pickle
rotation_container = {
    'stream stats': rotated_stream.select(channel="*Z")[0].stats,
    'angle': angle,
    'azimuth': azimuth,
    'gain factor': gain_factor,
}

# Serialize and save the data container using pickle
file_path = "/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Compliance/'
os.makedirs(file_path)

with open(file_path+'rotation_container.pkl', 'wb') as f1:
    pickle.dump(rotation_container, f1)
    
print("Data container saved.")

# Load the saved data container using pickle
with open('/Users/mohammadamin/Desktop/Data/YV/'+sta+'/Compliance/rotation_container.pkl', 'rb') as file:
    loaded_rotation_container = pickle.load(file)

#%% Saving compliance Container and loading 

com_container = {
    'Compliance': High_Com_c,
    'Frequency of Compliance': f_c,
    'Coherence': High_Czp,
    'Frequency of Coherence': f,
    'Stream': High_Com_Stream
}

# Serialize and save the data container using pickle
file_path = "/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Compliance/'
os.makedirs(file_path)

with open(file_path+'com_container.pkl', 'wb') as f1:
    pickle.dump(com_container, f1)
    
print("Data container saved.")

# Load the saved data container using pickle
with open('/Users/mohammadamin/Desktop/Data/YV/'+sta+'/Compliance/com_container.pkl', 'rb') as file:
    loaded = pickle.load(file)
