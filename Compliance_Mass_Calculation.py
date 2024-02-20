#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 14:45:21 2023

@author: mohammadamin
"""

# import numpy as np
from obspy import read
import os 
import compy

sta = ["RR28","RR29","RR34","RR36","RR38","RR40","RR50","RR52"]
gain_factor = [0.78,0.70,0.58,0.64,0.60,0.74,0.68,0.67]

stream = read()
stream.clear()
import pickle

file_path = "/Users/mohammadamin/Desktop/Compliance/"

for i in range(0,len(sta)):
    # Load compliance Container
    A = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta[i]+'/Rotated/')
    compliance = []
    coherence = []
    stream = []
    
    for ii in range(0, len(A)):
        try:
            rotated_stream = read("/Users/mohammadamin/Desktop/Data/YV/"+sta[i]+'/Rotated/'+A[ii])

            High_Com_c, High_Czp, High_Com_Stream, f_c, f, uncertainty = compy.Calculate_Compliance_beta(rotated_stream, gain_factor=gain_factor[i], time_window=1)
        
        # Add any additional code you want to execute when there's no error

        except Exception as e:
            print(f"An error occurred: {e}")
            # Optionally, add any error handling code here

            continue  # This will skip to the next iteration of the loop

        compliance.append(High_Com_c)
        coherence.append(High_Czp)
        stream.append(High_Com_Stream)
        
        
        
    com_container = {
            'Compliance': compliance,
            'Frequency of Compliance': f_c,
            'Uncertainty':uncertainty,
            'Coherence': coherence,
            'Frequency of Coherence': f,
            'Stream': stream
        }

    # Serialize and save the data container using pickle

    with open(file_path+f"com_container_{sta[i]}.pkl", 'wb') as f1:
        pickle.dump(com_container, f1)
