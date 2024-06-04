#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 11:13:05 2023

@author: mohammadamin
"""

import numpy as np
import matplotlib.pyplot as plt
import tiskitpy as tiskit
from obspy.clients.fdsn import Client
import scipy 
import obspy
nhnm = obspy.signal.spectral_estimation.get_nhnm()
nlnm = obspy.signal.spectral_estimation.get_nlnm()

def start():
    print(" Plot Parametres Loaded")
    plt_params()

def Calculate_Compliance(stream,f_min_com = 0.007,f_max_com = 0.017,gain_factor=0.66,time_window=1):
    nseg = 2**12
    TP = 1
    server="RESIF"
    client = Client(server)

    net = stream[0].stats.network
    sta = stream[0].stats.station
    
    invz = client.get_stations(
        network=net,
        station=sta,
        channel="BHZ",
        location="*",
        level="response")

    print("Splitting The stream into "+ str(time_window)+"-Hour" + ' Windows')
    print("...")
    # split_streams = split_stream(stream, duration = time_window*60*60)
    split_streams = cut_stream_with_overlap(stream, time_window*60*60, overlap = ((time_window*60)-1)*60)
    
    f, Dpp = scipy.signal.welch(split_streams[0].select(component='H')[0], fs=split_streams[0][0].stats.sampling_rate,
                                          nperseg=nseg, noverlap=(nseg*0.5),
                                          window=scipy.signal.windows.tukey(nseg,
                                                                            (TP*60*split_streams[0][0].stats.sampling_rate)/nseg),
                                          average='median')    
    k = wavenumber((2*np.pi*f), -invz[0][0][0].elevation)


    Dp = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    Dz = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    Czp = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    Dzp = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    # Com = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    # Com_Admitance = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    
    print('Calculating Coherence and power spectrum density funtions')
    
    for i in range(0,len(split_streams)):

        f, Dp[i] = scipy.signal.welch(split_streams[i].select(component='H')[0], fs=split_streams[i][0].stats.sampling_rate,
                                          nperseg=nseg, noverlap=(nseg*0.5),
                                          window=scipy.signal.windows.tukey(nseg,
                                                                            (TP*60*split_streams[i][0].stats.sampling_rate)/nseg),
                                          average='median')

        f, Dz[i] = scipy.signal.welch(split_streams[i].select(component='Z')[0], fs=split_streams[i][0].stats.sampling_rate,
                                      nperseg=nseg, noverlap=(nseg*0.5),
                                      window=scipy.signal.windows.tukey(nseg,
                                                                       (TP*60*split_streams[i][0].stats.sampling_rate)/nseg),
                                      average='median')

        f, Czp[i] = np.sqrt(scipy.signal.coherence(split_streams[i].select(component='Z')[0].data,
                                           split_streams[i].select(component='H')[
                                               0].data,
                                           fs=split_streams[i][0].stats.sampling_rate,
                                           nperseg=nseg, noverlap=(nseg*0.5),
                                           window=scipy.signal.windows.tukey(nseg,
                                                                             (TP*60*split_streams[i][0].stats.sampling_rate)/nseg)))


        f , Dzp[i] = scipy.signal.csd(split_streams[i].select(component='Z')[0],split_streams[i].select(component='H')[0],
                               fs=split_streams[i][0].stats.sampling_rate,
                               nperseg=nseg, noverlap=(nseg*0.5),
                               window=scipy.signal.windows.tukey(nseg,
                                                                 (TP*60*split_streams[i][0].stats.sampling_rate)/nseg))
        print( len( split_streams ) - i )
        
    High_Czp = [] 
    High_Dz = []
    High_Dp = [] 
    
    High_Com = []
    High_Com_ae = []
    High_Com_aw = []
    High_Com_gf = []
    High_Com_nothing = []
    
    High_Com_Admitance = []
    High_Com_Stream = []

    Fc1 = np.sqrt(9.8/(2*np.pi*0.5*-invz[0][0][0].elevation))
    Fc2 = np.sqrt(9.8/(2*np.pi*2*-invz[0][0][0].elevation))
    
    coherence_mask = (f >= f_min_com) & (f <= f_max_com)

    
    ppsd_mask = (f >= 0.03) & (f <= 0.08)
    
    mask_dz = (f >= Fc2) & (f <= Fc1)

    # Treshhold to remove exclude those compiance function that are lesser or greater than median of all funtion +- percentage
    # percentage = 0.1
    
    print("Data optipization ...")
    print("Removing gravitatinal Attraction of ocean surface waves...")
    print("Computing Compliance funtion and admittance funtion ")
    for i in range(0,len(Czp)):
        
        # if np.mean(Czp[i][coherence_mask]) > 0.95 and (1-percentage)*np.mean(Com_Admitance[:,coherence_mask]) < np.mean(Com_Admitance[i][coherence_mask] < (1+percentage)*np.mean(Com_Admitance[:,coherence_mask])) :
        if np.median(Czp[i][coherence_mask]) > 0.90 and np.mean(Czp[i][ppsd_mask]) > 0.80 and 0.5<np.mean(Dp[i][ppsd_mask]) < 1  and np.mean(Dz[i][ppsd_mask]) > 10e-17 and np.mean(Dz[i][mask_dz]) < 10e-12 :
        # if np.mean(Czp[i][coherence_mask]) > 0.80:
            
            pa_ratio,aw,hw = gravitational_attraction(np.sqrt((Dp[i])),-invz[0][0][0].elevation,f,pw=1025)
            
            #This term cause by change in distance of the sensor from the earth's center of mass Crawford report
            omega = 2* np.pi * f
            
            ad = (omega**2) /(omega**2 + 3.07* 10e-6)
            # Com = (k * Czp[i]* np.sqrt(np.abs(Dz[i]+aw**2) / np.abs(Dp[i]))) / gain_factor
            
            Com = (k * Czp[i])* ad*(pa_ratio + (np.sqrt(Dz[i])) / (np.sqrt(Dp[i]) * gain_factor))

            Com_ae = ad*((k * Czp[i]* (np.sqrt(Dz[i]))) / (np.sqrt(Dp[i]) * gain_factor))
            
            Com_aw = (k * Czp[i])* (pa_ratio + (np.sqrt(Dz[i])) / (np.sqrt(Dp[i]) * gain_factor))

            Com_gf = (k * Czp[i]* (np.sqrt(Dz[i]))) / (np.sqrt(Dp[i]) * gain_factor)
            
            Com_nothing = (k * Czp[i]* (np.sqrt(Dz[i]))) / (np.sqrt(Dp[i]) )


            Com_Admitance = k * (Dzp[i]/Dp[i]) * gain_factor
            
            
            High_Czp.append(Czp[i])
    
            High_Dz.append(Dz[i])
         
            High_Dp.append(Dp[i]*(gain_factor**2))
            
            High_Com.append(Com)
            High_Com_ae.append(Com_ae)
            High_Com_aw.append(Com_aw)
            High_Com_gf.append(Com_gf)
            High_Com_nothing.append(Com_nothing)
            
            High_Com_Admitance.append(Com_Admitance)
         
            High_Com_Stream.append(split_streams[i])
        
            print(i)
            
    plt.figure(dpi=300,figsize=(14,20))
    plt.title(stream[0].stats.network+"."+stream[0].stats.station)
    plt.subplot(411)                                   
    for i in range(0,len(High_Dz)):
        plt.semilogx(f,10*np.log10(High_Dz[i]*(2*np.pi*f)**4),linewidth = 0.5,color='r')
    plt.semilogx(f,10*np.log10(np.median(High_Dz*(2*np.pi*f)**4,axis=0)),linewidth = 2,color='b',label='Median')
    plt.vlines(x = f_min_com, ymin=-200, ymax=-80,color='black',linestyles="dashed",label="High Coherence Band")
    plt.vlines(x = f_max_com, ymin=-200, ymax=-80,color='black',linestyles="dashed")
    plt.vlines(x = Fc1, ymin=-200, ymax=80,color='green',linestyles="solid",label="Maximum Frequency of IG")
    plt.vlines(x = Fc2 , ymin=-200, ymax=80,color='green',linestyles="solid")
    
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Vertical Acc [m/s^2] dB')
    plt.text(0.01, 0.8, 'a)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
    plt.grid(True)
    plt.legend(loc='upper right',fontsize=17)
    plt.xlim([0.001,1])
    plt.ylim([-200,-80])
  
    plt.subplot(412)                                   
    for i in range(0,len(High_Dp)):
        plt.semilogx(f,10*np.log10(High_Dp[i]),linewidth = 0.5,color='r')
    plt.semilogx(f,10*np.log10(np.median(High_Dp,axis=0)),linewidth = 2,color='b',label='Median')
    plt.vlines(x = f_min_com, ymin=-40, ymax=100,color='black',linestyles="dashed",label="High Coherence Band")
    plt.vlines(x = f_max_com, ymin=-40, ymax=100,color='black',linestyles="dashed")
    plt.vlines(x = Fc1, ymin=-40, ymax=100,color='green',linestyles="solid",label="Maximum Frequency of IG")
    plt.vlines(x = Fc2 , ymin=-40, ymax=100,color='green',linestyles="solid")
    
    
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Pressure')
    plt.text(0.01, 0.8, 'b)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
    plt.grid(True)
    plt.legend(loc='upper right',fontsize=17)
    plt.xlim([0.001,1])
    plt.ylim([-40,80])

    plt.subplot(413)                                   
    for i in range(0,len(High_Czp)):
        plt.semilogx(f,High_Czp[i],linewidth = 0.5,color='r')
    plt.semilogx(f,np.median(High_Czp,axis=0),linewidth = 2,color='b',label='Median')    
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Coherence')
    plt.text(0.01, 0.8, 'c)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
    plt.grid(True)        
    plt.xlim([0.001,1])
    plt.vlines(x = f_min_com, ymin=0, ymax=1,color='black',linestyles="dashed",label="High Coherence Band")
    plt.vlines(x = f_max_com, ymin=0, ymax=1,color='black',linestyles="dashed")
    
    plt.vlines(x = Fc1, ymin=0, ymax=1,color='green',linestyles="solid",label="Maximum Frequency of IG")
    plt.vlines(x = Fc2 , ymin=0, ymax=1,color='green',linestyles="solid")
    
    plt.legend(loc='upper right',fontsize=17)

    plt.subplot(414)                                   
    for i in range(0,len(High_Com)):
        plt.loglog(f,High_Com[i],linewidth = 0.5,color='r')
    plt.loglog(f,np.median(High_Com,axis=0),linewidth = 2,color='b',label='Median')
    # plt.loglog(f,np.mean(High_Com,axis=0),linewidth = 2,color='r',label='Median')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Compliance')
    plt.text(0.01, 0.8, 'd)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
    plt.grid(True)        
    plt.xlim([f_min_com,f_max_com])
    plt.ylim([10e-13,10e-10])
    plt.tight_layout()
    plt.legend(loc='upper right',fontsize=17)

    indices = np.where((f >= 0.005) & (f <= 0.025))
    f_c = f[indices]
    High_Com_c = []
    High_Com_ae1 = []
    High_Com_aw1 = []
    High_Com_gf1 = []
    High_Com_nothing1 = []
    
    print("Filtering ...")
    print("Computing uncertainty of Compliance measurments ...")

    for i in range(0,len(High_Com)):
        High_Com_c.append(High_Com[i][indices])
        High_Com_ae1.append(High_Com_ae[i][indices])
        High_Com_aw1.append(High_Com_aw[i][indices])
        High_Com_gf1.append(High_Com_gf[i][indices])
        High_Com_nothing1.append(High_Com_nothing[i][indices])
            
    uncertainty = f_c.copy()
    High_Com_c = np.array(High_Com_c)
    High_Com_ae1 = np.array(High_Com_ae1)
    High_Com_aw1 = np.array(High_Com_aw1)
    High_Com_gf1 = np.array(High_Com_gf1)
    High_Com_nothing1 = np.array(High_Com_nothing1)
    
    overlap_points = int(nseg * 0.9)

    # Calculate the number of non-overlapping data points per window
    non_overlap_points = nseg - overlap_points

    # Determine the number of windows that can be obtained from the data
    number_of_window = int((split_streams[i][0].stats.npts - overlap_points) / non_overlap_points)

    # for i in range(0,len(f_c)):
    #     uncertainty[i] = np.max(High_Com_c[:,i]) - np.min(High_Com_c[:,i])
    
    plt.figure(dpi=300,figsize=(12,12))
    plt.subplot(211)
    plt.plot(f_c,np.median(High_Com_aw1,axis=0),linewidth=3,color='green',label="aw Correction")
    plt.plot(f_c,np.median(High_Com_ae1,axis=0),linewidth=3,color='blue',label="ae Correction")

    plt.plot(f_c,np.median(High_Com_gf1,axis=0),linewidth=3,color='red',label="Gain Factor Correction",linestyle="dashed")
    plt.plot(f_c,np.median(High_Com_nothing1,axis=0),linewidth=3,color='black',label="Before Correction",linestyle="dashed")

    plt.yscale('log')
    # plt.loglog(f,np.mean(High_Com,axis=0),linewidth = 2,color='r',label='Median')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Compliance')
    plt.grid(True)        
    plt.xlim([f_min_com,f_max_com])
    plt.xlim([f_min_com,0.025])
    # plt.ylim([10e-13,10e-11])
    plt.tight_layout()
    plt.legend(loc='lower right',fontsize=17)

    plt.subplot(212)
    plt.plot(f_c,np.median(High_Com_aw1,axis=0) - np.median(High_Com_ae1,axis=0),linewidth=3,color='green',label="aw effect")
    plt.plot(f_c, np.median(High_Com_gf1,axis=0)-np.median(High_Com_ae1,axis=0),linewidth=3,color='blue',label="ae effect")

    plt.yscale('log')
    # plt.loglog(f,np.mean(High_Com,axis=0),linewidth = 2,color='r',label='Median')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Compliance')
    plt.grid(True)        
    plt.xlim([f_min_com,f_max_com])
    plt.xlim([f_min_com,0.025])
    # plt.ylim([10e-13,10e-11])
    plt.tight_layout()
    plt.legend(loc='upper right',fontsize=17)
    
    uncertainty_theory = Comliance_uncertainty(High_Com_c[0],High_Czp[0][indices],number_of_window)
  
    uncertainty = np.max(High_Com_c,axis=0) - np.min(High_Com_c,axis=0)
    
    return(High_Com_c,High_Czp,High_Com_Stream,f_c,f,uncertainty,uncertainty_theory)
#%%

def sliding_window(a, ws, ss=None, hann=True):
    """
    Function to split a data array into overlapping, possibly tapered sub-windows

    Parameters
    ----------
    a : :class:`~numpy.ndarray`
        1D array of data to split
    ws : int
        Window size in samples
    ss : int
        Step size in samples. If not provided, window and step size
         are equal.

    Returns
    -------
    out : :class:`~numpy.ndarray`
        1D array of windowed data
    nd : int
        Number of windows

    """
    
    if ss is None:
        ss = ws

    if len(a) < ws:
        raise ValueError("The length of the array must be at least as long as the window size.")

    # Calculate the number of windows
    valid = len(a) - ws
    nd = 1 + valid // ss
    out = np.ndarray((nd, ws), dtype=a.dtype)

    for i in range(nd):
        start = i * ss
        stop = start + ws
        # Ensure we don't go past the end of the array
        stop = min(stop, len(a))
        window = a[start: stop]
        if len(window) < ws:
            # For the last window, if it's not full-length, pad it with zeros
            window = np.pad(window, (0, ws - len(window)), 'constant')
        if hann:
            out[i] = window * np.hanning(ws)
        else:
            out[i] = window

    return out, nd

#%%
def Calculate_Compliance_beta(stream,f_min_com = 0.007,f_max_com = 0.02,gain_factor=0.66,time_window = 2):
    nseg = 2**12
    TP = 2
    kk = 1
    times_step = 5 # Minute
    server="RESIF"
    client = Client(server)

    net = stream[0].stats.network
    sta = stream[0].stats.station
    
    invz = client.get_stations(
        network=net,
        station=sta,
        channel="BHZ",
        location="*",
        level="response")

    # ws = int(time_window*60*60 * stream[0].stats.sampling_rate)
    # step_size = 1
    # ss = int(ws*(step_size/60)) 
    split_streams = cut_stream_with_overlap(stream, time_window*60*60, overlap = ((time_window*60)-times_step)*60)
    
    f, Dpp = scipy.signal.welch(split_streams[0].select(component='H')[0], 
                                fs=stream[0].stats.sampling_rate,
                                nperseg=nseg,
                                noverlap=(nseg*0.5),
                                nfft = kk*nseg,
                                window=scipy.signal.windows.tukey(nseg,(TP*60*stream[0].stats.sampling_rate)/nseg),
                                average='median')    
    k = wavenumber((2*np.pi*f), -invz[0][0][0].elevation)

    # Dp = np.zeros([len(split_streams),int(kk*nseg / 2 + 1)])
    # Dz = np.zeros([len(split_streams),int(kk*nseg / 2 + 1)])
    Czp = np.zeros([len(split_streams),int(kk*nseg / 2 + 1)])
    # Dzp = np.zeros([len(split_streams),int(kk*nseg / 2 + 1)])
    # Com = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    # Com_Admitance = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    
    print('Calculating Coherence and power spectrum density funtions')
    
    High_Cohrence_index = []
    for i in range(0,len(split_streams)):
    
        f, Czp[i] = scipy.signal.coherence(split_streams[i].select(component='Z')[0].data,
                                                   split_streams[i].select(component='H')[0].data,
                                                   fs=stream[0].stats.sampling_rate,
                                                   nperseg=nseg,
                                                   noverlap=(nseg*0.5),
                                                   nfft = kk*nseg,
                                                   window=scipy.signal.windows.tukey(nseg,(TP*60*stream[0].stats.sampling_rate)/nseg))

        print( len( split_streams ) - i )
        
    Czp = np.sqrt(Czp)
    
    
    print("Data optipization ...")
    print("Removing gravitatinal Attraction of ocean surface waves...")
    print("Computing Compliance funtion and admittance funtion ")
    
    
    f1 = np.argmin(np.abs(f-f_min_com))
    f2 = np.argmin(np.abs(f-f_max_com))
    
    f11 = np.argmin(np.abs(f-0.02))
    f22 = np.argmin(np.abs(f-0.05))
    
    treshhold_coh = 0.90
    treshhold_coh2 = 0.90
    
    Dp = []
    Dz = []
    Dzp = []
    High_Czp = [] 
    High_Com_Stream = []

    for i in range(0,len(Czp)):
        
        if np.median(Czp[i][f1:f2]) > treshhold_coh and np.mean(Czp[i][f1:f2]) > treshhold_coh2 :

            Dp.append(scipy.signal.welch(split_streams[i].select(component='H')[0], 
                                          fs = stream[0].stats.sampling_rate,
                                          nperseg=nseg, 
                                          noverlap=(nseg*0.5),
                                          nfft = kk*nseg,
                                          window=scipy.signal.windows.tukey(nseg,(TP*60*stream[0].stats.sampling_rate)/nseg),
                                          average='median')[1]*(gain_factor**2))

            Dz.append(scipy.signal.welch(split_streams[i].select(component='Z')[0],
                                          fs=stream[0].stats.sampling_rate,
                                          nperseg=nseg,
                                          noverlap=(nseg*0.5),
                                          nfft = kk*nseg,
                                          window=scipy.signal.windows.tukey(nseg,(TP*60*stream[0].stats.sampling_rate)/nseg),
                                          average='median')[1])
            Dzp.append(scipy.signal.csd(split_streams[i].select(component='Z')[0],
                                            split_streams[i].select(component='H')[0],
                                            fs=stream[0].stats.sampling_rate,
                                            nperseg=nseg, noverlap=(nseg*0.5),
                                            window=scipy.signal.windows.tukey(nseg,
                                                                      (TP*60*stream[0].stats.sampling_rate)/nseg))[1])
            
            High_Czp.append(Czp[i])    
            
            High_Com_Stream.append(split_streams[i])

            
            High_Cohrence_index.append(i)
        print(len(Czp)-i)

    
    Good_Dp = []
    Good_Dz = [] 
    Good_Czp = [] 
    Good_Com = []
    Good_Com_Stream = []


    # Treshhold to remove exclude those compiance function that are lesser or greater than median of all funtion +- percentage
    # percentage = 0.1
    

    # loc = np.where((np.mean(Czp,axis=1)) == max(np.mean(Czp,axis=1)))

    for i in range(0,len(High_Czp)):

        if np.min(High_Czp[i][f1:f2]) > 0.60 and 20 < 10*np.log10(np.median(Dp[i][f1:f2])) < 50 and  10*np.log10(np.median(Dp[i][f11:f22])) < 50 and -180 < np.mean(10*np.log10(Dz[i][f1:f2]*(2*np.pi*f[f1:f2])**4)) < -110 :
            pa_ratio,aw,hw = gravitational_attraction(np.sqrt((Dp[i])),-invz[0][0][0].elevation,f,pw=1025)
            
            #This term cause by change in distance of the sensor from the earth's center of mass Crawford report
            omega = 2* np.pi * f
            
            ad = (omega**2) /(omega**2 + 3.07* 10e-6)
            
            Com = (k * High_Czp[i])* ad*(pa_ratio + (np.sqrt(Dz[i])) / (np.sqrt(Dp[i]) ))


            # Com_Admitance = k * (Dzp[i]/Dp[i]) 
            
            Good_Dp.append(Dp[i])
            Good_Dz.append(Dz[i]) 
            Good_Czp.append(High_Czp[i])
            Good_Com.append(Com)
            Good_Com_Stream.append(High_Com_Stream[i])
         
        
        print(len(High_Czp) - i )
        
    plt_params()
    plt.figure(dpi=300,figsize=(25,25))
    plt.suptitle(stream[0].stats.network+"."+stream[0].stats.station)
    plt.subplot(221)       
                            
    for i in range(0,len(Good_Dz)):
        plt.semilogx(f,10*np.log10(Good_Dz[i]*(2*np.pi*f)**4),linewidth = 0.5,color='g')
    plt.semilogx(f,10*np.log10(np.median(Good_Dz*(2*np.pi*f)**4,axis=0)),linewidth = 2,color='b',label='Median')
    plt.vlines(x = f_min_com, ymin=-200, ymax=-80,color='black',linestyles="dashed",label="High Coherence Band")
    plt.vlines(x = f_max_com, ymin=-200, ymax=-80,color='black',linestyles="dashed")
    plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
    plt.plot(1/nlnm[0],nlnm[1],'--k')
    # plt.vlines(x = Fc1, ymin=-200, ymax=80,color='green',linestyles="solid",label="Maximum Frequency of IG")
    # plt.vlines(x = Fc2 , ymin=-200, ymax=80,color='green',linestyles="solid")
    
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Vertical Acc [m/s^2] dB')
    # plt.text(0.01, 0.8, 'a)', transform=plt.gca().transAxes,fontweight='bold')
    plt.grid(True)
    plt.xlim([0.001,0.1])
    plt.ylim([-200,-80])
    plt.legend(loc='upper right',fontsize=25)

    plt.subplot(222)                                   
    for i in range(0,len(Good_Dp)):
        plt.semilogx(f,10*np.log10(Good_Dp[i]),linewidth = 0.5,color='g')
    plt.semilogx(f,10*np.log10(np.median(Good_Dp,axis=0)),linewidth = 2,color='b',label='Median')
    plt.vlines(x = f_min_com, ymin=-40, ymax=100,color='black',linestyles="dashed",label="High Coherence Band")
    plt.vlines(x = f_max_com, ymin=-40, ymax=100,color='black',linestyles="dashed")
    # plt.vlines(x = Fc1, ymin=-40, ymax=100,color='green',linestyles="solid",label="Maximum Frequency of IG")
    # plt.vlines(x = Fc2 , ymin=-40, ymax=100,color='green',linestyles="solid")
    
    
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Pressure')
    # plt.text(0.01, 0.8, 'b)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
    plt.grid(True)
    plt.xlim([0.001,0.1])
    plt.ylim([-40,80])
    plt.legend(loc='upper right',fontsize=25)

    plt.subplot(223)                                   
    for i in range(0,len(Good_Czp)):
        plt.semilogx(f,Good_Czp[i],linewidth = 0.5,color='g')
    plt.semilogx(f,np.mean(Good_Czp,axis=0),linewidth = 2,color='b',label='Median')    
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Coherence')
    # plt.text(0.01, 0.8, 'c)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
    plt.grid(True)        
    plt.xlim([0.001,0.1])
    plt.vlines(x = f_min_com, ymin=0, ymax=1,color='black',linestyles="dashed",label="High Coherence Band")
    plt.vlines(x = f_max_com, ymin=0, ymax=1,color='black',linestyles="dashed")
    
    # plt.vlines(x = Fc1, ymin=0, ymax=1,color='green',linestyles="solid",label="Maximum Frequency of IG")
    # plt.vlines(x = Fc2 , ymin=0, ymax=1,color='green',linestyles="solid")
    plt.ylim([0,1])
    
    plt.legend(loc='upper right',fontsize=25)

    plt.subplot(224)
    for i in range(0,len(Good_Com)):
        plt.plot(f,Good_Com[i],linewidth = 0.25,color='g')
        
    # plt.plot(f,np.median(High_Com,axis=0),linewidth = 2,color='b',label='Median')
    
    plt.errorbar(f,np.mean(Good_Com,axis=0), yerr=np.std(Good_Com,axis=0), 
                  fmt='o',markersize=10,color='black', ecolor='black', capsize=10,linewidth=5,label="Median")

    # plt.loglog(f,np.mean(High_Com,axis=0),linewidth = 2,color='r',label='Median')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Compliance')
    # plt.text(0.01, 0.8, 'd)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
    plt.grid(True)        
    plt.xlim([f_min_com,f_max_com])
    plt.ylim([-1e-11,1*10e-10])
    plt.tight_layout()
    plt.legend(loc='upper right',fontsize=25)


    indices = np.where((f >= 0.001) & (f <= 0.1))
    f_c = f[indices]
    Good_Com_c = []
    # High_Com_ae1 = []
    # High_Com_aw1 = []
    # High_Com_gf1 = []
    # High_Com_nothing1 = []
    
    print("Filtering ...")
    print("Computing uncertainty of Compliance measurments ...")

    for i in range(0,len(Good_Com)):
        Good_Com_c.append(Good_Com[i][indices])
        # High_Com_ae1.append(High_Com_ae[i][indices])
        # High_Com_aw1.append(High_Com_aw[i][indices])
        # High_Com_gf1.append(High_Com_gf[i][indices])
        # High_Com_nothing1.append(High_Com_nothing[i][indices])
            
    uncertainty = f_c.copy()
    Good_Com_c = np.array(Good_Com_c)
    # High_Com_ae1 = np.array(High_Com_ae1)
    # High_Com_aw1 = np.array(High_Com_aw1)
    # High_Com_gf1 = np.array(High_Com_gf1)
    # High_Com_nothing1 = np.array(High_Com_nothing1)
    
    # overlap_points = int(nseg * 0.5)

    # Calculate the number of non-overlapping data points per window
    # non_overlap_points = nseg - overlap_points

    # Determine the number of windows that can be obtained from the data
    # number_of_window = int((stream[0].stats.npts - overlap_points) / non_overlap_points)

    # for i in range(0,len(f_c)):
    #     uncertainty[i] = np.max(High_Com_c[:,i]) - np.min(High_Com_c[:,i])
    # plt.rcParams.update({'font.size': 25})
    # plt.figure(dpi=300,figsize=(12,12))
    # plt.subplot(211)
    # plt.plot(f_c,np.median(High_Com_aw1,axis=0),linewidth=3,color='green',label="aw Correction")
    # plt.plot(f_c,np.median(High_Com_ae1,axis=0),linewidth=3,color='blue',label="ae Correction")

    # plt.plot(f_c,np.median(High_Com_gf1,axis=0),linewidth=3,color='red',label="Gain Factor Correction",linestyle="dashed")
    # plt.plot(f_c,np.median(High_Com_nothing1,axis=0),linewidth=3,color='black',label="Before Correction",linestyle="dashed")

    # plt.yscale('log')
    # # plt.loglog(f,np.mean(High_Com,axis=0),linewidth = 2,color='r',label='Median')
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Compliance')
    # plt.grid(True)        
    # plt.xlim([f_min_com,f_max_com])
    # plt.xlim([f_min_com,0.025])
    # # plt.ylim([10e-13,10e-11])
    # plt.tight_layout()
    # plt.legend(loc='lower right',fontsize=17)

    # plt.subplot(212)
    # plt.plot(f_c,np.median(High_Com_aw1,axis=0) - np.median(High_Com_ae1,axis=0),linewidth=3,color='green',label="aw effect")
    # plt.plot(f_c, np.median(High_Com_gf1,axis=0)-np.median(High_Com_ae1,axis=0),linewidth=3,color='blue',label="ae effect")

    # plt.yscale('log')
    # # plt.loglog(f,np.mean(High_Com,axis=0),linewidth = 2,color='r',label='Median')
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Compliance')
    # plt.grid(True)        
    # plt.xlim([f_min_com,f_max_com])
    # plt.xlim([f_min_com,0.025])
    # # plt.ylim([10e-13,10e-11])
    # plt.tight_layout()
    # plt.legend(loc='upper right',fontsize=17)
    
    # uncertainty_theory = Comliance_uncertainty(High_Com_c[0],High_Czp[0][indices],number_of_window)
  
    uncertainty = np.std(Good_Com_c,axis=0)
    
    return(Good_Com_c,Good_Czp,Good_Com_Stream,f_c,f,uncertainty)

#%%
def gravitational_attraction(High_Dp,depth_s,f,pw=1025):

    G =  6.6743e-11 # Gravitational constant
    H = depth_s
    kw = wavenumber((2*np.pi*f), depth_s)
    g = 9.8
    psf = High_Dp
    
    hw = (psf * np.cosh(kw * H)) / (pw*g)      
          
    # aw = 2 * np.pi * G * pw * np.exp(-2 * np.pi * kw * H) * hw
    aw = 2 * np.pi * G * pw * np.exp( - kw * H) * hw
    
    pa_ratio = 2 * np.pi * ( G / g ) * np.exp( -kw * H ) * np.cosh( kw * H)
    
    # plt.figure(dpi=300,figsize=(16,12))
    # plt.loglog(f,aw/hw,linewidth=3,color='blue')
    # plt.ylim([10e-9,10e-7])
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Acceleration / sea surface displacement(m/s^2 / m)')
    # plt.tight_layout()
    
    return(pa_ratio,aw,hw)

#%%
def Comliance_uncertainty(compliance_function,coherence_function,number_of_window):
    sigma = compliance_function *((np.sqrt(1-coherence_function**2))/
                                  (np.abs(coherence_function)*np.sqrt(2*number_of_window)))
    return(sigma)
#%%
def Rotate(stream,time_window = 1):

    
    #Split the data into pieces    
    print("Splitting The stream into "+ str(time_window)+"-Hour" + ' Windows')
    print("...")
    split_streams = split_stream(stream, duration = time_window*60*60)

    azimuth = np.zeros([len(split_streams)])
    angle = np.zeros([len(split_streams)])
    variance_percentage = np.zeros([len(split_streams)])

    print("Reducing Tilt Effect")
    for i in range(0,len(split_streams)):
        try:
            print(len(split_streams) - i)
            D = tiskit.CleanRotator(split_streams[i], remove_eq=False)
    
            azimuth[i] = D.azimuth
            angle[i] =  D.angle
            variance_percentage[i] = np.var(split_streams[i].select(channel="*Z")[0].data)

            split_streams[i] = D.apply(split_streams[i])

    # Transfer Function
            F = tiskit.DataCleaner(split_streams[i], ["*1","*2"], n_to_reject=0)
            split_streams[i] = F.clean_stream(split_streams[i])
            
            variance_percentage[i] = variance_percentage[i] / np.var(split_streams[i].select(channel="*Z")[0].data)

            split_streams[i][0].stats.location = '00'
            split_streams[i][1].stats.location = '00'
            split_streams[i][2].stats.location = '00'
            split_streams[i][3].stats.location = '00'
        except Exception as e:
                print(f"Error occurred while processing item: {i}")
    
    time_steps = np.arange(len(angle))
    
    t1 = np.arange(0,len(angle))
    
    time_interval  =  24 # if you do for daily put 1 , if you doing it for hourly put 24 
    # dates = []
    # for i in range(0, int(len(t1) // (time_interval*7))):
    #     print(i)
    #     dates.append(str(stream[0].stats.starttime + (stream[0].stats.endtime - stream[0].stats.starttime) * i / int(len(t1) // (time_interval*30)))[0:10])
        
        
        # Assuming t1 is a list that contains the time values
        # Calculate the positions of the ticks to align with the dates
    # tick_positions = [int(i * len(t1) / (len(dates) - 1)) for i in range(len(dates))]
        
    


   
    for i in range(0,len(angle)):
        if angle[i] < 0 :
            angle[i] = -angle[i]
            azimuth[i] =azimuth[i]+180
            
    azimuth = np.remainder(azimuth, 360.)    
    
    
    #Outlayer Correction
    for i in range(0,len(angle)): 
        if angle[i] > 5 :
            angle[i] = angle[i-1]
            
    for i in range(0,len(azimuth)): 
        if np.mean(azimuth)/1.5> azimuth[i] > np.mean(azimuth)*1.5 :
            azimuth[i] = azimuth[i-1]
        
    variance = np.log10(variance_percentage+10e-0)
    # variance = ((variance - 1) // 10) * 10 + 1
    
    norm = plt.Normalize(variance.min(), variance.max())
    
    cmap = plt.cm.Greys  # Using 'viridis', but you can choose any colormap
    cmap = plt.cm.viridis_r # Using 'viridis', but you can choose any colormap
    # cmap = plt.cm.jet_r # Using 'viridis', but you can choose any colormap
    
    # Plot
    plt.figure(dpi=300,figsize=(30,20))
    plt.subplot(211)
    plt.title("YV." + str(stream[0].stats.station) +'  '+ str(stream[0].stats.starttime)[0:10]+'--'+str(stream[0].stats.endtime)[0:10]+" Tilt ["+str(time_window)+" Hour]")
    
    plt.scatter(time_steps, azimuth, c=variance, cmap=cmap, norm=norm,s=25)
    # plt.plot(time_steps,np.poly1d(np.polyfit(time_steps,azimuth,deg=deg,w=variance))(time_steps),color="red",linewidth=5)
    # for i in range(0,len(eq_spans)):
    #     plt.plot(int((eq_spans.start_times[i] - stream_stats.starttime) // (time_window*3600)),
    #               angle[int(eq_spans.start_times[i] - stream_stats.starttime) // (time_window*3600)],'o', color='white', markersize=5)
    plt.ylim([np.mean(azimuth)-60, np.mean(azimuth)+60])
    plt.grid(True)
    
    plt.colorbar(label='Log10 Variance Reduction')  # Shows the mapping of color to 'varia' values
    plt.xticks([])
    plt.ylabel("Azimuth ["u"\u00b0]")
    plt.text(0.025, 0.95, 'a)', transform=plt.gca().transAxes, fontsize=40, fontweight='bold', va='top')
    
    
    plt.subplot(212)
    plt.scatter(time_steps, angle, c=variance, cmap=cmap, norm=norm,s=25)
    # plt.plot(time_steps,np.poly1d(np.polyfit(time_steps,angle,deg=deg,w=variance))(time_steps),color="red",linewidth=5)
    
    # for i in range(0,len(eq_spans)):
    #     plt.plot(int((eq_spans.start_times[i] - stream_stats.starttime) // (time_window*3600)),
    #               angle[int(eq_spans.start_times[i] - stream_stats.starttime) // (time_window*3600)],'o', color='white', markersize=5)
    plt.grid(True)
    
    plt.colorbar(label='Log10 Variance Reduction')  # Shows the mapping of color to 'varia' values
    plt.xlabel("Time [Date]")
    plt.ylabel("Incident Angle ["u"\u00b0]")
    # plt.xticks(tick_positions, dates, rotation=65)
    plt.ylim([-0.1, 5])
    plt.text(0.025, 0.95, 'b)', transform=plt.gca().transAxes, fontsize=40, fontweight='bold', va='top')
    
    plt.tight_layout()


    
    print("Merging Stream ...")
    rotated_stream = stream.copy()
    rotated_stream.clear()
    for i in range(0,len(split_streams)) : 
        rotated_stream = rotated_stream + split_streams[i]
    rotated_stream.merge(fill_value='interpolate')
    return(rotated_stream,azimuth,angle,variance_percentage)
    # return(azimuth,angle)
#%%

def Rotate_angles(stream,time_window = 1):
    
    #Split the data into pieces    
    print("Splitting The stream into "+ str(time_window)+"-Hour" + ' Windows')
    print("...")
    split_streams = split_stream(stream, duration = time_window*60*60)
    
    for i in range(0,len(split_streams)):
        split_streams[i].detrend('simple')
        
    azimuth = np.zeros([len(split_streams)])
    angle = np.zeros([len(split_streams)])
    
    variance_percentage = np.zeros([len(split_streams)])
    
    print("Reducing Tilt Effect")
    for i in range(0,len(split_streams)):
        try:
            print(len(split_streams) - i)
            D = tiskit.CleanRotator(split_streams[i], remove_eq=False)
    
            azimuth[i] = D.azimuth
            angle[i] =  D.angle
            variance_percentage[i] = D.variance

            # variance_percentage[i] = 100* (1 - np.var(D.apply(split_streams[i]).select(channel="BHZ")) / np.var(split_streams[i].select(channel="BHZ")))

        except Exception as e:
                print(f"Error occurred while processing item: {i}")
    

    return(azimuth,angle,variance_percentage)
    # return(azimuth,angle)
    

#%%
def split_stream(stream, duration):
    split_streams = []
    start_time = stream[0].stats.starttime
    end_time = start_time + duration

    while end_time <= stream[0].stats.endtime:
        split_stream = stream.slice(start_time, end_time)
        split_streams.append(split_stream)

        start_time = end_time
        end_time += duration

    return split_streams
#%%
from obspy import read, Stream

def cut_stream_with_overlap(stream, window_length, overlap):
    """
    Cut an ObsPy stream with a given window length and overlap.

    Parameters:
        stream (obspy.Stream): ObsPy Stream object containing seismic data.
        window_length (float): Length of the window in seconds.
        overlap (float): Overlap between consecutive windows in seconds.

    Returns:
        list of obspy.Stream: List of overlapping segments.
    """
    # Calculate the step size based on window length and overlap
    step_size = window_length - overlap

    # Initialize an empty list to store the overlapping segments
    segments = []

    # Start at the beginning of the stream
    start_time = stream[0].stats.starttime

    # Loop through the stream with the specified step size
    while start_time + window_length <= stream[-1].stats.endtime:
        end_time = start_time + window_length
        segment = stream.slice(start_time, end_time)
        segments.append(segment)
        start_time += step_size

    return segments

#%%
def wavenumber(omega, H):
      """
      Function to approximate wavenumber from dispersion relation

      H is depth below the seafloor, in meters
      omega is a vector of positive angular frequencies

      Stephen G. Mosher, 2020

      """

      import numpy.polynomial as poly

      g = 9.79329
      N = len(omega)

      # Approximations for k when k*H is very large (deep case) or
      # very small (shallow case)
      k_deep = omega**2 / g
      k_shal = omega / np.sqrt(g * H)

      """
      Alternatively, we can use a rational approximation to
      tanh(x) to solve k for any omega. This approximation gives
      a quartic equation, we take the positive real roots as the
      value of k we're interested in. The rational approximation
      being used is always better than the shallow approximation.
      However, it's only better than the deep approximation if
      k*H < 2.96. Therefore, we keep the solutions to k we find,
      using the rational approximation for k*H < 2.96 and use the
      deep water approximation to solve for k otherwise. The
      average error is just under 1% and the maximum error is
      2.5%.
      """

      k = np.zeros(len(omega))

      for i, om in enumerate(omega):

          if i == 0:
              k[i] = 0.
          else:

              a0 = -27 * om**2 / g    # constant terms
              a1 = 0.                 # no linear terms
              a2 = 27 * H - (9 * om**2 * H**2)/g     # quadratic terms
              a3 = 0.                 # no cubic terms
              a4 = H**3           # quartic terms

              p = poly.Polynomial([a0, a1, a2, a3, a4])
              solu = poly.Polynomial.roots(p)
              positive_roots = solu[solu > 0]
              real_positive_root = \
                  positive_roots[positive_roots.imag == 0].real[0]
              k[i] = real_positive_root

      # For k*H >= 2.96, prefer the deep approximation above
      for i, wavenumber in enumerate(k_deep):
          if wavenumber * H > 2.96:
              k[i] = k_deep[i]

      return k
#%%
def rms(arr):
    """
    Calculate the Root Mean Square (RMS) of an array.

    Parameters:
        arr (numpy.ndarray): The input array.

    Returns:
        float: The RMS value of the array.
    """
    squared_values = np.square(arr)
    mean_squared = np.mean(squared_values)
    rms = np.sqrt(mean_squared)
    return rms

#%%
def split_and_save_stream(stream, interval_minutes, output_dir):
    start_time = stream[0].stats.starttime
    end_time = stream[0].stats.endtime
    station_name = stream[0].stats.station
    current_time = start_time
    while current_time <= end_time:
        interval_start = current_time
        interval_end = current_time + interval_minutes * 60  # Convert minutes to seconds

        if interval_end > end_time:
            interval_end = end_time
        
        interval_stream = stream.slice(interval_start, interval_end)
        
        # Save the interval stream to a file
        filename = f"{station_name}_{interval_start.date}.{interval_start.time}-{interval_end.time}.mseed"
        filepath = output_dir + '/' + filename
        interval_stream.write(filepath, format="MSEED")
        
        current_time = interval_end + 1  # Move to the next interval
    
    print("Splitting and saving complete.")

#%%
def optimizer(Com,Czp,stream,f,alpha=0.95,beta = 0.70,zeta = 0.4,f_min_com=0.007,f_max_com=0.018):
    coherence_mask = (f >= f_min_com) & (f <= f_max_com)
    coherence_mask2 = (f >= 0.025) & (f <= 0.03)
    High_Czp = [] 
    High_Com = []
    High_Com_Stream = []
    for i in range(0,len(Czp)):
        
        if np.mean(Czp[i][coherence_mask]) > alpha:
            if np.mean(Czp[i][coherence_mask2]) < beta:
                if np.min(Czp[i][coherence_mask]) > zeta:
            # Com1 = (k * Czp[i]* (np.sqrt(Dz[i]))) / (np.sqrt(Dp[i]) / gain_factor)
                        
                    High_Czp.append(Czp[i])
     
                    High_Com.append(Com[i])

                    High_Com_Stream.append(stream[i])
            
                    print(i)
            
    # plt.rcParams.update({'font.size': 25})
    # plt.figure(dpi=300,figsize=(12,8))
    # for i in range(0,len(Czp)):
    #     plt.semilogx(f,Czp[i],linewidth = 0.5,color='r')
    # for i in range(0,len(High_Czp)):
    #     plt.semilogx(f,High_Czp[i],linewidth = 0.5,color='g')
    # plt.semilogx(f,np.median(High_Czp,axis=0),linewidth = 2,color='b',label='Median of optimized')    
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Coherence')
    # plt.grid(True)        
    # plt.xlim([0.001,1])
    # plt.vlines(x = f_min_com, ymin=0, ymax=1,color='black',linestyles="dashed",label="High Coherence Band")
    # plt.vlines(x = f_max_com, ymin=0, ymax=1,color='black',linestyles="dashed")
    
    
    # plt.legend(loc='upper right',fontsize=17)

    return(High_Com,High_Czp,High_Com_Stream)
#%%
def optimizer_rms(Com,Czp,stream,f,a1,a2,percentage = 10,alpha=0.9,beta=0.5):

    High_Czp = [] 
    High_Com = []
    High_Com_Stream = []
    
    rms_data = rms(np.median(Com,axis=0)[a1:a2])
    for i in range(0,len(Com)):
        
        if rms_data *(1 - (percentage )/100)< rms(Com[i][a1:a2]) < rms_data *(1 + (percentage )/100):
            if 1 > np.median(Czp[i][a1:a2]) > alpha:
                if 1 > np.mean(Czp[i][a1:a2]) > beta:                

                    High_Czp.append(Czp[i])
     
                    High_Com.append(Com[i])

                    High_Com_Stream.append(stream[i])
            
                    print(i)
            
    # plt.rcParams.update({'font.size': 25})
    # plt.figure(dpi=300,figsize=(12,8))
    # for i in range(0,len(Czp)):
    #     plt.semilogx(f,Czp[i],linewidth = 0.5,color='r')
    # for i in range(0,len(High_Czp)):
    #     plt.semilogx(f,High_Czp[i],linewidth = 0.5,color='g')
    # plt.semilogx(f,np.median(High_Czp,axis=0),linewidth = 2,color='b',label='Median of optimized')    
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Coherence')
    # plt.grid(True)        
    # plt.xlim([0.001,1])
    # plt.vlines(x = f_min_com, ymin=0, ymax=1,color='black',linestyles="dashed",label="High Coherence Band")
    # plt.vlines(x = f_max_com, ymin=0, ymax=1,color='black',linestyles="dashed")
    
    
    # plt.legend(loc='upper right',fontsize=17)

    return(High_Com,High_Czp,High_Com_Stream)
#%%
from obspy import Stream

def trim_streams_to_same_length(stream, channels=['BH1', 'BH2', 'BDH', 'BHZ']):
    """
    Trims all channels in the stream to a common time window determined by the maximum start time
    and the minimum end time across all channels.

    :param stream: ObsPy Stream object containing the channels.
    :param channels: List of channel names (e.g., ['BH1', 'BH2', 'BDH', 'BHZ']).
    :return: Stream object with channels trimmed to the same time window.
    """
    # Extract the traces for the specified channels
    traces = [stream.select(channel=ch)[0] for ch in channels]

    # Find the maximum start time and minimum end time across all traces
    start_times = [tr.stats.starttime for tr in traces]
    end_times = [tr.stats.endtime for tr in traces]
    max_start_time = max(start_times)
    min_end_time = min(end_times)

    # Trim all traces to the common time window
    for tr in traces:
        tr.trim(max_start_time, min_end_time)

    # Create a new Stream object with the trimmed traces
    trimmed_stream = Stream(traces)

    return trimmed_stream


#%%
def overlap_checker(stream_splitted):
    start_times = np.zeros(len(stream_splitted))
    end_times = np.zeros(len(stream_splitted))
    
    for i in range(0,len(stream_splitted)):
        start_times[i] = stream_splitted[i][0].stats.starttime
        end_times[i] = stream_splitted[i][0].stats.endtime
    
    start_x_labels = []
    start_x_positions =[]
    start_x_labels.append(str(stream_splitted[0][0].stats.starttime)[0:10])
    start_x_labels.append(str(stream_splitted[20][0].stats.starttime)[0:10])
    start_x_labels.append(str(stream_splitted[-20][0].stats.starttime)[0:10])    
    start_x_labels.append(str(stream_splitted[-1][0].stats.starttime)[0:10])
    
    start_x_positions.append(start_times[0])
    start_x_positions.append(start_times[20])
    start_x_positions.append(start_times[-20])    
    start_x_positions.append(start_times[-1])
    
    
    plt.rcParams.update({'font.size': 25})
    plt.figure(dpi=300,figsize=(16,12))
    # plt.subplot(211)
    plt.hist(start_times,bins=len(stream_splitted))
    plt.xticks(start_x_positions, start_x_labels, rotation=90)  # You can adjust the rotation angle as needed

    plt.xlabel('Time')
    plt.ylabel("Number of Windows That Overlaped")
    plt.legend(loc="upper right")
    plt.grid(True)
    plt.tight_layout()



#%%
def plt_params():
    plt.rcParams['font.size'] = 40
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = ['STIXGeneral']
    plt.rcParams['font.weight'] = 'normal'
    plt.rcParams['mathtext.default'] = 'regular'
    
    plt.rcParams['axes.grid'] = True
    plt.rcParams['axes.grid.which'] = 'major'
    plt.rcParams['axes.labelcolor'] = 'black'
    plt.rcParams['axes.labelpad'] = 4.0
    plt.rcParams['axes.labelsize'] = 'medium'
    plt.rcParams['axes.labelweight'] = 'normal'
    plt.rcParams['axes.linewidth'] = 0.6
    plt.rcParams['axes.titlecolor'] = 'black'
    plt.rcParams['axes.titlepad'] = 5.0
    plt.rcParams['axes.titlesize'] = 'large'
    plt.rcParams['axes.titleweight'] = 'normal'
    plt.rcParams['axes.xmargin'] = 0.05
    plt.rcParams['axes.ymargin'] = 0.05
    
    plt.rcParams['grid.alpha'] = 1
    plt.rcParams['grid.color'] = '#b0b0b0'
    plt.rcParams['grid.linestyle'] = '--'
    plt.rcParams['grid.linewidth'] = 0.6
    
    plt.rcParams['xtick.color'] = 'black'
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['xtick.labelsize'] = 35
    plt.rcParams['xtick.major.pad'] = 2.0
    plt.rcParams['xtick.minor.pad'] = 2.0
    plt.rcParams['xtick.minor.visible'] = True
    plt.rcParams['ytick.color'] = 'black'
    plt.rcParams['ytick.direction'] = 'out'
    plt.rcParams['ytick.labelsize'] = 35
    plt.rcParams['ytick.major.pad'] = 2.0
    plt.rcParams['ytick.minor.pad'] = 2.0
    plt.rcParams['ytick.minor.visible'] = True
    plt.rcParams['xtick.minor.size'] = 10  # Length of minor ticks on x-axis
    plt.rcParams['ytick.minor.size'] = 10  # Length of minor ticks on y-axis
    plt.rcParams['ytick.minor.width'] = 1  # Width of minor ticks on y-axis
    
    # Setting parameters for major ticks
    plt.rcParams['xtick.major.size'] = 10  # Length of major ticks on x-axis
    plt.rcParams['ytick.major.size'] = 10  # Length of major ticks on y-axis
    plt.rcParams['xtick.major.width'] = 2  # Width of major ticks on x-axis
    plt.rcParams['ytick.major.width'] = 2  # Width of major ticks on y-axis
    plt.rcParams['xtick.minor.width'] = 1  # Width of minor ticks on x-axis

    plt.rcParams['legend.borderaxespad'] = 0.0
    plt.rcParams['legend.borderpad'] = 0.5
    plt.rcParams['legend.columnspacing'] = 1.5
    plt.rcParams['legend.edgecolor'] = 'gray'
    plt.rcParams['legend.facecolor'] = 'white'
    plt.rcParams['legend.fancybox'] = False
    plt.rcParams['legend.fontsize'] = 'small'
    plt.rcParams['legend.framealpha'] = 0.8
    plt.rcParams['legend.handleheight'] = 1.0
    plt.rcParams['legend.handlelength'] = 2.0
    plt.rcParams['legend.handletextpad'] = 0.5


        