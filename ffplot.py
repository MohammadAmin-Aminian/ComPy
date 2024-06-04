#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Mohammad Amin Aminian
"""
# import obstools as obs
import numpy as np
import warnings
import scipy
import scipy.stats
# import math as M
import matplotlib.pyplot as plt
import obspy
import obspy.signal
nhnm = obspy.signal.spectral_estimation.get_nhnm()
nlnm = obspy.signal.spectral_estimation.get_nlnm()
import matplotlib as mpl

warnings.simplefilter("ignore", np.ComplexWarning)
def start():
    print(" Plot Parametres Loaded")
    plt_params()

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
        # no step size was provided. Return non-overlapping windows
        ss = ws

    # Calculate the number of windows to return, ignoring leftover samples, and
    # allocate memory to contain the samples
    valid = len(a) - ss
    nd = (valid) // ss
    out = np.ndarray((nd, ws), dtype=a.dtype)

    if nd == 0:
        if hann:
            out = a * np.hanning(ws)
        else:
            out = a

    for i in range(nd):
        # "slide" the window along the samples
        start = i * ss
        stop = start + ws
        if hann:
            out[i] = a[start: stop] * np.hanning(ws)
        else:
            out[i] = a[start: stop]

    return out, nd
#%%

def psd(st,nseg=2**12):
        
    f,Dz = scipy.signal.welch(st.select(component='Z')[0],
                              fs=st[0].stats.sampling_rate,
                              nperseg = nseg, 
                              noverlap=(nseg/2),
                              nfft=nseg,
                              scaling='density',
                              average = 'median')
    
    f,Dh1 = scipy.signal.welch(st.select(component='1')[0],
                                fs=st[0].stats.sampling_rate,
                                nperseg = nseg, 
                                noverlap=(nseg/2),
                                nfft=nseg,
                                scaling='density',
                                average = 'median')
    
    f,Dh2 = scipy.signal.welch(st.select(component='2')[0],
                                fs=st[0].stats.sampling_rate,
                                nperseg = nseg, 
                                noverlap=(nseg/2),
                                nfft=nseg,
                                scaling='density',
                                average = 'median')
    
    f,Dp = scipy.signal.welch(st.select(component='H')[0],
                              fs=st[0].stats.sampling_rate,
                              nperseg = nseg, 
                              noverlap=(nseg/2),
                              nfft=nseg,
                              scaling='density',
                              average = 'median')

    plt.figure(dpi=300,figsize=(25,20))
    plt.plot(f,10*np.log10(Dz*(2*np.pi*f)**4),label="BHZ",linewidth = 5)             
    plt.plot(f,10*np.log10(Dh1*(2*np.pi*f)**4),label="BH1",linewidth = 5)             
    plt.plot(f,10*np.log10(Dh2*(2*np.pi*f)**4),label="BH2",linewidth = 5)             
    plt.plot(f,10*np.log10(Dp)-140,label="BDH - 140dB",linewidth = 5)
    
    plt.xscale('log')
    plt.title(st[0].stats.network+"."+st[0].stats.station+"  "+str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Pressure $(Pa^2/Hz)$[dB]\n or Acceleration$((m/s^2)^2/Hz)$ [dB] ")
    plt.grid(True)
    plt.ylim([-200 ,-80])
    # plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.xlim([0.005,1])

    # plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
    plt.plot(1/nhnm[0],nhnm[1],'--k')
    plt.plot(1/nlnm[0],nlnm[1],'--k')
    plt.legend(loc = "lower right")
    plt.tight_layout()

#%%

def psd_all(st,st1,st2,st3,nseg=2**13):
    f,Dz = scipy.signal.welch(st.select(component='Z')[0],fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    f,Dz1 = scipy.signal.welch(st1.select(component='Z')[0],fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    f,Dz2 = scipy.signal.welch(st2.select(component='Z')[0],fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    f,Dz3 = scipy.signal.welch(st3.select(component='Z')[0],fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))

    plt.rcParams.update({'font.size': 40})
    plt.figure(dpi=300,figsize=(20,16))
    plt.plot(f,10*np.log10(Dz*(2*np.pi*f)**4),label="Z (Raw - EQ)",linewidth = 3,color='red')             
    plt.plot(f,10*np.log10(Dz1*(2*np.pi*f)**4),label="Z' (Z - glitch)",linewidth = 3,color='blue')             
    plt.plot(f,10*np.log10(Dz2*(2*np.pi*f)**4),label="Z'' ( Rotated Z'')",linewidth = 3,color='orange')             
    plt.plot(f,10*np.log10(Dz3*(2*np.pi*f)**4),label="Z''' (Z''- H1 - H2)",linewidth = 3,color='green')             

    plt.xscale('log')
    plt.title(st[0].stats.network+"."+st[0].stats.station+".BHZ  "+str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Acceleration$((m/s^2)^2/Hz)$ [dB] ")
    plt.grid(True)
    # plt.ylim([0 ,80])
    plt.ylim([-200 ,-80])
    plt.xlim([0.005,1])

    plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
    plt.plot(1/nlnm[0],nlnm[1],'--k')
    plt.legend(loc = "upper left",fontsize=30)
    plt.tight_layout()
 
#%%
def psd_h(st,tw=6,nseg=2**11):
    
    ws = int(tw*60*60 * st[0].stats.sampling_rate)
    Z ,nd = sliding_window(st.select(component='Z')[0].data, ws = ws,hann=True)
    f,Dz = scipy.signal.welch(Z,fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    
    P ,nd = sliding_window(st.select(component='H')[0].data, ws = ws,hann=False)
    f,Dp = scipy.signal.welch(P,fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    
    
    plt.rcParams.update({'font.size': 40})
    plt.figure(dpi=300,figsize=(16,20))
    plt.subplot(211)
    for i in range(0,len(Z)):
        plt.plot(f,10*np.log10(Dz[i]*(2*np.pi*f)**4),linewidth = 0.5,color="r")             
    
    plt.xscale('log')
    plt.title(st[0].stats.network+"."+st[0].stats.station+"  " 
              +str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10] +" every " + str(tw) + " hours")    
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Acc [dB] ")
    plt.grid(True)
    plt.ylim([-200 ,-80])
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
    plt.plot(1/nlnm[0],nlnm[1],'--k')
    plt.legend(loc = "upper right")
    
    plt.subplot(212)
    for i in range(0,len(Z)):
        plt.plot(f,10*np.log10(Dp[i]),linewidth = 0.5,color="r")             
    
    plt.xscale('log')   
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Pressure [dB] ")
    plt.grid(True)
    # plt.ylim([-200 ,-80])
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.legend(loc = "upper right")
    plt.tight_layout()
    
    
    
    plt.rcParams.update({'font.size': 40})   
    plt.figure(dpi=300,figsize=(25,25))
    plt.subplot(121)
    # plt.subplot(211)
    plt.imshow(10*np.log10(Dz*(2*np.pi*f)**4),aspect = 'auto',extent=[f[0],f[len(f)-1],0,nd*tw])
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    plt.xscale('log')
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.ylabel("Time (Hour)")
    plt.xlabel('Frequency (Hz)')
    plt.title("Acceleration Spectrogram")
    plt.tight_layout()
    color_bar = plt.colorbar(format='%.0f dB')  # You can adjust the format as needed
    color_bar.set_label("ACC (dB)")
    
    plt.subplot(122)
    plt.imshow(10*np.log10(Dp),aspect = 'auto',extent=[f[0],f[len(f)-1],0,nd*tw])
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    plt.xscale('log')
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.xlabel('Frequency (Hz)')
    plt.title("Pressure Spectrogram")
    plt.tight_layout()
    color_bar = plt.colorbar(format='%.0f dB')  # You can adjust the format as needed
    color_bar.set_label("Pressure (dB)")


#%%

def psd_h_all(st, st1, st2, st3, tw=6, nseg=2**11, treshhold_high=1e-14, treshhold_low=1e-17):
    
    # Initialize lists for filtered data
    # Dz_filtered = []
    Dz1_filtered = []
    Dz2_filtered = []
    Dz3_filtered = []


    # Window size calculation
    ws = int(tw * 60 * 60 * st[0].stats.sampling_rate)
    
    # Compute for st
    Z, _ = sliding_window(st.select(component='Z')[0].data, ws=ws, hann=True)
    f, Dz = scipy.signal.welch(Z, fs=st[0].stats.sampling_rate, nperseg=nseg, noverlap=(nseg/2))
    Dz = Dz * (2 * np.pi * f)**4
    Dz1 = Dz.copy()
    
    Clip_indices = np.where((f >= 0.001) & (f <= 0.01))[0]

    # Dz_filtered = [dz for dz in Dz if treshhold_low < np.median(dz) < treshhold_high]
    
    # Compute for st1
    Z1, _ = sliding_window(st1.select(component='Z')[0].data, ws=ws, hann=True)
    _, Dz1 = scipy.signal.welch(Z1, fs=st1[0].stats.sampling_rate, nperseg=nseg, noverlap=(nseg/2))
    Dz1 = Dz1 * (2 * np.pi * f)**4
    
    for i in range(0,len(Dz1)):
        if treshhold_low < np.median(Dz1[i][Clip_indices]) < treshhold_high:
            Dz1_filtered.append(Dz1[i])
            
    # Compute for st2
    Z2, _ = sliding_window(st2.select(component='Z')[0].data, ws=ws, hann=True)
    _, Dz2 = scipy.signal.welch(Z2, fs=st2[0].stats.sampling_rate, nperseg=nseg, noverlap=(nseg/2))
    Dz2 = Dz2 * (2 * np.pi * f)**4
    
    for i in range(0,len(Dz2)):
        if treshhold_low < np.median(Dz2[i][Clip_indices]) < treshhold_high:
            Dz2_filtered.append(Dz2[i])    
    # Compute for st3
    Z3, _ = sliding_window(st3.select(component='Z')[0].data, ws=ws, hann=True)
    _, Dz3 = scipy.signal.welch(Z3, fs=st3[0].stats.sampling_rate, nperseg=nseg, noverlap=(nseg/2))
    Dz3 = Dz3 * (2 * np.pi * f)**4
    
    for i in range(0,len(Dz3)):
        if treshhold_low < np.median(Dz3[i][Clip_indices]) < treshhold_high:
            Dz3_filtered.append(Dz3[i])
            
    # PPSD Bug Cleaner
    buggy_indices = np.where((f >= 0.048) & (f <= 0.052))[0]
    if buggy_indices.size > 0:
        # Use interpolation to estimate the correct values for buggy data points
        good_indices = np.where((f < 0.048) | (f > 0.052))[0]
    
    
        # for ii in range(0,len(Dz_filtered)):
            # correct_values = np.interp(buggy_indices, good_indices, Dz_filtered[ii][good_indices])
            # Dz_filtered[ii][buggy_indices] = correct_values
            
        # for ii in range(0,len(Dz1_filtered)):
        #     correct_values1 = np.interp(buggy_indices, good_indices, Dz1_filtered[ii][good_indices])
        #     Dz1_filtered[ii][buggy_indices] = correct_values1
            
        for ii in range(0,len(Dz2_filtered)):
            correct_values2 = np.interp(buggy_indices, good_indices, Dz2_filtered[ii][good_indices])
            Dz2_filtered[ii][buggy_indices] = correct_values2
            
        for ii in range(0,len(Dz3_filtered)):            
            correct_values3 = np.interp(buggy_indices, good_indices, Dz3_filtered[ii][good_indices])
            Dz3_filtered[ii][buggy_indices] = correct_values3
            
    
    plt.figure(dpi=300,figsize=(25,30))
    
    plt.suptitle(st[0].stats.network+"."+st[0].stats.station+"  " 
              +str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10] +" [Hourly]")
    # plt.suptitle(st[0].stats.network+"."+st[0].stats.station+"  " 
    #           +str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10] +" ["+str(tw)+"-Hour]")   
    plt.subplot(321)
    for i in range(0,len(Z)):
        plt.plot(f,10*np.log10(Dz[i]),linewidth = 0.5,color="gray")   
          
    plt.plot(f,10*np.log10(Dz[i-100]),linewidth = 2,color="gray",label="BHZ [Raw]")             
    plt.plot(f,10*np.log10(np.median(Dz,axis=0)),linewidth = 5.5,color="red",label="Median")             

    plt.xscale('log') 
    # plt.xlabel("Frequency (Hz)")
    plt.ylabel("Acc [dB] ")
    plt.grid(True)
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
    plt.plot(1/nlnm[0],nlnm[1],'--k')
    plt.legend(loc = "upper right")
    

    plt.ylim([-200 ,-80])
    plt.xlim([0.005 ,0.1])
    plt.text(0.05, 0.95, 'a)', transform=plt.gca().transAxes, fontsize=40, fontweight='bold', va='top')

    plt.subplot(322)
    
    for i in range(0,len(Dz1_filtered)):
        plt.plot(f,10*np.log10(Dz1_filtered[i]),linewidth = 0.5,color="gray")             
    
    plt.plot(f,10*np.log10(Dz1_filtered[i]),linewidth =2,color="gray",label="BHZ'[Rotatation - Events]")             
    plt.plot(f,10*np.log10(np.median(Dz1_filtered,axis=0)),linewidth = 5.5,color="Orange",label="Median")             

    plt.xscale('log')   
    # plt.xlabel("Frequency (Hz)")
    # plt.ylabel("Acc [dB] ")
    plt.grid(True)
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.plot(1/nhnm[0],nhnm[1],'--k')
    plt.plot(1/nlnm[0],nlnm[1],'--k')
    plt.legend(loc = "upper right")
    plt.ylim([-200 ,-80])
    plt.xlim([0.005 ,0.1])
    plt.text(0.05, 0.95, 'b)', transform=plt.gca().transAxes, fontsize=40, fontweight='bold', va='top')

    
    plt.subplot(323)
    for i in range(0,len(Dz2_filtered)):
        plt.plot(f,10*np.log10(Dz2_filtered[i]),linewidth = 0.5,color="gray")             

    plt.plot(f,10*np.log10(Dz2_filtered[i]),linewidth = 2,color="gray",label="BHZ''[Deglitch]")             
  
    plt.plot(f,10*np.log10(np.median(Dz2_filtered,axis=0)),linewidth = 5.5,color="blue",label="Median")             

    plt.xscale('log')   
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Acc [dB] ")
    plt.grid(True)
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.plot(1/nhnm[0],nhnm[1],'--k')
    plt.plot(1/nlnm[0],nlnm[1],'--k')
    plt.legend(loc = "upper right")
    plt.ylim([-200 ,-80])
    plt.xlim([0.005 ,0.1])
    plt.text(0.05, 0.95, 'c)', transform=plt.gca().transAxes, fontsize=40, fontweight='bold', va='top')

    plt.subplot(324)
    for i in range(0,len(Dz3_filtered)):
        plt.plot(f,10*np.log10(Dz3_filtered[i]),linewidth = 0.5,color="gray")             
    
    plt.plot(f,10*np.log10(Dz3_filtered[i]),linewidth = 2,color="gray",label="BHZ'''[Coherent Noise Removal]")             
    plt.plot(f,10*np.log10(np.median(Dz3_filtered,axis=0)),linewidth = 5.5,color="green",label="Median")             

    plt.xscale('log')  
    plt.xlabel("Frequency [Hz]")
    # plt.ylabel("Acc [dB] ")
    plt.grid(True)
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.plot(1/nhnm[0],nhnm[1],'--k')
    plt.plot(1/nlnm[0],nlnm[1],'--k')
    plt.legend(loc = "upper right")
    plt.ylim([-200 ,-80])
    plt.xlim([0.005 ,0.1])
    plt.text(0.05, 0.95, 'd)', transform=plt.gca().transAxes, fontsize=40, fontweight='bold', va='top')

    plt.subplot(3,2,(5,6))
    
    # plt.figure(dpi=300,figsize=(35,40))
    plt.plot(f,10*np.log10(np.median(Dz,axis=0)),linewidth = 5,label = "BHZ [Raw]" ,color="red")             
    # plt.plot(f,10*np.log10(np.median(Dz_filtered,axis=0)),linewidth = 3,label = "BHZ-High Energy")         
    plt.plot(f,10*np.log10(np.median(Dz1_filtered,axis=0)),linewidth = 5,label = "BHZ'[Rotatation - Events]" ,color="orange")     
    plt.plot(f,10*np.log10(np.median(Dz2_filtered,axis=0)),linewidth = 5,label = "BHZ''[Deglitch]",color="blue")
    plt.plot(f,10*np.log10(np.median(Dz3_filtered,axis=0)),linewidth = 5,label = "BHZ'''[Coherent Noise Removal]",color="green")

    plt.xscale('log')   
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Acc [dB] ")
    plt.grid(True)
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.plot(1/nhnm[0],nhnm[1],'--k')
    # plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
    plt.plot(1/nlnm[0],nlnm[1],'--k')
    plt.legend(loc = "upper right")
    plt.ylim([-170 ,-130])
    plt.xlim([0.005 ,0.1])
    plt.text(0.025, 0.95, 'e)', transform=plt.gca().transAxes, fontsize=40, fontweight='bold', va='top')

    # plt.ylim([-200 ,-80])
    # plt.xlim([0.005 ,0.1])

    plt.legend(loc='upper right')
    plt.tight_layout()
    # plt.savefig("PSD.pdf",dpi = 300)
    plt.savefig("/Users/mohammadamin/Desktop/preprocessing.pdf")

    plt.show()
    
#%%

def psd_h_all_beta(st, st1, st2, st3, nseg=2**12):
    
    # Initialize lists for filtered data

    # Window size calculation
    
    # Compute for st
    Z =st.select(component='Z')[0].data
    f, Dz = scipy.signal.welch(Z, fs=st[0].stats.sampling_rate, nperseg=nseg, noverlap=(nseg/2))
    Dz = Dz * (2 * np.pi * f)**4
    Dz1 = Dz.copy()
    
    # Dz_filtered = [dz for dz in Dz if treshhold_low < np.median(dz) < treshhold_high]
    
    # Compute for st1
    Z1 =st1.select(component='Z')[0].data
    _, Dz1 = scipy.signal.welch(Z1, fs=st1[0].stats.sampling_rate, nperseg=nseg, noverlap=(nseg/2))
    Dz1 = Dz1 * (2 * np.pi * f)**4
    

            
    # Compute for st2
    Z2 =st2.select(component='Z')[0].data
    _, Dz2 = scipy.signal.welch(Z2, fs=st2[0].stats.sampling_rate, nperseg=nseg, noverlap=(nseg/2))
    Dz2 = Dz2 * (2 * np.pi * f)**4
  
    # Compute for st3
    Z3 =st3.select(component='Z')[0].data
    _, Dz3 = scipy.signal.welch(Z3, fs=st3[0].stats.sampling_rate, nperseg=nseg, noverlap=(nseg/2))
    Dz3 = Dz3 * (2 * np.pi * f)**4


    plt.figure(dpi=300,figsize=(20,20))
    
    # plt.text(0.0052, -137, 'e)', fontsize=50,weight="bold" ,va='top', ha='left')
    # plt.figure(dpi=300,figsize=(35,40))
    plt.plot(f,10*np.log10(Dz),linewidth = 5,label = "BHZ [Raw]" ,color="red")             
    # plt.plot(f,10*np.log10(np.median(Dz_filtered,axis=0)),linewidth = 3,label = "BHZ-High Energy")         
    plt.plot(f,10*np.log10(Dz1),linewidth = 5,label = "BHZ'" ,color="Orange")     
    plt.plot(f,10*np.log10(Dz2),linewidth = 5,label = "BHZ''",color="blue")
    plt.plot(f,10*np.log10(Dz3),linewidth = 5,label = "BHZ'''",color="green")

    plt.xscale('log')   
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Acc [dB] ")
    plt.grid(True)
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    # plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
    # plt.plot(1/nlnm[0],nlnm[1],'--k')
    plt.legend(loc = "upper right",facecolor='lightgreen')
    plt.tick_params(axis='both', which='minor', length=10, width=1)
    plt.tick_params(axis='both', which='major', length=10, width=2)
    plt.ylim([-200 ,-80])
    plt.xlim([0.001 ,1])
    
    # plt.ylim([-200 ,-80])
    # plt.xlim([0.005 ,0.1])

    plt.legend(loc='upper right',facecolor='lightgreen')
    plt.tight_layout()
    # plt.savefig("PSD.pdf",dpi = 300)
    plt.show()
#%%
import numpy as np
import scipy.signal

def plot_transfer_function(st, nseg=2**12, TP=5):
    '''
    Plot the transfer function between channels.
    
    Parameters
    ----------
    st : Stream
        Seismic data stream.
    nseg : int, optional
        Number of segments for FFT. The default is 2**12.
    TP : int, optional
        Time for taper of 2 sides of each segment, uses Tukey window. 
        The default is 5 minutes.
    '''
    fs = st[0].stats.sampling_rate
    window = scipy.signal.windows.tukey(nseg, (TP*60*fs)/nseg)
    
    # Compute CSD and PSD
    f, Pxx = scipy.signal.welch(st.select(component='Z')[0].data, fs=fs, nperseg=nseg, noverlap=(nseg*0.9), window=window)
    f, Pxy = scipy.signal.csd(st.select(component='Z')[0].data, st.select(component='H')[0].data, fs=fs, nperseg=nseg, noverlap=(nseg*0.9), window=window)
    
    # Compute the transfer function H(f) = Pxy(f) / Pxx(f)
    H = Pxy / Pxx
    
    # Plot the magnitude and phase of the transfer function
    plt.figure(figsize=(15, 6))
    
    # Magnitude plot
    plt.subplot(1, 2, 1)
    plt.semilogx(f, 20 * np.log10(np.abs(H)), label='Magnitude')
    plt.title('Transfer Function Magnitude')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude (dB)')
    plt.grid(True)
    plt.legend()
    
    # Phase plot
    plt.subplot(1, 2, 2)
    plt.semilogx(f, np.angle(H, deg=True), label='Phase')
    plt.title('Transfer Function Phase')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Phase (degrees)')
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    plt.show()


#%%
def coh(st,nseg=2**12,TP=5):
    '''
     Coherence between channels
     overlap = 0.5
     
    Parameters
    ----------
    st : Stream
    nseg : number of segment to calculate fft The default is 2**12.
    TP : time for taper of 2 side of eatch segments, uses tukey window. 
    The default is 5 minutes.

    Returns
    -------
    None.

    '''
    f,Czp = scipy.signal.coherence(st.select(component='Z')[0].data,
                               st.select(component='H')[0].data,
                               fs=st[0].stats.sampling_rate,
                               nperseg =nseg,noverlap=(nseg*0.5),
                               window=scipy.signal.windows.tukey(nseg,
                               (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,C12 = scipy.signal.coherence(st.select(component='1')[0].data,
                                st.select(component='2')[0].data,
                                fs=st[0].stats.sampling_rate,
                                nperseg =nseg,noverlap=(nseg*0.5),
                                window=scipy.signal.windows.tukey(nseg,
                                (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,Cz1 = scipy.signal.coherence(st.select(component='Z')[0].data,
                                st.select(component='1')[0].data,
                                fs=st[0].stats.sampling_rate,
                                nperseg =nseg,noverlap=(nseg*0.5),
                                window=scipy.signal.windows.tukey(nseg,
                                (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,Cz2 = scipy.signal.coherence(st.select(component='Z')[0].data,
                                st.select(component='2')[0].data,
                                fs=st[0].stats.sampling_rate,
                                nperseg =nseg,noverlap=(nseg*0.5),
                                window=scipy.signal.windows.tukey(nseg,
                                (TP*60*st[0].stats.sampling_rate)/nseg))
    
    plt.rcParams.update({'font.size': 25})
    plt.figure(dpi=300,figsize=((15,15)))

    plt.subplot(221)
    plt.semilogx(f,Czp,color = 'b',linewidth = 3 , label = 'P & Z')
    plt.title(""+st[0].stats.network+"."+
              st[0].stats.station+"  "+
              str(st[0].stats.starttime)[0:10]+"--"+
              str(st[0].stats.endtime)[0:10])    
    plt.ylabel("Coherence ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([(2/nseg),1])
    plt.legend(loc="upper right")

    plt.subplot(222)
    plt.semilogx(f,C12,color = 'b',linewidth = 3 , label = 'H1 & H2')
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([(2/nseg),1])
    plt.legend(loc="upper right")

    plt.subplot(223)
    plt.semilogx(f,Cz1,color = 'b',linewidth = 3 , label = 'H1 & Z')
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Coherence ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([(2/nseg),1])
    plt.legend(loc="upper right")

    
    plt.subplot(224)
    plt.semilogx(f,Cz2,color = 'b',linewidth = 3 , label = 'H2 & Z')
    plt.ylabel("Coherence ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([(2/nseg),1])
    plt.legend(loc="upper right")

    plt.tight_layout()
#%%

def coh_compliance(st,nseg=2**12,TP=5):
    '''
     Coherence between channels
     overlap = 0.9
     
    Parameters
    ----------
    st : Stream
    nseg : number of segment to calculate fft The default is 2**12.
    TP : time for taper of 2 side of eatch segments, uses tukey window. 
    The default is 5 minutes.

    Returns
    -------
    None.

    '''
    f,Czp = scipy.signal.coherence(st.select(component='Z')[0].data,
                               st.select(component='H')[0].data,
                               fs=st[0].stats.sampling_rate,
                               nperseg =nseg,noverlap=(nseg*0.9),
                               window=scipy.signal.windows.tukey(nseg,
                               (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,C1p = scipy.signal.coherence(st.select(component='1')[0].data,
                               st.select(component='H')[0].data,
                               fs=st[0].stats.sampling_rate,
                               nperseg =nseg,noverlap=(nseg*0.9),
                               window=scipy.signal.windows.tukey(nseg,
                               (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,C2p = scipy.signal.coherence(st.select(component='2')[0].data,
                               st.select(component='H')[0].data,
                               fs=st[0].stats.sampling_rate,
                               nperseg =nseg,noverlap=(nseg*0.9),
                               window=scipy.signal.windows.tukey(nseg,
                               (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,C12 = scipy.signal.coherence(st.select(component='1')[0].data,
                               st.select(component='2')[0].data,
                               fs=st[0].stats.sampling_rate,
                               nperseg =nseg,noverlap=(nseg*0.9),
                               window=scipy.signal.windows.tukey(nseg,
                               (TP*60*st[0].stats.sampling_rate)/nseg))
    
    plt.figure(dpi=300,figsize=((20,20)))
    plt.suptitle(""+st[0].stats.network+"."+
              st[0].stats.station+"  "+
              str(st[0].stats.starttime)[0:10]+"--"+
              str(st[0].stats.endtime)[0:10],y=0.95)  
    plt.subplot(221)
    plt.semilogx(f,Czp,color = 'b',linewidth = 3 , label = 'Z & P')
  
    plt.ylabel("Coherence ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([(2/nseg),1])
    plt.legend(loc="upper right")

    plt.subplot(222)
    plt.semilogx(f,C12,color = 'b',linewidth = 3 , label = 'H1 & H2')
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([(2/nseg),1])
    plt.legend(loc="upper right")

    plt.subplot(223)
    plt.semilogx(f,C1p,color = 'b',linewidth = 3 , label = 'H1 & P')
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Coherence ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([(2/nseg),1])
    plt.legend(loc="upper right")

    
    plt.subplot(224)
    plt.semilogx(f,C2p,color = 'b',linewidth = 3 , label = 'H2 & P')
    plt.ylabel("Coherence ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([(2/nseg),1])
    plt.legend(loc="upper right")

    plt.tight_layout()
    
#%%  
def coh_h(st,tw=1,nseg=2**12,lim=[0.005,0.1],TP = 5):
    '''
     Coherence Over time
     Coherogram
     overlap = 0.5
     
    Parameters
    ----------
    st : Stream
    tw : Time in hour for sliding window. The default is 6 hours.
    nseg : number of segment to calculate fft The default is 2**12.
    TP : time for taper of 2 side of eatch segments, uses tukey window. 
    The default is 5 minutes.

    Returns
    -------
    None.

    '''
    jj = 1
    ws = int(tw*60*60 * st[0].stats.sampling_rate)
    Z ,nd = sliding_window(st.select(component='Z')[0].data, ws = ws,hann=True)
    P ,nd = sliding_window(st.select(component='H')[0].data, ws = ws,hann=True)

    f,Czp = scipy.signal.coherence(Z,
                                   P,
                                   fs=st[0].stats.sampling_rate,
                                   nperseg =nseg,
                                   noverlap=(nseg*0.5),
                                   nfft=nseg*jj,
                                   window=scipy.signal.windows.tukey(nseg,
                                   (TP*60*st[0].stats.sampling_rate)/nseg))
    number_month = (st[0].stats.endtime - st[0].stats.starttime) // (30*24*3600/7)
    
    dates = []
    for i in range(0, int(number_month)+1):
        # print(i)
        dates.append(str(st[0].stats.starttime + (st[0].stats.endtime - st[0].stats.starttime) * i / number_month)[0:10])
        
    
    # tick_positions = [int(i * len(Czp[0]) / (len(dates) - 1)) for i in range(len(dates))]

    f1 = np.argmin(np.abs(f-0.007))
    f2 = np.argmin(np.abs(f-0.015))
    treshhold_coh = 0.9
    Czp_High = []
    # loc = np.where((np.mean(Czp,axis=1)) == max(np.mean(Czp,axis=1)))
    plt.rcParams.update({'font.size': 40})   
    plt.figure(dpi=300,figsize=(25,20))
    for i in range(0,len(Czp)):
        if  np.mean(Czp[i][f1:f2]) < treshhold_coh:
            plt.semilogx(f,Czp[i],linewidth = 0.5,color="r")
            
    for i in range(0,len(Czp)):
        if np.mean(Czp[i][f1:f2]) > treshhold_coh:
            plt.semilogx(f,Czp[i],linewidth = 1 ,color="g")
            Czp_High.append(Czp[i])
    # plt.semilogx(f,Czp[loc[0][0]],linewidth = 0.75,color="b")
    plt.semilogx(f,np.mean(Czp,axis=0),linewidth = 5,color="b",label = 'Mean')
    plt.semilogx(f,np.mean(Czp_High,axis=0),linewidth = 5,color="black",label="Mean-High")

    plt.title("Coherence (P&Z)"+"  for every " + str(tw) + " hours  "+st[0].stats.network+"."+st[0].stats.station+"  "+
              str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Coherence ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim(lim[0],lim[1])
    plt.tick_params(axis='both', which='minor', length=10, width=1)
    plt.tick_params(axis='both', which='major', length=10, width=2)
    
    plt.legend(loc="lower right",facecolor="lightgreen")
    
    plt.tight_layout()
    

    plt.rcParams.update({'font.size': 40})   
    plt.figure(dpi=300,figsize=(20,20))
    plt.imshow(Czp,aspect = 'auto',extent=[f[0],f[len(f)-1],0,nd*tw])
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    # plt.yticks(tick_positions, dates, rotation=0)

    plt.tick_params(axis='both', which='minor', length=10, width=1)
    plt.tick_params(axis='both', which='major', length=10, width=2)
    
    plt.xscale('log')
    plt.xlim(lim[0],lim[1])
    plt.ylabel("Time (Hour)")
    plt.xlabel('Frequency (Hz)')
    plt.title("Coherence Between P and Z  As Function of Time")
    color_bar = plt.colorbar(format='%.1f ')  # You can adjust the format as needed
    color_bar.set_label("Coherency")
    color_bar.ax.tick_params(labelsize=25)
    plt.tight_layout()

#%%
def spectrogram(st,nseg=2**12):

    f,t,Sh1 = scipy.signal.spectrogram(st.select(component='1')[0].data,fs=st[1].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2))
    f,t,Sh2 = scipy.signal.spectrogram(st.select(component='2')[0].data,fs=st[2].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2))
    f,t,Sz = scipy.signal.spectrogram(st.select(component='Z')[0].data,fs=st[3].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2))
    f,t,Sp = scipy.signal.spectrogram(st.select(component='H')[0].data,fs=st[0].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2))
    
    plt.set_cmap('jet')
    plt.figure(dpi = 300,figsize=(8,8))
    
    plt.subplot(4,1,1)
    # plt.pcolormesh(t, f, Sh1, shading='gouraud')
    plt.pcolormesh(t, f, 10*np.log10(Sh1))
    plt.ylabel('Frequency [Hz]')
    # plt.xlabel('Time [sec]')
    plt.yscale('log')
    plt.ylim([(st[0].stats.sampling_rate*2/nseg),1])
    plt.title("BH1")
    plt.colorbar()

    plt.subplot(4,1,2)
    # plt.pcolormesh(t, f, Sh2, shading='gouraud')
    plt.pcolormesh(t, f, 10*np.log10(Sh2))
    plt.ylabel('Frequency [Hz]')
    # plt.xlabel('Time [sec]')
    plt.yscale('log')
    plt.ylim([(st[0].stats.sampling_rate*2/nseg),1])
    plt.title("BH2")
    plt.colorbar()

    plt.subplot(4,1,3)
    # plt.pcolormesh(t, f, Sz, shading='gouraud')
    plt.pcolormesh(t, f, 10*np.log10(Sz))
    plt.ylabel('Frequency [Hz]')
    # plt.xlabel('Time [sec]')
    plt.yscale('log')
    plt.ylim([(st[0].stats.sampling_rate*2/nseg),1])
    plt.title("BHZ")
    plt.colorbar()

    
    plt.subplot(4,1,4)
    # plt.pcolormesh(t, f, Sz, shading='gouraud')
    plt.pcolormesh(t, f, 10*np.log10(Sp))
    plt.ylabel('Frequency [Hz]')
    # plt.xlabel('Time [sec]')
    plt.yscale('log')
    plt.ylim([(st[0].stats.sampling_rate*2/nseg),1])
    plt.title("BHZ")
    plt.colorbar()
    plt.tight_layout()

    ###############################################################################
    
    #for the second and forth plot i put manually, i wanted to chose 10e-2. it can be automatically.
    
    plt.figure(dpi = 300,figsize=(12,16))
    plt.subplot(132)
    plt.pcolormesh(f, t, 10*np.log10(Sz.T))
    
    plt.ylabel('Frequency [Hz]')
    plt.xscale('log')
    plt.xlim([0.005,1])
    plt.title("BHZ")
    # plt.colorbar()
    
    plt.subplot(133)
    plt.pcolormesh(f, t, 10*np.log10(Sp.T))
    plt.ylabel('Frequency [Hz]')
    plt.xscale('log')
    plt.xlim([0.005,1])
    plt.title("BDH")
    # plt.colorbar()
    plt.tight_layout()
    
    
    plt.subplot(4,1,2)
    plt.plot(t,Sz[20])
    
    plt.xlim([0,t[len(t)-1]])
    plt.ylabel('Frequency [Hz]')
    plt.yscale('log')
    # plt.ylim([1e-18,1e-14])
    plt.grid(True)


    
    plt.subplot(4,1,4)
    plt.plot(t,Sp[10]) 
    plt.xlim([0,t[len(t)-1]])
    plt.ylabel('Frequency [Hz]')
    plt.xlabel('Time [sec]')
    plt.yscale('log')
    # plt.ylim([10,1e5])
    plt.grid(True)
    plt.tight_layout()

#%%
def coherogram_spectrogram(st,nseg=2**12,tw =1):
    f_min = 0.007
    f_max = 0.018
    Tresh_coh = 0.8
    Tresh_Dz = -170
    Tresh_Dp = 25
    
    TP=5
    f,t,Sz = scipy.signal.spectrogram(st.select(component='Z')[0].data,fs=st[3].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2),window='hann')
    f,t,Sp = scipy.signal.spectrogram(st.select(component='H')[0].data,fs=st[0].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2),window='hann')
    

    ws = int(tw*60*60 * st[0].stats.sampling_rate)
    Z ,nd = sliding_window(st.select(component='Z')[0].data, ws = ws,hann=True)
    P ,nd = sliding_window(st.select(component='H')[0].data, ws = ws,hann=True)

    f,Czp = scipy.signal.coherence(Z,P,fs=st[0].stats.sampling_rate,nperseg =nseg,
                                   noverlap=(nseg*0.5),
                                   window=scipy.signal.windows.tukey(nseg,
                                   (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,Dzz = scipy.signal.welch(Z,fs=st[0].stats.sampling_rate,nperseg =nseg,
                                   noverlap=(nseg*0.5),
                                   window=scipy.signal.windows.tukey(nseg,
                                   (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,Dpp = scipy.signal.welch(P,fs=st[0].stats.sampling_rate,nperseg =nseg,
                                   noverlap=(nseg*0.5),
                                   window=scipy.signal.windows.tukey(nseg,
                                   (TP*60*st[0].stats.sampling_rate)/nseg))
    
    number_month = (st[0].stats.endtime - st[0].stats.starttime) // (30*24*3600)
    
    
    f1 = np.argmin(np.abs(f-f_min))
    f2 = np.argmin(np.abs(f-f_max))
    
    
    dates = []
    for i in range(0, int(number_month)+1):
        # print(i)
        dates.append(str(st[0].stats.starttime + (st[0].stats.endtime - st[0].stats.starttime) * i / number_month)[0:10])

    tick_positions = [int(i * len(Sz[0]) / (len(dates) - 1)) for i in range(len(dates))]

    
    t2 = np.arange(0,len(t))
    t1 = np.arange(0,len(t))
    
    plt.set_cmap("viridis")
    norm_z = mpl.colors.Normalize(vmin=-190, vmax=-150)
    norm_p = mpl.colors.Normalize(vmin=-60, vmax=50)
    norm_coh = mpl.colors.Normalize(vmin=0, vmax=1)
    

    
    plt.figure(dpi = 300,figsize=(30,30))
    plt.suptitle(str(st[0].stats.network)+'.'+str(st[0].stats.station))
    plt.subplot(131)
    plt.pcolormesh(f, t2, 10*np.log10(Sp.T),norm=norm_p)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Time [Date]')
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.title("Spectrogram BDH", y=1.025)
    plt.vlines(f[f1], 0, t2[-1],linewidth=5,linestyle='dashed',color='black')
    plt.vlines(f[f2], 0, t2[-1],linewidth=5,linestyle='dashed',color='black')

    plt.yticks(tick_positions, dates, rotation=0)

    plt.tick_params(axis='both', which='minor', length=10, width=1)
    plt.tick_params(axis='both', which='major', length=10, width=2)
    
    # plt.colorbar()
    plt.colorbar(label="Pressure $(Pa^2/Hz)$[dB]",orientation='vertical')
    plt.tight_layout()
    
    plt.subplot(132)
    plt.pcolormesh(f, t2, 10*np.log10(Sz.T*(2*np.pi*f)**4),norm=norm_z)  
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.title("Spectrogram BHZ", y=1.025)
    plt.yticks([])
    plt.colorbar(label="Acceleration$((m/s^2)^2/Hz)$ [dB] ",orientation='vertical')
    plt.vlines(f[f1], 0, t2[-1],linewidth=5,linestyle='dashed',color='black')
    plt.vlines(f[f2], 0, t2[-1],linewidth=5,linestyle='dashed',color='black')

    plt.tick_params(axis='both', which='minor', length=10, width=1)
    plt.tick_params(axis='both', which='major', length=10, width=2)
    plt.subplot(133)
    plt.pcolormesh(f,np.arange(0,len(Czp)),Czp,norm = norm_coh)
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.title("Coherogram", y=1.025)
    cbar = plt.colorbar(label='Coherence', orientation='vertical')
    plt.yticks([])
    plt.tight_layout()
    cbar.set_ticks([round(tick, 1) for tick in cbar.get_ticks()])
    plt.vlines(f[f1], 0, nd,linewidth=5,linestyle='dashed',color='black')
    plt.vlines(f[f2], 0, nd,linewidth=5,linestyle='dashed',color='black')
    plt.ylim([0,nd])
    plt.tick_params(axis='both', which='minor', length=10, width=1)
    plt.tick_params(axis='both', which='major', length=10, width=2)
    # plt.subplot(144)
    # # plt.plot(Czp_filtered[:,20], np.arange(0, len(Czp)), 'b')   

    # f1 = np.argmin(np.abs(f-f_min))
    # f2 = np.argmin(np.abs(f-f_max))
    # plt.plot(np.median(Czp[:,f1:f2],axis=1), np.arange(0, len(Czp)),'b')
    # plt.plot(scipy.signal.savgol_filter(np.median(Czp[:,f1:f2],axis=1), 24, 1)*1, np.arange(0, len(Czp)),'black',linewidth=3)
    # plt.vlines(x=0.80, ymin=0, ymax=len(Czp),linestyles='dashed',color='r',label='0.80 Threshold',linewidth =3)
    # # plt.xlim([0.75,1])
    # plt.xlabel('Coherence')
    # plt.title(str(f_min) +' to ' +str(f_max) + ' Hz', y=1.025)
    # plt.grid(True)
    # plt.ylim([0,len(Czp)])
    # plt.legend(loc='lower left',fontsize=25)
    # plt.yticks([])
    # plt.tight_layout()
    
    Dzz_smoothed = 10*np.log10((scipy.signal.savgol_filter(np.median(Dzz[:,f1:f2]*(2*np.pi*f[f1:f2])**4,axis=1), 10, 1)))
    Dpp_smoothed = 10*np.log10(scipy.signal.savgol_filter(np.median(Dpp[:,f1:f2],axis=1), 10, 1))
    Czp_smoothed =scipy.signal.savgol_filter(np.median(Czp[:,f1:f2],axis=1), 10, 1)
    
    good_windows = []
    bad_windows = []
    
    
    for i in range(0,len(Dzz)):
        if 10*np.log10(np.median(Dzz[i,f1:f2]*(2*np.pi*f[f1:f2])**4)) > Tresh_Dz and 10*np.log10(np.median(Dpp[i,f1:f2])) > Tresh_Dp and np.median(Czp[i,f1:f2]) > Tresh_coh:
            good_windows.append(i)
        else:
            bad_windows.append(i)
            
    tick_positions_2 = [int(i * len(Dpp) / (len(dates) - 1)) for i in range(len(dates))]

    plt.figure(dpi = 300,figsize=(30,25))
    # plt.suptitle("Median [" +str(f_min) +' to ' +str(f_max) + ' Hz ]')
    plt.subplot(131)
    # plt.plot(Czp_filtered[:,20], np.arange(0, len(Czp)), 'b')   

    f1 = np.argmin(np.abs(f-f_min))
    f2 = np.argmin(np.abs(f-f_max))
    plt.plot(10*np.log10(np.median(Dpp[:,f1:f2],axis=1)), np.arange(0, len(Czp)),'b')
    plt.plot(10*np.log10(scipy.signal.savgol_filter(np.median(Dpp[:,f1:f2],axis=1), 10, 1)), np.arange(0, len(Czp)),'black',linewidth=5)
    # plt.vlines(x=Tresh_Dp, ymin=0, ymax=len(Czp),linestyles='dashed',color='r',label='10 dB Threshold',linewidth = 5)
    # plt.vlines(x=0.80, ymin=0, ymax=len(Czp),linestyles='dashed',color='r',label='0.80 Threshold',linewidth =3)
    # plt.xlim([0.75,1])
    plt.xlabel('$Pa^2/Hz$ [dB]')
    plt.ylabel('Time [Date]')
    plt.title("Spectrogram Pressure Data")
    # plt.legend(loc='upper right',fontsize=30)
    plt.grid(True)
    plt.ylim([0,len(Czp)])
    # plt.legend(loc='lower left',fontsize=25)
    plt.yticks(tick_positions_2, dates, rotation=60)
    plt.tight_layout()
    plt.xlim([20,40])

    plt.subplot(132)
    # plt.plot(Czp_filtered[:,20], np.arange(0, len(Czp)), 'b')   


    plt.plot(10*np.log10((np.median(Dzz[:,f1:f2]*(2*np.pi*f[f1:f2])**4,axis=1))), np.arange(0, len(Czp)),'b')
    plt.plot(10*np.log10((scipy.signal.savgol_filter(np.median(Dzz[:,f1:f2]*(2*np.pi*f[f1:f2])**4,axis=1), 10, 1))), np.arange(0, len(Czp)),'black',linewidth=5)
    # plt.vlines(x=0.80, ymin=0, ymax=len(Czp),linestyles='dashed',color='r',label='0.80 Threshold',linewidth =3)
    # plt.xlim([0.75,1])
    # plt.vlines(x=Tresh_Dz, ymin=0, ymax=len(Czp),linestyles='dashed',color='r',label='-160 dB Threshold',linewidth = 5)
    plt.xlabel('$(m/s^2)^2/Hz$ [dB]')

    plt.title("Spectrogram Seismic Data")
    plt.grid(True)
    plt.ylim([0,len(Czp)])
    # plt.legend(loc='lower left',fontsize=25)
    
    plt.yticks([])
    # plt.legend(loc='upper right',fontsize=30)

    # plt.yticks([])
    plt.tight_layout()
    plt.xlim([-180,-120])

    plt.subplot(133)
    # plt.plot(Czp_filtered[:,20], np.arange(0, len(Czp)), 'b')   

    f1 = np.argmin(np.abs(f-f_min))
    f2 = np.argmin(np.abs(f-f_max))
    plt.plot(np.median(Czp[:,f1:f2],axis=1), np.arange(0, len(Czp)),'b')
    plt.plot(scipy.signal.savgol_filter(np.median(Czp[:,f1:f2],axis=1), 10, 1), np.arange(0, len(Czp)),'black',linewidth=5)
    plt.vlines(x=Tresh_coh, ymin=0, ymax=len(Czp),linestyles='dashed',color='r',label='0.80 Threshold',linewidth = 5)
    for ii in range(0,len(good_windows)):
        plt.axhspan(good_windows[ii],good_windows[ii]+1, color='green', alpha=0.5)
    for ii in range(0,len(bad_windows)):
        plt.axhspan(bad_windows[ii],bad_windows[ii]+1, color='lightcoral', alpha=0.5)
    plt.xlabel('Coherency')
    plt.title('Coherogram')
    plt.grid(True)
    plt.ylim([0,len(Czp)])
    plt.xlim([0,1])
    # plt.legend(loc='lower left',fontsize=25)
    plt.yticks([])
    # plt.legend(loc='upper right',fontsize=30)

    plt.tight_layout()
    
    # Set the number of decimal places in the colorbar ticks
  
    # plt.figure(dpi = 300,figsize=(25,20))
    # plt.subplot(131)
    # plt.pcolormesh(f, t2, 10*np.log10(Sp.T))
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Time')
    # plt.xscale('log')
    # plt.xlim(0.005,0.1)    
    # plt.title("Spectrogram BDH")
    # # plt.yticks(np.arange(0,len(t),int(len(t)/3)), dates, rotation=45)  # Set x-axis names with rotation
    # # plt.colorbar()
    # plt.colorbar(label="Pressure $(Pa^2/Hz)$[dB]",orientation='vertical')
    # plt.tight_layout()
    
    # plt.subplot(132)
    # plt.pcolormesh(f, t2, 10*np.log10(Sz.T*(2*np.pi*f)**4))  
    # plt.xlabel('Frequency [Hz]')
    # plt.xscale('log')
    # plt.xlim(0.005,0.1)    
    # plt.title("Spectrogram BHZ")
    # # plt.yticks(np.arange(0,len(t),int(len(t)/3)), dates, rotation=45)  # Set x-axis names with rotation
    # plt.colorbar(label="Acceleration$((m/s^2)^2/Hz)$ [dB] ",orientation='vertical')

    # # plt.colorbar()
     
    # plt.subplot(133)
    # plt.pcolormesh(f,np.arange(0,len(Czp)),Czp)
    # plt.xscale('log')
    # plt.xlim(0.005,0.1)    
    # plt.xlabel('Frequency [Hz]')
    # plt.xscale('log')
    # plt.title("Coherogram")
    # plt.colorbar(label='Coherence', orientation='vertical')
    # plt.yticks(np.arange(0,len(t1),int(len(t1)/3)), dates, rotation=45)  # Set x-axis names with rotation
    plt.tight_layout()
#%%

def coherogram_spectrogram_all(st,nseg=2**12,data_indetvarl = 7):
    
    f,t,Sz = scipy.signal.spectrogram(st.select(component='Z')[0].data,fs=st[3].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2),window='hann')
    f,t,Sh1 = scipy.signal.spectrogram(st.select(component='1')[0].data,fs=st[3].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2),window='hann')
    f,t,Sh2 = scipy.signal.spectrogram(st.select(component='2')[0].data,fs=st[3].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2),window='hann')
    f,t,Sp = scipy.signal.spectrogram(st.select(component='H')[0].data,fs=st[0].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2),window='hann')
    
    
    number_month = (st[0].stats.endtime - st[0].stats.starttime) // (data_indetvarl*24*3600)
    
    
    dates = []
    for i in range(0, int(number_month)+1):
        # print(i)
        dates.append(str(st[0].stats.starttime + (st[0].stats.endtime - st[0].stats.starttime) * i / number_month)[0:10])

    tick_positions = [int(i * len(Sz[0]) / (len(dates) - 1)) for i in range(len(dates))]

    
    t2 = np.arange(0,len(t))
    
    plt.set_cmap("viridis_r")
    norm_z = mpl.colors.Normalize(vmin=-190, vmax=-150)
    norm_p = mpl.colors.Normalize(vmin=-60, vmax=50)    
    norm_h = mpl.colors.Normalize(vmin=-180, vmax=-80)    

    
    plt.figure(dpi = 300,figsize=(40,30))
    plt.suptitle(str(st[0].stats.network)+'.'+str(st[0].stats.station))
    plt.subplot(141)
    plt.pcolormesh(f, t2, 10*np.log10(Sp.T),norm=norm_p)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Time [Date]')
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.title("Spectrogram BDH", y=1.025)

    plt.yticks(tick_positions, dates, rotation=0)

    plt.tick_params(axis='both', which='minor', length=10, width=1)
    plt.tick_params(axis='both', which='major', length=10, width=2)
    
    # plt.colorbar()
    plt.colorbar(label="Pressure $(Pa^2/Hz)$[dB]",orientation='vertical')
    plt.tight_layout()
    
    plt.subplot(142)
    plt.pcolormesh(f, t2, 10*np.log10(Sz.T*(2*np.pi*f)**4),norm=norm_z)  
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.title("Spectrogram BHZ", y=1.025)
    plt.yticks([])
    plt.colorbar(label="Acceleration$((m/s^2)^2/Hz)$ [dB] ",orientation='vertical')

    plt.subplot(143)
    plt.pcolormesh(f, t2, 10*np.log10(Sh1.T*(2*np.pi*f)**4),norm=norm_h)  
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.title("Spectrogram BH1", y=1.025)
    plt.yticks([])
    plt.colorbar(label="Acceleration$((m/s^2)^2/Hz)$ [dB] ",orientation='vertical')
    
    plt.subplot(144)
    plt.pcolormesh(f, t2, 10*np.log10(Sh2.T*(2*np.pi*f)**4),norm=norm_h)  
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.title("Spectrogram BH2", y=1.025)
    plt.yticks([])
    plt.colorbar(label="Acceleration$((m/s^2)^2/Hz)$ [dB] ",orientation='vertical')

#%%

def coherogram_spectrogram_all(st,st1,st2,st3,nseg=2**11,tw =1):
    f_min = 0.008
    f_max = 0.020
    Tresh_coh = 0.8
    Tresh_Dz = -170
    Tresh_Dp = 25
    
    TP=5
    
    f,t,Sz = scipy.signal.spectrogram(st.select(component='Z')[0].data,fs=st[3].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2),window='hann')

    ws = int(tw*60*60 * st[0].stats.sampling_rate)
    Z ,nd = sliding_window(st.select(component='Z')[0].data, ws = ws,hann=True)
    Z1 ,nd = sliding_window(st1.select(component='Z')[0].data, ws = ws,hann=True)
    Z2 ,nd = sliding_window(st2.select(component='Z')[0].data, ws = ws,hann=True)
    Z3 ,nd = sliding_window(st3.select(component='Z')[0].data, ws = ws,hann=True)

    P ,nd = sliding_window(st.select(component='H')[0].data, ws = ws,hann=True)

    f,Czp = scipy.signal.coherence(Z,P,fs=st[0].stats.sampling_rate,nperseg =nseg,
                                   noverlap=(nseg*0.5),
                                   window=scipy.signal.windows.tukey(nseg,
                                   (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,Czp1 = scipy.signal.coherence(Z1,P,fs=st[0].stats.sampling_rate,nperseg =nseg,
                                   noverlap=(nseg*0.5),
                                   window=scipy.signal.windows.tukey(nseg,
                                   (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,Czp2 = scipy.signal.coherence(Z2,P,fs=st[0].stats.sampling_rate,nperseg =nseg,
                                   noverlap=(nseg*0.5),
                                   window=scipy.signal.windows.tukey(nseg,
                                   (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,Czp3 = scipy.signal.coherence(Z3,P,fs=st[0].stats.sampling_rate,nperseg =nseg,
                                   noverlap=(nseg*0.5),
                                   window=scipy.signal.windows.tukey(nseg,
                                   (TP*60*st[0].stats.sampling_rate)/nseg))
        
    
    number_month = (st[0].stats.endtime - st[0].stats.starttime) // (1*24*3600)
    
    
    f1 = np.argmin(np.abs(f-f_min))
    f2 = np.argmin(np.abs(f-f_max))
    
    
    dates = []
    for i in range(0, int(number_month)+1):
        # print(i)
        dates.append(str(st[0].stats.starttime + (st[0].stats.endtime - st[0].stats.starttime) * i / number_month)[0:10])

    tick_positions = [int(i * len(Czp) / (len(dates) - 1)) for i in range(len(dates))]

    
    t2 = np.arange(0,len(t))
    t1 = np.arange(0,len(t))
    
    plt.set_cmap("viridis")
    
    
    norm_z = mpl.colors.Normalize(vmin=-190, vmax=-150)
    norm_p = mpl.colors.Normalize(vmin=-60, vmax=50)
    norm_coh = mpl.colors.Normalize(vmin=0, vmax=1)
 

    plt.figure(dpi=300, figsize=(40, 30))  # Adjusted figsize to better suit the new layout
    plt.suptitle("Coherence Between Vertical Acceleration and Pressure")
    plt.subplot(141)
    plt.pcolormesh(f,np.arange(0,len(Czp)),Czp,norm = norm_coh)
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Time [Date]')

    plt.xscale('log')
    plt.title("Raw", y=1.025)
    # cbar = plt.colorbar(label='Coherence', orientation='vertical')
    cbar = plt.colorbar()

    plt.yticks([])
    plt.tight_layout()
    plt.yticks(tick_positions, dates, rotation=60)

    cbar.set_ticks([round(tick, 1) for tick in cbar.get_ticks()])
    plt.vlines(f[f1], 0, nd,linewidth=7,linestyle='dashed',color='red')
    plt.vlines(f[f2], 0, nd,linewidth=7,linestyle='dashed',color='red')
    plt.ylim([0,nd])
    
    plt.subplot(142)
    plt.pcolormesh(f,np.arange(0,len(Czp1)),Czp1,norm = norm_coh)
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.title("Tilt Reduction - Glitch Removal", y=1.025)
    # cbar = plt.colorbar(label='Coherence', orientation='vertical')
    cbar = plt.colorbar()

    plt.yticks([])
    plt.tight_layout()
    cbar.set_ticks([round(tick, 1) for tick in cbar.get_ticks()])
    plt.vlines(f[f1], 0, nd,linewidth=7,linestyle='dashed',color='red')
    plt.vlines(f[f2], 0, nd,linewidth=7,linestyle='dashed',color='red')
    plt.ylim([0,nd])
    
    
    plt.subplot(143)
    plt.pcolormesh(f,np.arange(0,len(Czp2)),Czp2,norm = norm_coh)
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.title("Event Removal", y=1.025)
    # cbar = plt.colorbar(label='Coherence', orientation='vertical')
    cbar = plt.colorbar()
    plt.yticks([])
    plt.tight_layout()
    cbar.set_ticks([round(tick, 1) for tick in cbar.get_ticks()])
    plt.vlines(f[f1], 0, nd,linewidth=7,linestyle='dashed',color='red')
    plt.vlines(f[f2], 0, nd,linewidth=7,linestyle='dashed',color='red')
    plt.ylim([0,nd])
    
    plt.subplot(144)
    plt.pcolormesh(f,np.arange(0,len(Czp3)),Czp3,norm = norm_coh)
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.title("Coherent Noise Removal", y=1.025)
    # cbar = plt.colorbar(label='Coherence', orientation='vertical')
    cbar = plt.colorbar()
    plt.yticks([])
    plt.tight_layout()
    cbar.set_ticks([round(tick, 1) for tick in cbar.get_ticks()])
    plt.vlines(f[f1], 0, nd,linewidth=7,linestyle='dashed',color='red')
    plt.vlines(f[f2], 0, nd,linewidth=7,linestyle='dashed',color='red')
    plt.ylim([0,nd])
    
    
    plt.subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.07, wspace=0.02, hspace=0.2)
    plt.tight_layout(pad=1.0, w_pad=0, h_pad=2.0)
    
    
#%%

def coherogram_spectrogram_alpha(st,nseg=2**12,tw =1):
    f_min = 0.007
    f_max = 0.020
    Tresh_coh = 0.8
    Tresh_Dz = -170
    Tresh_Dp = 25
    
    TP=5
    f,t,Sz = scipy.signal.spectrogram(st.select(component='Z')[0].data,fs=st[3].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2),window='hann')
    f,t,Sp = scipy.signal.spectrogram(st.select(component='H')[0].data,fs=st[0].stats.sampling_rate,nperseg = nseg,noverlap=(nseg/2),window='hann')
    

    ws = int(tw*60*60 * st[0].stats.sampling_rate)
    Z ,nd = sliding_window(st.select(component='Z')[0].data, ws = ws,hann=True)
    P ,nd = sliding_window(st.select(component='H')[0].data, ws = ws,hann=True)

    f,Czp = scipy.signal.coherence(Z,P,fs=st[0].stats.sampling_rate,nperseg =nseg,
                                   noverlap=(nseg*0.5),
                                   window=scipy.signal.windows.tukey(nseg,
                                   (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,Dzz = scipy.signal.welch(Z,fs=st[0].stats.sampling_rate,nperseg =nseg,
                                   noverlap=(nseg*0.5),
                                   window=scipy.signal.windows.tukey(nseg,
                                   (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,Dpp = scipy.signal.welch(P,fs=st[0].stats.sampling_rate,nperseg =nseg,
                                   noverlap=(nseg*0.5),
                                   window=scipy.signal.windows.tukey(nseg,
                                   (TP*60*st[0].stats.sampling_rate)/nseg))
    
    number_month = (st[0].stats.endtime - st[0].stats.starttime) // (7*24*3600)
    
    
    f1 = np.argmin(np.abs(f-f_min))
    f2 = np.argmin(np.abs(f-f_max))
    
    
    dates = []
    for i in range(0, int(number_month)+1):
        # print(i)
        dates.append(str(st[0].stats.starttime + (st[0].stats.endtime - st[0].stats.starttime) * i / number_month)[0:10])

    tick_positions = [int(i * len(Sz[0]) / (len(dates) - 1)) for i in range(len(dates))]

    
    t2 = np.arange(0,len(t))
    t1 = np.arange(0,len(t))
    
    plt.set_cmap("viridis")
    # plt.set_cmap("magma")
    # plt.set_cmap("inferno")
    # plt.set_cmap("plasma")
    
    
    norm_z = mpl.colors.Normalize(vmin=-190, vmax=-150)
    norm_p = mpl.colors.Normalize(vmin=-60, vmax=50)
    norm_coh = mpl.colors.Normalize(vmin=0, vmax=1)
 
    # plt.yticks([])
    # plt.tight_layout()
    
    Dzz_smoothed = 10*np.log10((scipy.signal.savgol_filter(np.median(Dzz[:,f1:f2]*(2*np.pi*f[f1:f2])**4,axis=1), 10, 1)))
    Dpp_smoothed = 10*np.log10(scipy.signal.savgol_filter(np.median(Dpp[:,f1:f2],axis=1), 10, 1))
    Czp_smoothed =scipy.signal.savgol_filter(np.median(Czp[:,f1:f2],axis=1), 10, 1)
    
    good_windows = []
    bad_windows = []
    
    
    for i in range(0,len(Dzz)):
        if 10*np.log10(np.median(Dzz[i,f1:f2]*(2*np.pi*f[f1:f2])**4)) > Tresh_Dz and 10*np.log10(np.median(Dpp[i,f1:f2])) > Tresh_Dp and np.median(Czp[i,f1:f2]) > Tresh_coh:
            good_windows.append(i)
        else:
            bad_windows.append(i)
            
    tick_positions_2 = [int(i * len(Dpp) / (len(dates) - 1)) for i in range(len(dates))]

   

    import matplotlib.gridspec as gridspec

    plt.figure(dpi=300, figsize=(40, 20))  # Adjusted figsize to better suit the new layout
    gs = gridspec.GridSpec(1, 6, width_ratios=[3, 1, 3, 1, 3, 1])  # Adjust ratios as needed
    
    # Subplot 1 (twice the width)
    ax1 = plt.subplot(gs[0])
    # Your plotting commands for subplot 1...
    plt.pcolormesh(f, t2, 10*np.log10(Sp.T),norm=norm_p)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Time [Date]')
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.title("Pressure Spectrogram", y=1.025)
    plt.vlines(f[f1], 0, t2[-1],linewidth=7,linestyle='dashed',color='red')
    plt.vlines(f[f2], 0, t2[-1],linewidth=7,linestyle='dashed',color='red')
    
    plt.yticks(tick_positions, dates, rotation=60)
     
     # plt.colorbar()
    plt.colorbar(label="Pressure $(Pa^2/Hz)$[dB]",orientation='vertical')
    plt.tight_layout()
    # Subplot 2 (half the width of subplot 1)
    ax2 = plt.subplot(gs[1])
    # Your plotting commands for subplot 2...
    f1 = np.argmin(np.abs(f-f_min))
    f2 = np.argmin(np.abs(f-f_max))
    plt.plot(10*np.log10(np.median(Dpp[:,f1:f2],axis=1)), np.arange(0, len(Czp)),'b')
    plt.plot(10*np.log10(scipy.signal.savgol_filter(np.median(Dpp[:,f1:f2],axis=1), 10, 1)), np.arange(0, len(Czp)),'black',linewidth=5)
    
    plt.xlabel('$Pa^2/Hz$ [dB]')
    plt.grid(True)
    plt.ylim([0,len(Czp)])
    plt.tight_layout()
    plt.xlim([20,40])
    plt.yticks([])
    # Subplot 3 (twice the width, next to subplot 2)
    ax3 = plt.subplot(gs[2])
    # Your plotting commands for subplot 3...
    plt.pcolormesh(f, t2, 10*np.log10(Sz.T*(2*np.pi*f)**4),norm=norm_z)  
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.title("Vertical Acceleration Spectrogram", y=1.025)
    plt.yticks([])
    plt.colorbar(label="Acceleration$((m/s^2)^2/Hz)$ [dB] ",orientation='vertical')
    plt.vlines(f[f1], 0, t2[-1],linewidth=7,linestyle='dashed',color='red')
    plt.vlines(f[f2], 0, t2[-1],linewidth=7,linestyle='dashed',color='red')
    
    
    # Subplot 4 (half the width of subplot 3)
    ax4 = plt.subplot(gs[3])
    # Your plotting commands for subplot 4...
    plt.plot(10*np.log10((np.median(Dzz[:,f1:f2]*(2*np.pi*f[f1:f2])**4,axis=1))), np.arange(0, len(Czp)),'b')
    plt.plot(10*np.log10((scipy.signal.savgol_filter(np.median(Dzz[:,f1:f2]*(2*np.pi*f[f1:f2])**4,axis=1), 10, 1))), np.arange(0, len(Czp)),'black',linewidth=5)
    plt.xlabel('$(m/s^2)^2/Hz$ [dB]')
    plt.grid(True)
    plt.ylim([0,len(Czp)])
    plt.yticks([])
    plt.tight_layout()
    plt.xlim([-180,-140])
    
    # Subplot 5 (twice the width, next to subplot 4)
    ax5 = plt.subplot(gs[4])
    # Your plotting commands for subplot 5...
    plt.pcolormesh(f,np.arange(0,len(Czp)),Czp,norm = norm_coh)
    plt.xscale('log')
    plt.xlim(0.005,0.1)    
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.title("Coherogram ( Pressure / Acceleration )", y=1.025)
    cbar = plt.colorbar(label='Coherence', orientation='vertical')
    plt.yticks([])
    plt.tight_layout()
    cbar.set_ticks([round(tick, 1) for tick in cbar.get_ticks()])
    plt.vlines(f[f1], 0, t2[-1],linewidth=7,linestyle='dashed',color='red')
    plt.vlines(f[f2], 0, t2[-1],linewidth=7,linestyle='dashed',color='red')
    
    plt.ylim([0,nd])
    # Subplot 6 (half the width of subplot 5)
    ax6 = plt.subplot(gs[5])
    # Your plotting commands for subplot 6...
    f1 = np.argmin(np.abs(f-f_min))
    f2 = np.argmin(np.abs(f-f_max))
    plt.plot(np.median(Czp[:,f1:f2],axis=1), np.arange(0, len(Czp)),'b')
    plt.plot(scipy.signal.savgol_filter(np.median(Czp[:,f1:f2],axis=1), 10, 1), np.arange(0, len(Czp)),'black',linewidth=5)
    plt.vlines(x=Tresh_coh, ymin=0, ymax=len(Czp),linestyles='dashed',color='r',label='0.80 Threshold',linewidth = 5)
    for ii in range(0,len(good_windows)):
        plt.axhspan(good_windows[ii],good_windows[ii]+1, color='green', alpha=0.6)
    for ii in range(0,len(bad_windows)):
        plt.axhspan(bad_windows[ii],bad_windows[ii]+1, color='white', alpha=0.5)
    plt.xlabel('Coherency')
    plt.grid(True)
    plt.ylim([0,len(Czp)])
    plt.xlim([0,1])
    plt.yticks([])
    plt.subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.07, wspace=0.02, hspace=0.2)
    plt.tight_layout(pad=1.0, w_pad=0, h_pad=2.0)

    
    # plt.savefig('/Users/mohammadamin/Desktop/spectrogram_coherogram.pdf' , format='pdf')
      
    
#%%
def coherence_significance_level(n_windows, prob=0.95):
    """
    Definition: L_1(alpha, q) = sqrt(1-alpha**(1/q))

    where alpha = 1-prob and 2(q+1) = nwinds (degree of freedom)

    For nwinds >> 1, L1 ~ sqrt(1-alpha**(2/nwinds))
    For a 95% signif level this comes out to
        sqrt(1-.05**(2/nwinds)) for nwinds >> 1.
    I previously used sqrt(2/nwinds) for the 95% signif level (alpha=0.05),
    but L1 is much closer to sqrt(6/nwinds).

    Args:
        n_windows (int): number of windows
        prob (float): significance level (between 0 and 1)
    """
    assert prob >= 0 and prob <= 1
    alpha = 1 - prob
    q = n_windows/2 - 1
    return np.sqrt(1 - alpha ** (1. / q))

#%%
def compliance_uncertainty(com, coh,nd):
    """
    Definition: compliance *[sqr(1-coh(w))/(abs(coh(w))*sqr(2nd))]
    return the uncertainty of compliance function that can be used to plot error bars 
    for cumpliance function.

    Args:
        coh: Coherence function between vertical channel and pressure channel.
        com = Compliance function
        nd (int): number of windows
    """
    
    return(com *[np.sqrt(1-coh)/(abs(coh*np.sqrt(2*nd)))])
    
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
def compl(st,depth,nseg=2**13):
    f,Czp = scipy.signal.coherence(st.select(component='Z')[0].data,
                               st.select(component='H')[0].data,
                               fs=st[0].stats.sampling_rate,
                               nperseg =nseg,noverlap=(nseg/8))
    
    f,Dz = scipy.signal.welch(st.select(component='Z')[0],fs=st[0].stats.sampling_rate,
                              nperseg = nseg, noverlap=(nseg/2))
    f,Dp = scipy.signal.welch(st.select(component='H')[0],fs=st[0].stats.sampling_rate,
                              nperseg = nseg, noverlap=(nseg/2))

    k = wavenumber(2*np.pi*f, depth)
    
    com = k * Czp * np.sqrt(Dz/Dp)

    return com
    
#%%
def plot(stream1,stream2):
    
    plt.figure(dpi=300)
    plt.plot(stream1.select(channel='*Z')[0],color='gray')
    plt.plot(stream2.select(channel='*Z')[0],color='black')
#%%
def rms(d):
    '''
    Compute RMS of Data vector
    
    srt(mean(data^2))

    Parameters
    ----------
    d : Data

    Returns
    -------
    RMS of Data

    '''
    return(np.sqrt(np.mean(d**2)))    
#%%

def velp(vs , p = 0.25):
    """
    Caclculate Primary Velocity from possion ratio and Shear Velocity
    vp = vs * np.sqrt((1-p)/(0.5-p))
    https://wiki.seg.org/wiki/Poisson%27s_ratio
    Parameters
    ----------
    vs : Shear Velocity (m/s)
    
    p : poisson ratio . The default is 0.25.

    Returns
    -------
    vp : primary veocity m/s

    """
    vp = vs * np.sqrt((1-p)/(0.5-p))
    
    return(vp)
    

#%%

def density(vs, p = 0.25):
    """
    Calculate Density from shear velocity by using Gardner Eq.
    
    Density  = 0.31*(Vp**0.25)
    
    https://www.subsurfwiki.org/wiki/Gardner%27s_equation
    Parameters
    ----------
    vs : shear velocity
    
    p = poisson ratio. The default is 0.25

    Returns
    -------
    
    Density : g/cc

    """
    
    Density = 1000* 0.31 *((vs * np.sqrt((1-p)/(0.5-p)))**0.25)
    
    return(Density)
    
#%%
def misfit(d,m,l=2,s=1):
    '''
        
    Calculate misfit (L2)

    Parameters
    ----------
    d : Measured Data
    
    m : Modeled Data
    
    l : power of the norm, default = 2.
       
    s : Estimated uncertainty
        The default is 1.    
    Returns
    -------
    misfit
    '''
    misfit = np.sum(((d-m)/s)**l)
    return(misfit)

#%%
def liklihood(d,m,k=1,s=1):
    '''
    Calculate Linkihood by using gaussian misfit
    Tarantola paper-monte carlo
    Parameters
    ----------
    d : Measured Data
    
    m : Modeled Data
    
    k : I don'y know what is it! 
        The default is 1.
        
    s : Estimated uncertainty
        The default is 1.
    Returns
    -------
    likilihood
    '''
    L =   k * np.exp(-0.5*np.sum((d-m)**2/(s)**2))
    
    return(L)
    
#%%
def plt_params():
    plt.rcParams['font.size'] = 50
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = ['STIXGeneral']
    plt.rcParams['font.weight'] = 'normal'
    plt.rcParams['mathtext.default'] = 'regular'
    
    plt.rcParams['axes.grid'] = True
    plt.rcParams['axes.grid.which'] = 'major'
    plt.rcParams['axes.labelcolor'] = 'black'
    plt.rcParams['axes.labelpad'] = 4.0
    plt.rcParams['axes.labelsize'] = 50
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
    plt.rcParams['xtick.labelsize'] = 50
    plt.rcParams['xtick.major.pad'] = 2.0
    plt.rcParams['xtick.minor.pad'] = 2.0
    plt.rcParams['xtick.minor.visible'] = True
    plt.rcParams['ytick.color'] = 'black'
    plt.rcParams['ytick.direction'] = 'out'
    plt.rcParams['ytick.labelsize'] = 50
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


        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    