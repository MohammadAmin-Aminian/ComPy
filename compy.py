#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 11:13:05 2023

@author: mohammadamin
"""

import numpy as np
# import ffplot as fp
import matplotlib.pyplot as plt
import tiskit
# from obspy import read
from obspy.clients.fdsn import Client
import scipy 

def Calculate_Compliance(stream,f_min_com = 0.008,f_max_com = 0.015,gain_factor=0.66,time_window=1):
    nseg = 2**12
    TP = 5
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
    
    f, Dp = scipy.signal.welch(split_streams[0].select(component='H')[0], fs=split_streams[0][0].stats.sampling_rate,
                                          nperseg=nseg, noverlap=(nseg*0.9),
                                          window=scipy.signal.windows.tukey(nseg,
                                                                            (TP*60*split_streams[0][0].stats.sampling_rate)/nseg),
                                          average='median')    
    k = wavenumber((2*np.pi*f), -invz[0][0][0].elevation)



    Dp = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    Dz = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    Czp = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    Dzp = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    Com = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    Com_Admitance = np.zeros([len(split_streams),int(nseg / 2 + 1)])
    
    print('Calculating Compliance Funtion')
    
    for i in range(0,len(split_streams)):

        f, Dp[i] = scipy.signal.welch(split_streams[i].select(component='H')[0], fs=split_streams[i][0].stats.sampling_rate,
                                          nperseg=nseg, noverlap=(nseg*0.9),
                                          window=scipy.signal.windows.tukey(nseg,
                                                                            (TP*60*split_streams[i][0].stats.sampling_rate)/nseg),
                                          average='median')

        f, Dz[i] = scipy.signal.welch(split_streams[i].select(component='Z')[0], fs=split_streams[i][0].stats.sampling_rate,
                                      nperseg=nseg, noverlap=(nseg*0.9),
                                      window=scipy.signal.windows.tukey(nseg,
                                                                       (TP*60*split_streams[i][0].stats.sampling_rate)/nseg),
                                      average='median')

        f, Czp[i] = np.sqrt(scipy.signal.coherence(split_streams[i].select(component='Z')[0].data,
                                           split_streams[i].select(component='H')[
                                               0].data,
                                           fs=split_streams[i][0].stats.sampling_rate,
                                           nperseg=nseg, noverlap=(nseg*0.9),
                                           window=scipy.signal.windows.tukey(nseg,
                                                                             (TP*60*split_streams[i][0].stats.sampling_rate)/nseg)))


        f , Dzp[i] = scipy.signal.csd(split_streams[i].select(component='Z')[0],split_streams[i].select(component='H')[0],
                               fs=split_streams[i][0].stats.sampling_rate,
                               nperseg=nseg, noverlap=(nseg*0.9),
                               window=scipy.signal.windows.tukey(nseg,
                                                                 (TP*60*split_streams[i][0].stats.sampling_rate)/nseg))
        Com[i] = (k * Czp[i]* np.sqrt(np.abs(Dz[i]) / np.abs(Dp[i]))) / gain_factor
        Com_Admitance[i] = k * (Dzp[i]/Dp[i]) / gain_factor
        
        print( len( split_streams ) - i )
        
    High_Czp = [] 
    High_Dz = []
    High_Dp = [] 
    High_Com = []
    High_Com_Admitance = []
    High_Com_Stream = []


    coherence_mask = (f >= f_min_com) & (f <= f_max_com)
    
    coherence_mask_dp = (f >= 0.03) & (f <= 0.08)

    # Treshhold to remove exclude those compiance function that are lesser or greater than median of all funtion +- percentage
    percentage = 0.1
    
    for i in range(0,len(Czp)):
        
        # if np.mean(Czp[i][coherence_mask]) > 0.95 and (1-percentage)*np.mean(Com_Admitance[:,coherence_mask]) < np.mean(Com_Admitance[i][coherence_mask] < (1+percentage)*np.mean(Com_Admitance[:,coherence_mask])) :
        if np.mean(Czp[i][coherence_mask]) > 0.96 and np.mean(Dp[i][coherence_mask_dp]) < 1 :
         
            High_Czp.append(Czp[i])
    
            High_Dz.append(Dz[i])
         
            High_Dp.append(Dp[i]*(gain_factor**2))
            
            High_Com.append(Com[i])

            High_Com_Admitance.append(Com_Admitance[i])
         
            High_Com_Stream.append(split_streams[i])
        
            print(i)
    Fc1 = np.sqrt(9.8/(2*np.pi*0.5*-invz[0][0][0].elevation))
    Fc2 = np.sqrt(9.8/(2*np.pi*2*-invz[0][0][0].elevation))
    
    plt.rcParams.update({'font.size': 25})
    plt.figure(dpi=300,figsize=(14,16))
    plt.subplot(411)                                   
    for i in range(0,len(High_Dz)):
        plt.semilogx(f,10*np.log10(High_Dz[i]*(2*np.pi*f)**4),linewidth = 0.5,color='r')
    plt.semilogx(f,10*np.log10(np.median(High_Dz*(2*np.pi*f)**4,axis=0)),linewidth = 2,color='b',label='Median')
    plt.vlines(x = f_min_com, ymin=-180, ymax=-80,color='black',linestyles="dashed",label="High Coherence Band")
    plt.vlines(x = f_max_com, ymin=-180, ymax=-80,color='black',linestyles="dashed")
    plt.vlines(x = Fc1, ymin=-180, ymax=80,color='green',linestyles="solid",label="Maximum Frequency of IG")
    plt.vlines(x = Fc2 , ymin=-180, ymax=80,color='green',linestyles="solid")
    
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Vertical Acc [m/s^2] dB')
    plt.text(0.01, 0.8, 'a)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
    plt.grid(True)
    plt.legend(loc='upper right',fontsize=17)
    plt.xlim([0.001,1])
    plt.ylim([-170,-95])
  
    plt.subplot(412)                                   
    for i in range(0,len(High_Dp)):
        plt.semilogx(f,10*np.log10(High_Dp[i]),linewidth = 0.5,color='r')
    plt.semilogx(f,10*np.log10(np.median(High_Dp,axis=0)),linewidth = 2,color='b',label='Median')
    plt.vlines(x = f_min_com, ymin=-20, ymax=100,color='black',linestyles="dashed",label="High Coherence Band")
    plt.vlines(x = f_max_com, ymin=-20, ymax=100,color='black',linestyles="dashed")
    plt.vlines(x = Fc1, ymin=-20, ymax=100,color='green',linestyles="solid",label="Maximum Frequency of IG")
    plt.vlines(x = Fc2 , ymin=-20, ymax=100,color='green',linestyles="solid")
    
    
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Pressure')
    plt.text(0.01, 0.8, 'b)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
    plt.grid(True)
    plt.legend(loc='upper right',fontsize=17)
    plt.xlim([0.001,1])
    plt.ylim([-20,80])

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
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Compliance')
    plt.text(0.01, 0.8, 'd)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
    plt.grid(True)        
    plt.xlim([f_min_com,f_max_com])
    plt.ylim([10e-13,10e-10])
    plt.tight_layout()
    plt.legend(loc='upper right',fontsize=17)

    indices = np.where((f >= f_min_com) & (f <= f_max_com))
    f_c = f[indices]
    High_Com_c = []

    for i in range(0,len(High_Com)):
        High_Com_c.append(High_Com[i][indices])
    
    uncertainty = f_c.copy()
    High_Com_c = np.array(High_Com_c)
    

    for i in range(0,len(f_c)):
        uncertainty[i] = np.max(High_Com_c[:,i]) - np.min(High_Com_c[:,i])
    
    return(High_Com_c,f_c,uncertainty)
#%%
def Rotate(stream,time_window = 2):
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
    
    invp = client.get_stations(
        network=net,
        station=sta,
        channel="BDH",
        location="*",
        level="response")
#Split the data into pieces    
    print("Splitting The stream into "+ str(time_window)+"-Hour" + ' Windows')
    print("...")
    split_streams = split_stream(stream, duration = time_window*60*60)

    azimuth = np.zeros([len(split_streams)])
    angle = np.zeros([len(split_streams)])
    
    print("Reducing Tilt Effect")
    for i in range(0,len(split_streams)):
        try:
            print(len(split_streams) - i)
            D = tiskit.CleanRotator(split_streams[i], remove_eq=False)
    
            azimuth[i] = D.azimuth
            angle[i] =  D.angle
    
            split_streams[i] = D.apply(split_streams[i])

    # Transfer Function
            F = tiskit.DataCleaner(split_streams[i], ["*1","*2"], n_to_reject=0)
            split_streams[i] = F.clean_stream(split_streams[i])
            
            split_streams[i][0].stats.location = '00'
            split_streams[i][1].stats.location = '00'
            split_streams[i][2].stats.location = '00'
            split_streams[i][3].stats.location = '00'
        except Exception as e:
                print(f"Error occurred while processing item: {i}")
    
    print("Removing Instrument Response")
    print('...')
    for i in range(0,len(split_streams)):
    
        split_streams[i].select(channel="*Z").remove_response(inventory=invz,
                                                              output="DISP", plot=False)

        split_streams[i].select(channel="*H").remove_response(inventory=invp,
                                                              output="DEF", plot=False)
    plt.rcParams.update({'font.size': 25})
    plt.figure(dpi=300,figsize=(25,10))
    plt.subplot(211)
    plt.title("Tilt Correction Every " +str(time_window)+"H [" +str(stream[0].stats.station)+']')

    plt.plot(azimuth,color = 'b',linewidth = 3 , label = 'Azimuth')
    plt.xlabel("Time [Hour]")
    plt.ylabel("Degree ["u"\u00b0]")
    plt.legend(loc="upper right")
    plt.grid(True)
    
    plt.subplot(212)
    plt.plot(angle,color = 'b',linewidth = 3 , label = 'Angle')
    plt.xlabel("Time [Hour]")
    plt.ylabel("Degree ["u"\u00b0]")
    plt.grid(True)
    plt.legend(loc="upper right")
    plt.tight_layout()
    
    print("Merging Stream ...")
    rotated_stream = stream.copy()
    rotated_stream.clear()
    for i in range(0,len(split_streams)) : 
        rotated_stream = rotated_stream + split_streams[i]
    rotated_stream.merge(fill_value='interpolate')
    return(rotated_stream,azimuth,angle)
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

