#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 11:11:57 2023

@author: mohammadamin
"""

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

def Calculate_Compliance(stream,f_min_com = 0.007,f_max_com = 0.017,gain_factor=0.66,time_window=1):
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
    split_streams = cut_stream_with_overlap(stream, time_window*60*60, overlap = ((time_window*60)-0.5)*60)
    
    f, Dpp = scipy.signal.welch(split_streams[0].select(component='H')[0], fs=split_streams[0][0].stats.sampling_rate,
                                          nperseg=nseg, noverlap=(nseg*0.9),
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

    coherence_mask = (f >= 0.008) & (f <= 0.016)
    
    ppsd_mask = (f >= 0.03) & (f <= 0.08)
    
    mask_dz = (f >= Fc2) & (f <= Fc1)

    # Treshhold to remove exclude those compiance function that are lesser or greater than median of all funtion +- percentage
    # percentage = 0.1
    
    print("Data optipization ...")
    print("Removing gravitatinal Attraction of ocean surface waves...")
    print("Computing Compliance funtion and admittance funtion ")
    for i in range(0,len(Czp)):
        
        # if np.mean(Czp[i][coherence_mask]) > 0.95 and (1-percentage)*np.mean(Com_Admitance[:,coherence_mask]) < np.mean(Com_Admitance[i][coherence_mask] < (1+percentage)*np.mean(Com_Admitance[:,coherence_mask])) :
        if np.median(Czp[i][coherence_mask]) > 0.80 and np.mean(Czp[i][ppsd_mask]) < 0.8 and np.mean(Dp[i][ppsd_mask]) < 1  and np.mean(Dz[i][ppsd_mask]) > 10e-17 and np.mean(Dz[i][mask_dz]) < 10e-12 :
        # if np.mean(Czp[i][coherence_mask]) > 0.99:
            
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
            
    plt.rcParams.update({'font.size': 25})
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
    plt.rcParams.update({'font.size': 25})
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
    
    eq_spans = tiskit.TimeSpans.from_eqs(stream.select(channel='*Z')[0].stats.starttime,
                                         stream.select(
                                             channel='*Z')[0].stats.endtime,
                                         minmag=6, days_per_magnitude=1.5)
    
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
    st = stream
    
    t1 = np.arange(0,len(angle))
    
    dates = []
    for i in range(0, int(len(t1) // (24*30))):
        print(i)
        dates.append(str(st[0].stats.starttime + (st[0].stats.endtime - st[0].stats.starttime) * i / int(len(t1) // (24*30)))[0:10])
    
    
    # Assuming t1 is a list that contains the time values
    # Calculate the positions of the ticks to align with the dates
    tick_positions = [int(i * len(t1) / (len(dates) - 1)) for i in range(len(dates))]
    
    plt.rcParams.update({'font.size': 40})
    plt.figure(dpi=300, figsize=(30, 20))
    plt.subplot(211)
    plt.title("YV." + str(stream[0].stats.station) +'  '+ str(stream[0].stats.starttime)[0:10]+'--'+str(stream[0].stats.endtime)[0:10]+" Tilt ")
    plt.plot(azimuth, '.', color='black')
    for i in range(0,len(eq_spans)):
        plt.plot(int((eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)),
                  azimuth[int(eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)],'o', color='red', markersize=10)

    plt.plot(int((eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)),
                  azimuth[int(eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)],'o', color='red', markersize=10,label="Event")
    # plt.xlabel("Time [Hour]")
    plt.xticks([])
    plt.ylabel("Azimuth ["u"\u00b0]")
    plt.legend(loc="lower right", facecolor='lightgreen')
    plt.grid(True)
    
    plt.subplot(212)
    plt.plot(angle, '.', color='black')
    
    for i in range(0,len(eq_spans)):
        plt.plot(int((eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)),
                  angle[int(eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)],'o', color='red', markersize=10)

    plt.plot(int((eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)),
                  angle[int(eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)],'o', color='red', markersize=10,label="Event")       
    plt.xlabel("Time [Date]")
    plt.ylabel("Incident Angle ["u"\u00b0]")
    plt.grid(True)
    plt.legend(loc="lower right", facecolor='lightgreen')
    
    # Set x-axis names with rotation and align ticks with the dates generated
    plt.xticks(tick_positions, dates, rotation=65)
    plt.ylim([-5, 5])
    plt.tight_layout()

    
    # print("Merging Stream ...")
    # rotated_stream = stream.copy()
    # rotated_stream.clear()
    # for i in range(0,len(split_streams)) : 
    #     rotated_stream = rotated_stream + split_streams[i]
    # rotated_stream.merge(fill_value='interpolate')
    # return(rotated_stream,azimuth,angle)
    return(azimuth,angle)

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
def optimizer(Com,Czp,Stream,f,alpha=0.95,beta = 0.70,f_min_com=0.008,f_max_com=0.015):
    coherence_mask = (f >= f_min_com) & (f <= f_max_com)
    coherence_mask2 = (f >= 0.02) & (f <= 0.025)
    High_Czp = [] 
    High_Com = []
    High_Com_Stream = []
    for i in range(0,len(Czp)):
        
        if np.mean(Czp[i][coherence_mask]) > alpha:
            if np.mean(Czp[i][coherence_mask2]) < beta:

            # Com1 = (k * Czp[i]* (np.sqrt(Dz[i]))) / (np.sqrt(Dp[i]) / gain_factor)
                        
                High_Czp.append(Czp[i])
     
                High_Com.append(Com[i])

                High_Com_Stream.append(Stream[i])
            
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











        