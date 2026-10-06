#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 16:33:13 2023

@author: Mohammad-Amin Aminian

DFG Calibration

"""
from obspy.signal.trigger import plot_trigger
# from obspy.signal.trigger import coincidence_trigger
# from obspy.signal.trigger import classic_sta_lta
import matplotlib.pyplot as plt
import scipy
from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
import numpy as np
import obspy
import tiskit
from disba import PhaseDispersion

def calculate_spectral_ratio(stream,mag = 7,coh_trsh=0.97,mean_trsh = 0.97,f_min=0.02, f_max=0.06,plot_condition = False):
    '''
    Parameters
    ----------
    stream : TYPE
        Stream Before Removeing instrument response [raw Stream].
    mag : TYPE, optional
        Minimumm Magnitude for earthquakes to use in calulation. The default is 7.
    f_min : TYPE, optional
        Low frequency cornerof band of intrest . The default is 0.02.
    f_max : TYPE, optional
        High frequency cornerof band of intrest. The default is 0.06.

    coh_trsh = Coherence Treshhold to accept eqrthquakes. The default is 0.97.
    
    mean_trsh = Mean Treshhold to accept eqrthquakes. The default is 0.97.

    Returns
    -------
    ratio : TYPE

    Calibration Ratio

    '''
    # Frequncy to filter the data
    freq1=0.005 
    freq2=0.1
    
    # rho = Density of water 
    rho=1028
    
    #nseg = number of segment for coherence and transfer function calculation,overlap is 90%
    nseg = 2**9
    
    #TP = Tapering paramer 5 minuntes
    TP = 5 
    client = Client("RESIF")

    stream2 = stream.copy()
    stream2.clear()
    stream22 = stream.copy()
    stream22.clear()
    
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
    
    print ("Downloading Earthquakes with magnitude greater than " + str(mag) +" Mw \n ...."  )
    
    eq_spans = tiskit.TimeSpans.from_eqs(stream.select(channel='*Z')[0].stats.starttime,
                                     stream.select(
                                         channel='*Z')[0].stats.endtime,
                                     minmag=mag, days_per_magnitude=0.5)
    
    print ( str(len(eq_spans)) +" Earthquakes with magnitude greater than " + str(mag) +" Mw has been found"  )

    for i in range(0,len(eq_spans)):
        stream1= stream.copy()
        stream22 =stream22 + stream1.trim(eq_spans.start_times[i],eq_spans.start_times[i]+2*3600)
    
    stream22.sort(['starttime','channel'])
    
    t1 = np.zeros([int(len(stream22)/4),1])
    t2 = np.zeros([int(len(stream22)/4),1])


    print ("Removing the instrument response ...."  )

    stream2.select(channel="*Z").remove_response(inventory=invz,
                                            output="ACC", plot=False)

    stream2.select(channel="*H").remove_response(inventory=invp,
                                            output="DEF", plot=False)
    
    stream22.select(channel="*Z").remove_response(inventory=invz,
                                            output="ACC", plot=False)

    stream22.select(channel="*H").remove_response(inventory=invp,
                                            output="DEF", plot=False)
    
    # Retain the intended 0.005--0.1 Hz band.
    # Applying lowpass(freq1) followed by highpass(freq2) leaves no passband.
    stream22.filter("highpass", freq=freq1)
    stream22.filter("lowpass", freq=freq2)

    stream22.sort(['starttime','channel'])

    print(" Estimating rayleigh wave arrival time ....")
    
    for i in range(0,int(len(stream22)/4)):
        sst,t1[i],t2[i] = rayleigh_arrival(stream22[i*4:(i+1)*4],timelag=-5,window=20,plot_condition=False)
        stream2 = sst + stream2
    stream2.sort(['starttime','channel'])
    
    # Tapring  Hann with 5 minutes as default value for tapring lenght
    stream2.taper(max_percentage=0.1,type='hann',max_length=(TP*60))
    
    #detrend in 3 steps 
    stream2.detrend('linear')
    stream2.detrend('demean')
    stream2.detrend('simple')
    
    ratio_psd = np.zeros([len(eq_spans),int(nseg/2 + 1)])

    Czp = np.zeros([len(eq_spans),int(nseg/2 + 1)])
    Gpp = np.zeros([len(eq_spans),int(nseg/2 + 1)])
    Gzz = np.zeros([len(eq_spans),int(nseg/2 + 1)])

    # freqs = np.fft.fftfreq(len(stream2.select(channel="*Z")[2].data),d=stream2.select(channel="*Z")[2].stats.sampling_rate)
    
    # Calculation of spectral ratio (Transfer Function) and coherence with wlech method 
    
    print("Calculating Coherence and Spectral ratio")
    for i in range(0,len(eq_spans)-1):
        
        # You should find a solution for this line [ if statement!!], its a mess maybe add nfft to the caclulation ...!!?
        # or add zeros to the matrix !!!????
        if len(stream2.select(component='Z')[i].data) == 2501:
            f,Czp[i] = scipy.signal.coherence(stream2.select(component='Z')[i].data,
                                         stream2.select(component='H')[i].data,
                                         fs=stream2[i].stats.sampling_rate,
                                         nperseg =nseg,noverlap=(nseg*0.5),
                                         window=scipy.signal.windows.tukey(nseg,
                                         (TP*60*stream2[i].stats.sampling_rate)/nseg))

            f,Gzz[i] = scipy.signal.welch(stream2.select(component='Z')[i].data,
                                     fs=stream2[i].stats.sampling_rate,
                                     nperseg =nseg,noverlap=(nseg*0.5),
                                     window=scipy.signal.windows.tukey(nseg,
                                     (TP*60*stream2[i].stats.sampling_rate)/nseg))
    
            f,Gpp[i] = scipy.signal.welch(stream2.select(component='H')[i].data,
                                     fs=stream2[i].stats.sampling_rate,
                                     nperseg =nseg,noverlap=(nseg*0.5),
                                     window=scipy.signal.windows.tukey(nseg,
                                     (TP*60*stream2[i].stats.sampling_rate)/nseg))
                        
            ratio_psd[i] = Gpp[i] / Gzz[i]

            # print(i)
        else:
            pass    # Calculate the spectral ratio
            
    
    # Finding good windows [High coherence earthquakes]
    
    coherence_mask = (f >= f_min) & (f <= f_max)
    
    print ("Selecting Earthquakes with median of Coherence greater than "+ str(coh_trsh) +" and mean greater than " + str(mean_trsh) )
   
    High_Czp = [] 
    High_ratio_psd = []
    for i in range(0,len(eq_spans)):
                if np.median(Czp[i][coherence_mask]) >  coh_trsh and np.mean(Czp[i][coherence_mask]) >  mean_trsh:
                    High_Czp.append(Czp[i])
                    High_ratio_psd.append(ratio_psd[i])
                    # print(i)
    # i = 22 pakistan earthquake in RR29 data
    if plot_condition:
        for i in range(0,len(eq_spans)):
            if np.median(Czp[i][coherence_mask]) >  coh_trsh :         
                plt.rcParams.update({'font.size': 25})
                plt.figure(dpi=300,figsize=(25,20))
        
                # Plot 2: Acceleration
                plt.subplot(511)
                plt.title(stream22.select(channel="*Z")[i].stats.network+"."+stream22.select(channel="*Z")[i].stats.station+"  "+
                          str(stream22.select(channel="*Z")[i].stats.starttime)[0:19])    
                plt.plot(stream22.select(channel="*Z")[i].times(), stream22.select(channel="*Z")[i].data,linewidth=3,color='blue')
                plt.vlines(t1[i], np.min(stream22.select(channel="*Z")[i].data),
                           np.max(stream22.select(channel="*Z")[i].data),color='r',linewidth=3,linestyles='dashed')
                plt.vlines(t2[i], np.min(stream22.select(channel="*Z")[i].data),
                           np.max(stream22.select(channel="*Z")[i].data),color='r',linewidth=3,linestyles='dashed')
                
                plt.xlabel('Time (s)')
                plt.ylabel('Vertical Acc[m/s^2]')
                plt.grid(True)
                plt.text(0.01, 0.8, 'a)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
                
                # Plot 2: Pressure
                plt.subplot(512)
                plt.plot(stream22.select(channel="*H")[i].times(), stream22.select(channel="*H")[i].data,linewidth=3,color='blue')
                plt.vlines(t1[i], np.min(stream22.select(channel="*H")[i].data), 
                           np.max(stream22.select(channel="*H")[i].data),color='r',linewidth=3,linestyles='dashed')
                plt.vlines(t2[i], np.min(stream22.select(channel="*H")[i].data), 
                           np.max(stream22.select(channel="*H")[i].data),color='r',linewidth=3,linestyles='dashed')
                
                plt.xlabel('Time (s)')
                plt.ylabel('Pressure [pa] ')
                # plt.title('Stream2')
                plt.grid(True)
                plt.text(0.01, 0.8, 'b)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
                
            
                plt.subplot(513)
                plt.plot(stream2.select(channel="*H")[i].times(),-stream2.select(channel="*H")[i].data/np.max(stream2.select(channel="*H")[i].data)
                         ,label='Vertical Acc',color='blue')
                plt.plot(stream2.select(channel="*Z")[i].times(),stream2.select(channel="*Z")[i].data/np.max(stream2.select(channel="*Z")[i].data)
                         ,label='pressure',color='r',linestyle='dashed')
                # plt.xlim([9.7*10e4,9.8*10e4])
                plt.xlabel('Time[s]')
                plt.ylabel('Normilized')
                plt.legend(loc='upper right')
                plt.text(0.01, 0.8, 'c)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
                
                #    Plot 3: Coherence
                plt.subplot(514)
                plt.semilogx(f,Czp[i],linewidth=3,color='blue')
                plt.vlines(x = f_min, ymin=0, ymax=1,color='r',linestyles="dashed",label="Frequency limits")
                plt.vlines(x = f_max, ymin=0, ymax=1,color='r',linestyles="dashed")
                
                plt.ylabel("Coherence ")
                plt.xlabel("Frequency [Hz] ")
                plt.grid(True)
                plt.ylim([0,1])
                plt.xlim([freq1,freq2])
                plt.grid(True)
                plt.text(0.01, 0.8, 'd)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
                plt.legend(loc='upper right')
                
                #    Calculate spectral ratio
                
                # Plot 4: Spectral Ratio/Frequency
                plt.subplot(515)
                # positive_freq_indices = np.where((freqs) >= 0)
                
                #    Plot the right-sided spectrum
                # plt.semilogx(freqs[positive_freq_indices],(ratio[i][positive_freq_indices]),linewidth=2,color='blue')
                plt.semilogx(f, np.sqrt(ratio_psd[i]),linewidth=3,color='blue',label="Measured")
                plt.vlines(x = f_min, ymin=-10e10, ymax=10e10,color='r',linestyles="dashed",label="Frequency limits")
                plt.vlines(x = f_max, ymin=-10e10, ymax=10e10,color='r',linestyles="dashed")
                
                plt.xlim([freq1,freq2])
                plt.ylim([-10e6,10e7])
                plt.legend(loc='upper right')
                plt.xlabel('Frequency (Hz)')
                plt.ylabel('Spectral Ratio')
                plt.title('Spectral Ratio ')
                plt.grid(True)
                plt.text(0.01, 0.8, 'e)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
                
                plt.tight_layout()
                plt.show()
                # plt.savefig( str(stream22.select(channel="*Z")[i].stats.starttime)[0:19] + '.png',dpi=300)
                # plt.clf()
    
    # Calculate the spectral ratio
    Data = np.median(np.sqrt(High_ratio_psd)/(-invz[0][0][0].elevation*rho),axis=0)
    
    Data_zero = np.where((f >= f_min) & (f <= f_max), Data, 0)
    
    pvel = calculate_speed_of_sound_in_water(depth= - invz[0][0].elevation)
    
    f_dispersion_curve, Model = theoretical_p_a_ratio(alpha=pvel,h=-invz[0][0][0].elevation,t=np.sort(1/f[::-1][0:256]),plot_condition=False)
    
    # f_dispersion_curve, Model = phase_dispersion(t=np.sort(1/f[::-1][0:512]))
    
    Model_zero = np.where((f_dispersion_curve >= f_min) & (f_dispersion_curve <= f_max), Model, 0)

    gain_factor = grid_search(Data_zero[1:257],Model_zero[::-1])
    
    # plt.rcParams.update({'font.size': 35})
    import compy
    compy.plt_params()
    plt.figure(dpi=300,figsize=(30,25))
    
    plt.subplot(211)

    # for i in range(0,len(Czp)):
    #     plt.semilogx(f,Czp[i],linewidth=0.25,color='r')
        
    for i in range(0,len(High_Czp)):
        plt.semilogx(f,High_Czp[i],linewidth=1,color='g')
      
    plt.title("High-Coherence Teleseismic Events")
    plt.semilogx(f,np.median(High_Czp,axis=0),color='blue',linewidth=3,label= "Median of EQs")
    plt.vlines(x = f_min, ymin=0, ymax=1,color='r',linestyles='dashed',linewidth=3,label="Frequency limits ")
    plt.vlines(x = f_max, ymin=0, ymax=1,color='r',linestyles='dashed',linewidth=3)
    plt.ylabel("Coherence ")
    # plt.xlabel("Frequency [Hz] ")
    plt.xlim([freq1,freq2])
    plt.grid(True)
    plt.legend(loc='lower left')
    # plt.text(0.02, 0.9, 'a)', transform=plt.gca().transAxes, fontsize=40, fontweight='bold')

    plt.subplot(212)
    # for i in range(0,len(ratio_psd)):
    #     plt.semilogx(f,np.sqrt(ratio_psd[i])/(-invz[0][0][0].elevation*rho),linewidth=0.5,color='r')
        
        
    for i in range(0,len(High_ratio_psd)):
        plt.semilogx(f,np.sqrt(High_ratio_psd[i])/(-invz[0][0][0].elevation*rho),linewidth=1,color='g')
        
      
    plt.semilogx(f, Data,linewidth=3,color='blue',label="Measured P/a ratio")
    # plt.semilogx(f,Data,color='r',linewidth=5)
    
    plt.semilogx(f_dispersion_curve,Model,color='purple',label='Theoretical P/a Ratio',linewidth=3)
    
    plt.semilogx(f,Data*gain_factor,color='black',label='Corrected P/a ratio',linewidth=3,linestyle='--')

    plt.vlines(x = f_min, ymin=-1, ymax=10,color='r',linestyles='dashed',label="Frequency limits ",linewidth=3)
    plt.vlines(x = f_max, ymin=-1, ymax=10,color='r',linestyles='dashed',linewidth=3)
    
    plt.xlim([freq1,freq2])
    plt.ylim([-1,10])
    plt.grid(True)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Spectral ratio / ''\u03C1''h')
    plt.title('Spectral Ratio ' + "[Gain factor of "+str(stream[0].stats.station) +" is "+ "%.2f" % gain_factor+"]")
    plt.legend(loc='lower left')

    
    plt.tight_layout()
    # plt.text(0.02, 0.9, 'b)', transform=plt.gca().transAxes, fontsize=40, fontweight='bold')
    # plt.show()
    
    # file_path_save = "/Users/mohammadamin/Desktop/figures_1Feb"
    # plt.savefig(file_path_save + f"Calibration.pdf")

    return(gain_factor)
#%%
def rayleigh_arrival(stream,window = 20 , timelag = - 2,plot_condition = False):
    
    max_index = np.argmax(stream.select(channel="*Z")[0].data)  # Find the maximum value in the trace
    
    fs = stream.select(channel="*Z")[0].stats.sampling_rate
    
    max_time = max_index / fs  # Calculate the time of the maximum value
    
    ray_arr = stream.select(channel="*Z")[0].stats.starttime + max_time + (timelag*60)
    end_time = stream.select(channel="*Z")[0].stats.starttime + max_time + ((timelag+window)*60)

    print(f"Fundamental Rayleigh arrival at {ray_arr} seconds.")
    
    st = stream.copy()
    
    st.trim(starttime=ray_arr,endtime = end_time)
    
    if plot_condition:
        plt.figure(dpi=300,figsize=(12,6))
        plt.subplot(211)
        plt.plot(stream.select(channel="*Z")[0].data)
        plt.vlines(max_index + (timelag*60)*fs, np.min(stream.select(channel="*Z")[0].data), np.max(stream.select(channel="*Z")[0].data),color='r')
        plt.vlines(max_index + ((timelag+window)*60*fs), np.min(stream.select(channel="*Z")[0].data), np.max(stream.select(channel="*Z")[0].data),color='r')
    
        plt.subplot(212)
        plt.plot(st.select(channel="*Z")[0].data)

    return(st,(max_time + (timelag*60)),(max_time + ((timelag+window)*60)))
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
    # misfit = np.sum(((d-m)/s)**l)
    misfit = np.sum(((d-m)**l)/(s**2))

    return(misfit)
#%%
# it can be better by writing the code for step size, you did before somewhere!!!
def grid_search(d,m):
    lower_limit = 0
    upper_limit = 500
    sensitivity = 0.01
    
    grided_d = np.zeros([(abs(lower_limit) + abs(upper_limit)),len(d)])
    misfit_value = np.zeros([(abs(lower_limit) + abs(upper_limit)),1])

    for i in range(lower_limit,upper_limit):
    
        multi_factor = i*sensitivity
        # print(multi_factor)
        grided_d[i] = multi_factor * d
        misfit_value[i] = misfit(grided_d[i],m,l=2,s=1)

    minimum_index = np.argmin(misfit_value)
    
    gain_factor = 1 / (lower_limit + (minimum_index*sensitivity))

    print("Gain factor is " + str(1/gain_factor))
    
    return(1/gain_factor)
#%%
def theoretical_p_a_ratio(alpha = 1500, rho = 1028,h = 4760,plot_condition=True,t = None,velocity_model = None):
    '''
    Carefull it is not angular frequency !!! does it matter???

    Parameters
    ----------
    alpha : TYPE, optional
        DESCRIPTION. The default is 1500.
    rho : TYPE, optional
        DESCRIPTION. The default is 1028.
    h : TYPE, optional
        DESCRIPTION. The default is 4760.
    plot_condition : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.

    '''
    f_dispersion_curve , fun_phase_vel = phase_dispersion(t=t)
    
    r = 2*np.pi * f_dispersion_curve * np.sqrt(alpha**-2 - (fun_phase_vel*1000)**-2)
        
    p_a_ratio = ( rho * np.sin( r * h))/(r * np.cos( r * h))
 
    if plot_condition:
        plt.rcParams.update({'font.size': 20})
        plt.figure(dpi=300,figsize=(12,8))
        plt.title('Normalized Theoretical P/a ratio')
        plt.semilogx(f_dispersion_curve,p_a_ratio/(rho*h),linewidth=3,color='blue')
        plt.xlim([0.01,0.10])
        plt.ylim([-1,20])
        plt.ylabel("Ratio / \u03C1H")
        plt.xlabel('Frequency [Hz]')
        plt.grid(True)
        plt.legend(loc='upper right',fontsize=15)
        plt.tight_layout()
        plt.show()
        
    return(f_dispersion_curve,p_a_ratio/(rho*h))

#%%
def phase_dispersion(velocity_model = None,plot_condition=True ,t = None):
    # Velocity model
    # thickness, Vp, Vs, density
    # km, km/s, km/s, g/cm3
    if velocity_model is None:
       velocity_model = np.array([
           [0.26, 1.75, 0.34, 1.84],
           
           [0.33, 5.00, 2.70, 2.55],
           [0.33, 5.00, 2.70, 2.55],
           
           [0.60, 6.50, 3.70, 2.85],
           [0.60, 6.50, 3.70, 2.85],
           [0.60, 6.50, 3.70, 2.85],
           
           [1, 7.10, 4.05, 3.05],
           [1, 7.10, 4.05, 3.05],
           [1, 7.10, 4.05, 3.05],
           [1, 7.10, 4.05, 3.05],
           [1, 7.10, 4.05, 3.05],
           
           [2, 7.60, 4.35, 3.25],
           [2, 7.60, 4.35, 3.25],
           [2, 7.60, 4.35, 3.25],
           [2, 7.60, 4.35, 3.25],
           [2, 7.60, 4.35, 3.25],

           [10, 8.02, 4.39, 3.38],
           [10, 8.02, 4.39, 3.38],
           [10, 8.02, 4.39, 3.38],
           [10, 8.02, 4.39, 3.38],
           [10, 8.02, 4.39, 3.38],

           [20, 8.1, 4.4, 3.4],
           [20, 8.1, 4.4, 3.4],
           [20, 8.1, 4.4, 3.4],
           [20, 8.1, 4.4, 3.4],
           [20, 8.1, 4.4, 3.4],

           [25, 8.05, 4.43, 3.36],
           [25, 8.05, 4.43, 3.36],

           [50, 8.5, 4.6, 3.43],
           [50, 8.7, 4.7, 3.49],
           [50, 8.8, 4.75, 3.53],

           [100, 9.38, 5.07, 3.78],
           [100, 9.9, 5.37, 3.92],
           
           [200, 10.1, 5.53, 3.98],
           [400, 11.6, 6.45, 4.68],

           ])

    # velocity_model = np.array([[0.7, 5.00, 2.70, 2.55],
    #                             [1.54, 6.50, 3.70, 2.85],
    #                             [4.75, 7.10, 4.05, 3.05],
    #                             [450.0, 7.55, 4.19, 3.10]
    #                               ])
    
    #CRUST2 + PREM
 
    # plotting 
    vs1 = np.zeros([int(np.sum(velocity_model[:, 0]))+1, 1])
    # vs1 = np.zeros([int(np.sum(velocity_model[:, 0])), 1])

    for i in range(0, len(velocity_model)):
        vs1[int(np.sum(velocity_model[:, 0][0:i])):int(
            np.sum(velocity_model[:, 0][0:i+1]))] = velocity_model[:, 2][i]

    vp1 = np.zeros([int(np.sum(velocity_model[:, 0]))+1, 1])
    # vp1 = np.zeros([int(np.sum(velocity_model[:, 0])), 1])

    for i in range(0, len(velocity_model)):
        vp1[int(np.sum(velocity_model[:, 0][0:i])):int(
            np.sum(velocity_model[:, 0][0:i+1]))] = velocity_model[:, 1][i]
   
    density = np.zeros([int(np.sum(velocity_model[:, 0]))+1, 1])
    # vs1 = np.zeros([int(np.sum(velocity_model[:, 0])), 1])

    for i in range(0, len(velocity_model)):
        density[int(np.sum(velocity_model[:, 0][0:i])):int(
            np.sum(velocity_model[:, 0][0:i+1]))] = velocity_model[:, 3][i]

    depth = np.arange(0, -np.sum(velocity_model[:,0]), -1)
    
    # Periods must be sorted starting with low periods
    if t is None:
        t = np.logspace(0, 3.0, 1000)

    # Compute the 3 first Rayleigh- and Love- wave modal dispersion curves
    # Fundamental mode corresponds to mode 0
    pd = PhaseDispersion(*velocity_model.T)
    cpr = [pd(t, mode=i, wave="rayleigh") for i in range(5)]
    # cpl = [pd(t, mode=i, wave="love") for i in range(3)]
    
    if plot_condition:
        plt.rcParams.update({'font.size': 25})
        plt.figure(dpi=300,figsize=(16,16))
        plt.subplot(221)
        plt.plot(vs1, depth, '-', label="VS ", color='red',linewidth=3)
        plt.plot(vp1, depth, ':', label="VP ", color='blue',linewidth=3)
        plt.grid(True)
        plt.legend(loc='lower left')
        plt.title("Velocity Model")
        plt.xlabel("Velocity [Km/s]")
        plt.ylabel("Depth [Km]")
        plt.text(0.02, 0.9, 'a)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
        # plt.ylim([min(depth)+1,0])
        plt.ylim([-500,0])
        
        plt.subplot(222)
        plt.plot(density, depth, '-', label="Density ", color='blue',linewidth=3)
        plt.grid(True)
        plt.legend(loc='lower left')
        plt.title("Density Model")
        plt.xlabel("Density [g/m^3]")
        plt.ylabel("Depth [Km]")
        plt.text(0.02, 0.9, 'b)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')

        # plt.ylim([min(depth)+1,0])
        plt.ylim([-500,0])

        
        plt.subplot(212)
        plt.title('Rayleigh wave dispersion curve')
        plt.semilogx(1/cpr[0][0],cpr[0][1],linewidth=3,color='blue',label='Mode 0 (Fundamental)')
        plt.semilogx(1/cpr[1][0],cpr[1][1],linewidth=3,color='green',label=' Mode 1')
        # plt.semilogx(1/cpr[2][0],cpr[2][1],linewidth=3,color='red',label=' Mode 2')
        # plt.semilogx(1/cpr[3][0],cpr[3][1],linewidth=3,color='orange',label=' Mode 3')
        # plt.semilogx(1/cpr[4][0],cpr[4][1],linewidth=3,color='yellow',label=' Mode 4')
        plt.ylabel("Phase velocity [km/s]")
        plt.xlabel('Frequency [Hz]')
        plt.text(0.02, 0.9, 'c)', transform=plt.gca().transAxes, fontsize=25, fontweight='bold')
        
        plt.grid(True)
        plt.legend(loc='lower left')
        plt.tight_layout()
        plt.show()
        
    return(1/cpr[0][0],cpr[0][1])
 
#%%
def cut_signal_above_zero(signal):
    """
    Cuts a signal from where it has a value above zero and removes trailing zeros.

    Parameters:
    - signal (list): A list of values representing the signal.

    Returns:
    - cut_signal (list): The cut signal without trailing zeros.
    """
    start_index = 0

    for i, value in enumerate(signal):
        if value > 0:
            start_index = i
            break

    cut_signal = signal[start_index:]

    # Remove trailing zeros
    end_index = len(cut_signal)
    for i, value in enumerate(cut_signal[::-1]):
        if value != 0:
            end_index = len(cut_signal) - i
            break

    cut_signal = cut_signal[:end_index]
    return cut_signal
#%%
def calculate_speed_of_sound_in_water(temperature=4, salinity=35, depth=4760):

    """Calculate the speed of sound in water using Mackenzie's formula."""
    # Constants
    a1 = 1448.96
    a2 = 4.591
    a3 = -5.304 * 10**(-2)
    a4 = 2.374 * 10**(-4)
    a5 = 1.340
    a6 = 1.630 * 10**(-2)
    a7 = 1.675 * 10**(-7)
    a8 = -1.025 * 10**(-2)
    a9 = -7.139 * 10**(-13)
    
    t = temperature
    s = salinity
    d = depth

    # Calculate speed of sound
    speed_of_sound = (
        a1 + a2 * t - a3 * t**2 + a4 * t**3 + a5 * (s - 35) + a6 * d + a7 * d**2 -
        a8 * t * (s - 35) - a9 * t * d**3
    )
    return speed_of_sound

#%%
from math import sin, cos, sqrt, atan2, radians

def calculate_distance(lat1, lon1, lat2, lon2):
    # Convert coordinates to radians
    lat1_rad = radians(lat1)
    lon1_rad = radians(lon1)
    lat2_rad = radians(lat2)
    lon2_rad = radians(lon2)
    
    # Earth's radius in kilometers
    radius = 6371.0
    
    # Calculate the differences between the coordinates
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad
    
    # Haversine formula
    a = sin(dlat / 2) ** 2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    distance = radius * c
    
    return distance


#%%
# Calibrating via P-wave Arrival
client = Client("RESIF")
# start =UTCDateTime("2012-10-12")
net = "YV"
sta = "RR38"

def pressure_calibration(stream,mag=7,i=1):
    '''
    Pressure gauge calibration using p wave arrival
    
    '''
    invz = client.get_stations(
        network=stream[0].stats.network,
        station=stream[0].stats.station,
        channel="BHZ",
        location="*",
        level="response")

    invp = client.get_stations(
        network=stream[0].stats.network,
        station=stream[0].stats.station,
        channel="BDH",
        location="*",
        level="response")
    eq_spans = tiskit.TimeSpans.from_eqs(stream.select(channel='*Z')[0].stats.starttime,
                                     stream.select(
                                         channel='*Z')[0].stats.endtime,
                                     minmag=mag, days_per_magnitude=0.5)
    
    stream.trim(eq_spans.start_times[i],eq_spans.start_times[i]+10*3600)
    
    stream.select(channel="*Z").remove_response(inventory=invz,
                                            output="ACC", plot=False)

    stream.select(channel="*H").remove_response(inventory=invp,
                                            output="DEF", plot=False)


    df = stream[0].stats.sampling_rate

    cft = obspy.signal.trigger.recursive_sta_lta(
        stream.select(channel="*Z")[0].data, int(5 * df), int(10 * df))
    
    # trig = coincidence_trigger("recstalta", 7, 2, stream, 1, sta=int(5 * df), lta=int(10 * df))
    
    p_arrival_times = obspy.signal.trigger.trigger_onset(cft,1,0.2) / df   
    plot_trigger(stream.select(channel="*Z")[0], cft,1, 0.5)

    
    stream.trim(stream[0].stats.starttime + p_arrival_times[i-1][0] -2 ,stream[0].stats.starttime + p_arrival_times[i-1][0]+2)
    stream.detrend()
    rho = 1028
    Acc = stream[0].data / (rho*-invz[0][0][0].elevation)
    
    plt.plot(stream[3].data),plt.plot(Acc)
    
    plt.plot(stream[3].data / -Acc),plt.plot(stream[3].data)
    
    scipy.stats.pearsonr(stream[0].data,stream[3].data)
#%%
def p_calibration(stream,gain_factor,rho=1025,mag=6):
    
    nseg=2**11
    TP=5
    
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

    stream.filter('bandpass',freqmin=0.01,freqmax=0.05,corners=4,zerophase=True)
    
    st1 = stream.copy()
    st2 = stream.copy()
    st3 = stream.copy()
    
    stream22 = read()
    stream22.clear()
    
    stream2 = read()
    stream2.clear()
    

    st1.select(channel="*Z").remove_response(inventory=invz,
                                                              output="ACC", plot=False)
    # st2.select(channel="*Z").remove_response(inventory=invz,
    #                                                           output="VEL", plot=False)
    # st3.select(channel="*Z").remove_response(inventory=invz,
    #                                                           output="DISP", plot=False)

    st1.select(channel="*H").remove_response(inventory=invp,
                                                              output="DEF", plot=False)
        
    stream.sort()
    
    eq_spans = tiskit.TimeSpans.from_eqs(stream.select(channel='*Z')[0].stats.starttime,
                                     stream.select(
                                         channel='*Z')[0].stats.endtime,
                                     minmag=mag, days_per_magnitude=0.5)
    
    print ( str(len(eq_spans)) +" Earthquakes with magnitude greater than " + str(mag) +" Mw has been found"  )

    for i in range(0,len(eq_spans)):
        stream1= stream.copy()
        stream1= st1.copy()
        stream22 =stream22 + stream1.trim(eq_spans.start_times[i],eq_spans.start_times[i]+2*3600)
    
    Czp = np.zeros([len(eq_spans),int(nseg/2 + 1)])

   # Calculation coherence
   
    print("Calculating Coherence and Spectral ratio")
    for i in range(0,len(eq_spans)):
        f,Czp[i] = scipy.signal.coherence(stream22[i*4:(i+1)*4].select(component='Z')[0].data,
                                         stream22[i*4:(i+1)*4].select(component='H')[0].data,
                                         fs=stream22[i*4:(i+1)*4][0].stats.sampling_rate,
                                         nperseg =nseg,noverlap=(nseg*0.5),
                                         window=scipy.signal.windows.tukey(nseg,
                                         (TP*60*stream22[i*4:(i+1)*4][0].stats.sampling_rate)/nseg))


    coherence_mask = (f >= 0.01) & (f <= 0.1)
    
    print ("Selecting Earthquakes with median of Coherence greater than "+ str(0.95) +" and mean greater than " + str(0.8) )
   
    High_Czp = [] 
    coh_trsh = 0.95
    mean_trsh = 0.80
    for i in range(0,len(eq_spans)):
                if np.median(Czp[i][coherence_mask]) >  coh_trsh and np.mean(Czp[i][coherence_mask]) >  mean_trsh:
                    High_Czp.append(Czp[i])
                    # print(i)
                    
    plt.rcParams.update({'font.size': 25})
    plt.figure(dpi=300,figsize=(16,12))
    for i in range(0,len(Czp)):
        plt.semilogx(f,Czp[i],linewidth=1.5,color='r')
        
    for i in range(0,len(High_Czp)):
        plt.semilogx(f,High_Czp[i],linewidth=3,color='g')
      
    plt.title("Coherence Of Teleseismic events")
    plt.vlines(x = 0.01, ymin=0, ymax=1,color='black',linestyles='dashed',label="Frequency limits ",linewidth=3)
    plt.vlines(x = 0.1, ymin=0, ymax=1,color='black',linestyles='dashed',linewidth=3)
    plt.ylabel("Coherence ")
    plt.xlabel("Frequency [Hz] ")
    plt.xlim([0.005,0.2])
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.tight_layout()


    plt.figure(dpi=300,figsize=(20,12))
    plt.subplot(211)
    plt.plot(stream22[i*4:(i+1)*4].select(component='Z')[0].data[1500:2000],label="Acceleration Vertical",linewidth=3)
    plt.plot(stream22[i*4:(i+1)*4].select(component='H')[0].data[1500:2000]/ (rho * invz[0][0][0].elevation),
             linestyle='dotted',label="Pressure /ρh",linewidth=3)
    plt.title("station " + sta)
    
    plt.legend(loc="upper right")
    plt.subplot(212)
    plt.plot(stream22[i*4:(i+1)*4].select(component='Z')[0].data[1500:2000],label="Acceleration Vertical",linewidth=3)
    plt.plot(stream22[i*4:(i+1)*4].select(component='H')[0].data[1500:2000]*gain_factor/ (rho * invz[0][0][0].elevation),
             linestyle='dotted',label="Pressure * gain factor /ρh",linewidth=3)
    plt.title("Calibrated, Gain Factor = " + str(gain_factor))

    plt.legend(loc="upper right")
    plt.tight_layout()
    
    plt.figure(dpi=300,figsize=(12,6))    
    plt.plot((st1.select(channel="*Z")[0].data / np.max(st1.select(channel="*Z")[0].data)),label="Acceleration Vertical")
    # plt.plot((st1.select(channel="*H1")[0].data / np.max(st1.select(channel="*H1")[0].data)),label="Acceleration H1")
    # plt.plot((st1.select(channel="*H2")[0].data / np.max(st1.select(channel="*H2")[0].data)),label="Acceleration H2")
    plt.plot((-st1.select(channel="*H")[0].data / np.max(st1.select(channel="*H")[0].data)),label="Pressure")

    plt.legend(loc="upper left")

    
    
    plt.figure(dpi=300,figsize=(12,6))    
    plt.plot((st1.select(channel="*Z")[0].data),label="Acceleration Vertical")
    # plt.plot((st1.select(channel="*H1")[0].data / np.max(st1.select(channel="*H1")[0].data)),label="Acceleration H1")
    # plt.plot((st1.select(channel="*H2")[0].data / np.max(st1.select(channel="*H2")[0].data)),label="Acceleration H2")
    plt.plot(gain_factor*(st1.select(channel="*H")[0].data / (rho * invz[0][0][0].elevation)),label="Pressure")

    plt.legend(loc="upper left")

#%%
import matplotlib as mpl
import scipy.signal

def plot_spectrogram(raw_stream):
    """
    Plots the spectrogram of a raw stream and highlights specific frequency bands.
    
    Parameters:
    - raw_stream: The raw data stream containing seismic or other time series data.
    """
    # Define the number of segments for the spectrogram
    nseg = 2**14

    # Compute the spectrogram
    f, t, Sp = scipy.signal.spectrogram(raw_stream.select(component='H')[0].data,
                                        fs=raw_stream[0].stats.sampling_rate,
                                        nperseg=nseg, noverlap=(nseg/2), window='hann',
                                        scaling="density")

    # Calculate the time difference between consecutive spectrogram points in hours
    # time_diff_hours = (t[1] - t[0]) / 3600

    # Calculate the number of months spanned by the data
    number_month = (raw_stream[0].stats.endtime - raw_stream[0].stats.starttime) // (7*24*3600)

    # Generate date labels for plotting
    dates = []
    for i in range(int(number_month) + 1):
        dates.append(str(raw_stream[0].stats.starttime + 
                         (raw_stream[0].stats.endtime - raw_stream[0].stats.starttime) * i / number_month)[0:10])

    # Determine tick positions for plotting
    tick_positions = [int(i * len(Sp[0]) / (len(dates) - 1)) for i in range(len(dates))]

    # Create a new time array for plotting
    t2 = np.arange(0, len(t))

    # Normalize pressure values for color mapping in the plot
    norm_p = mpl.colors.Normalize(vmin=0, vmax=60)

    # Define frequency bands for infra-gravity and microseismic analyses
    f_min_ig = 0.005
    f_max_ig = 0.02
    f1_ig = np.argmin(np.abs(f - f_min_ig))
    f2_ig = np.argmin(np.abs(f - f_max_ig))

    f_min_ms = 0.1
    f_max_ms = 0.5
    f1_ms = np.argmin(np.abs(f - f_min_ms))
    f2_ms = np.argmin(np.abs(f - f_max_ms))

    # Plot the spectrogram and frequency band analyses
    plt.figure(dpi=300, figsize=(25, 20))
    plt.suptitle(f"{raw_stream[0].stats.network}.{raw_stream[0].stats.station}")
    plt.subplot(121)

    # Spectrogram plot
    plt.pcolormesh(f, t2, 10*np.log10(Sp.T), norm=norm_p)
    plt.vlines(f[f1_ig], 0, t2[-1], linewidth=3, linestyle='dashed', color='blue')
    plt.vlines(f[f2_ig], 0, t2[-1], linewidth=3, linestyle='dashed', color='blue')
    plt.vlines(f[f1_ms], 0, t2[-1], linewidth=3, linestyle='dashed', color='red')
    plt.vlines(f[f2_ms], 0, t2[-1], linewidth=3, linestyle='dashed', color='red')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Date')
    plt.xlim([0.001, 1])
    plt.yticks(tick_positions, dates, rotation=45)
    plt.colorbar(label="Pressure $(Pa^2/Hz)$[dB]", orientation='vertical')
    plt.xscale("log")

    # Frequency band analysis plots
    plt.subplot(122)
    plt.plot(10*np.log10(np.mean(Sp[f1_ig:f2_ig], axis=0)), np.arange(0, len(Sp[0])), 'b', label='Infra-Gravity')
    plt.plot(15*np.log10(np.mean(Sp[f1_ms:f2_ms], axis=0)), np.arange(0, len(Sp[0])), 'red', label="Microsiesmic")
    plt.ylim([0, len(t)])
    plt.xlabel('Pressure $(Pa^2/Hz)$[dB]')
    plt.yticks([])
    plt.xlim([0, 60])

    plt.legend(loc="upper left")
    plt.tight_layout()
    
    f, t, Sp = scipy.signal.spectrogram(raw_stream.select(component='Z')[0].data,
                                        fs=raw_stream[0].stats.sampling_rate,
                                        nperseg=nseg, noverlap=(nseg/2), window='hann',
                                        scaling="density")

    # Calculate the time difference between consecutive spectrogram points in hours
    # time_diff_hours = (t[1] - t[0]) / 3600

    # Calculate the number of months spanned by the data
    number_month = (raw_stream[0].stats.endtime - raw_stream[0].stats.starttime) // (7*24*3600)

    # Generate date labels for plotting
    dates = []
    for i in range(int(number_month) + 1):
        dates.append(str(raw_stream[0].stats.starttime + 
                         (raw_stream[0].stats.endtime - raw_stream[0].stats.starttime) * i / number_month)[0:10])

    # Determine tick positions for plotting
    tick_positions = [int(i * len(Sp[0]) / (len(dates) - 1)) for i in range(len(dates))]

    # Create a new time array for plotting
    t2 = np.arange(0, len(t))

    # Normalize pressure values for color mapping in the plot
    norm_p = mpl.colors.Normalize(vmin=0, vmax=60)

    # Define frequency bands for infra-gravity and microseismic analyses
    f_min_ig = 0.005
    f_max_ig = 0.02
    f1_ig = np.argmin(np.abs(f - f_min_ig))
    f2_ig = np.argmin(np.abs(f - f_max_ig))

    f_min_ms = 0.1
    f_max_ms = 0.5
    f1_ms = np.argmin(np.abs(f - f_min_ms))
    f2_ms = np.argmin(np.abs(f - f_max_ms))

    # Plot the spectrogram and frequency band analyses
    plt.figure(dpi=300, figsize=(25, 20))
    plt.suptitle(f"{raw_stream[0].stats.network}.{raw_stream[0].stats.station}")
    plt.subplot(121)

    # Spectrogram plot
    plt.pcolormesh(f, t2, 10*np.log10(Sp.T))
    plt.vlines(f[f1_ig], 0, t2[-1], linewidth=3, linestyle='dashed', color='blue')
    plt.vlines(f[f2_ig], 0, t2[-1], linewidth=3, linestyle='dashed', color='blue')
    plt.vlines(f[f1_ms], 0, t2[-1], linewidth=3, linestyle='dashed', color='red')
    plt.vlines(f[f2_ms], 0, t2[-1], linewidth=3, linestyle='dashed', color='red')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Date')
    plt.xlim([0.001, 1])
    plt.yticks(tick_positions, dates, rotation=45)
    plt.colorbar(label="Acceleration $(m/s^2)^2$[dB]", orientation='vertical')
    plt.xscale("log")

    # Frequency band analysis plots
    plt.subplot(122)
    plt.plot(10*np.log10(np.mean(Sp[f1_ig:f2_ig], axis=0)), np.arange(0, len(Sp[0])), 'b', label='Infra-Gravity')
    plt.plot(10*np.log10(np.mean(Sp[f1_ms:f2_ms], axis=0)), np.arange(0, len(Sp[0])), 'red', label="Microsiesmic")
    plt.ylim([0, len(t)])
    plt.xlabel('Acceleration $(m/s^2)^2$[dB]')
    plt.yticks([])
    # plt.xlim([0, 60])

    plt.legend(loc="upper left")
    plt.tight_layout()
    
    

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
def coherogram_spectrogram_alpha(st,nseg=2**12,tw =1,f_min = 0.005,f_max = 0.02):
    
    Tresh_coh = 0.8
    Tresh_Dz = -170
    Tresh_Dp = -70
    
    
    plt.set_cmap("jet")
    norm_z = mpl.colors.Normalize(vmin=-190, vmax=-120)
    norm_p = mpl.colors.Normalize(vmin=0, vmax=60)
    norm_coh = mpl.colors.Normalize(vmin=0, vmax=1)
    
    
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
    plt.suptitle(st[0].stats.station)
    # Subplot 1 (twice the width)
    ax1 = plt.subplot(gs[0])
    # Your plotting commands for subplot 1...
    plt.pcolormesh(f, t2, 10*np.log10(Sp.T),norm=norm_p)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Time [Date]')
    plt.xscale('log')
    plt.xlim(0.001,1)    
    plt.title("Spectrogram BDH", y=1.025)
    plt.vlines(f[f1], 0, t2[-1],linewidth=5,linestyle='dashed',color='black')
    plt.vlines(f[f2], 0, t2[-1],linewidth=5,linestyle='dashed',color='black')
    
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
    # plt.xlim([20,40])
    plt.yticks([])
    # Subplot 3 (twice the width, next to subplot 2)
    ax3 = plt.subplot(gs[2])
    # Your plotting commands for subplot 3...
    plt.pcolormesh(f, t2, 10*np.log10(Sz.T*(2*np.pi*f)**4))  
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.xlim(0.001,1)    
    plt.title("Spectrogram BHZ", y=1.025)
    plt.yticks([])
    plt.colorbar(label="Acceleration$((m/s^2)^2/Hz)$ [dB] ",orientation='vertical')
    plt.vlines(f[f1], 0, t2[-1],linewidth=5,linestyle='dashed',color='black')
    plt.vlines(f[f2], 0, t2[-1],linewidth=5,linestyle='dashed',color='black')
    
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
    # plt.xlim([-180,-140])
    
    # Subplot 5 (twice the width, next to subplot 4)
    ax5 = plt.subplot(gs[4])
    # Your plotting commands for subplot 5...
    plt.pcolormesh(f,np.arange(0,len(Czp)),Czp,norm = norm_coh)
    plt.xscale('log')
    plt.xlim(0.001,1)    
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
    # Subplot 6 (half the width of subplot 5)
    ax6 = plt.subplot(gs[5])
    # Your plotting commands for subplot 6...
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
    plt.grid(True)
    plt.ylim([0,len(Czp)])
    # plt.xlim([0.5,1])
    plt.yticks([])
    plt.subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.07, wspace=0.02, hspace=0.2)
    plt.tight_layout(pad=1.0, w_pad=0, h_pad=2.0)

#%%
from scipy.signal import stft

def plot_stft(stream,nperseg=2**14):
  fs = stream[0].stats.sampling_rate
  
  """
  Performs a Short-Time Fourier Transform (STFT) on the provided signal and plots the results.

  Parameters:
  - signal: The input signal (a NumPy array).
  - fs: Sampling frequency of the signal (default is 1000Hz).
  - nperseg: Length of each segment (default is 256).
  - noverlap: Number of points to overlap between segments (default is None, which defaults to nperseg // 2).
  - nfft: Number of points in the FFT (default is None, which defaults to nperseg).
  - figsize: Tuple indicating the size of the plot (default is (10, 6)).
  """
  # Calculate the STFT
  f, t, Zxx = stft(stream[0], fs, nperseg=nperseg)
  norm_p = mpl.colors.Normalize(vmin=0, vmax=1)

  # Plotting
  plt.figure(dpi=300,figsize=(20,15))
  plt.pcolormesh(t, f, np.abs(Zxx),norm=norm_p)
  plt.colorbar(label='Magnitude')
  plt.ylabel('Frequency [Hz]')
  plt.xlabel('Time [sec]')
  plt.yscale("log")
  plt.title('STFT Magnitude')
  plt.ylim(0.001,1)
  plt.show()

#%%
def phase_frequency(stream,n = 2**14):

    fft_result_z = scipy.fft.fft(stream.select(channel="*Z")[0].data,n=n)
    fft_result_p = scipy.fft.fft(stream.select(channel="*H")[0].data,n=n)
    fft_result_h1 = scipy.fft.fft(stream.select(channel="*1")[0].data,n=n)
    fft_result_h2 = scipy.fft.fft(stream.select(channel="*2")[0].data,n=n)
    
    # Phase Calculation: Compute the phase for each frequency component
    phase_z = np.angle(fft_result_z,deg=True)
    phase_p = np.angle(fft_result_p,deg=True)
    phase_h1 = np.angle(fft_result_h1,deg=True)
    phase_h2 = np.angle(fft_result_h2,deg=True)
    
    # Frequency axis (for plotting)
    sampling_rate = stream[0].stats.sampling_rate  # Example: 1000 samples per second
    
    freq = np.linspace(0, sampling_rate/2, n//2)
    
    # Plot the phase spectrum
    plt.figure(dpi=300, figsize=(35, 25))
    # First subplot
    plt.suptitle("Phase Spectrum--"+ str(stream[0].stats.station)+"--"+
                 str(stream[0].stats.starttime)[0:10]+"--"+str(stream[0].stats.endtime)[0:10])
    plt.subplot(2, 2, 1) # (rows, columns, subplot number)
    plt.plot(freq, phase_z[:n//2], "o", color='black', markersize=3) # Smaller markers
    plt.xlim(0.001,1)
    plt.title('Vertical')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Phase [°]')
    plt.xscale('log')
    plt.grid(True)
    
    # Second subplot
    plt.subplot(2, 2, 2)
    plt.plot(freq, phase_h1[:n//2], "o", color='black', markersize=3) # Smaller markers
    plt.xlim(0.001,1)
    plt.title('Horizontal 1')
    plt.xscale('log')
    plt.grid(True)
    
    # Third subplot
    plt.subplot(2, 2, 3)
    plt.plot(freq, phase_p[:n//2], "o", color='black', markersize=3) # Smaller markers
    plt.xlim(0.001,1)
    plt.title('Pressure')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Phase [°]')
    plt.xscale('log')
    plt.grid(True)
    
    # Fourth subplot
    plt.subplot(2, 2, 4)
    plt.plot(freq, phase_h2[:n//2], "o", color='black', markersize=3) # Smaller markers
    plt.xlim(0.001,1)
    plt.title('Horizontal 2')
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.grid(True)
    plt.tight_layout()
    plt.show()


    
    