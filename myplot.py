#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Mohammad Amin Aminian
"""
import obstools as obs
import numpy as np
import warnings
import scipy
import scipy.stats
import math as M
import matplotlib.pyplot as plt
import obspy
import obspy.signal
nhnm = obspy.signal.spectral_estimation.get_nhnm()
nlnm = obspy.signal.spectral_estimation.get_nlnm()

warnings.simplefilter("ignore", np.ComplexWarning)

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

def psd(st,nseg=2**13):
    f,Dz = scipy.signal.welch(st.select(component='Z')[0],fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    f,Dh1 = scipy.signal.welch(st.select(component='1')[0],fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    f,Dh2 = scipy.signal.welch(st.select(component='2')[0],fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    f,Dp = scipy.signal.welch(st.select(component='H')[0],fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))

    plt.figure(dpi=300,figsize=(8,6))
    plt.plot(f,10*np.log10(Dz),label="BHZ")             
    plt.plot(f,10*np.log10(Dh1),label="BH1")             
    plt.plot(f,10*np.log10(Dh2),label="BH2")             
    plt.plot(f,10*np.log10(Dp)-140,label="BDH - 140dB")
    
    plt.xscale('log')
    plt.title(st[0].stats.network+"."+st[0].stats.station+"  "+str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Pressure $(Pa^2/Hz)$[dB]\n or Acceleration$((m/s^2)^2/Hz)$ [dB] ")
    plt.grid(True)
    plt.ylim([-200 ,-80])
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.tight_layout()

    plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
    plt.plot(1/nlnm[0],nlnm[1],'--k')
    plt.legend(loc = "upper right",fontsize=8)
#%%

def psd_all(st,st1,st2,st3,nseg=2**13):
    f,Dz = scipy.signal.welch(st.select(component='Z')[0],fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    f,Dz1 = scipy.signal.welch(st1.select(component='Z')[0],fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    f,Dz2 = scipy.signal.welch(st2.select(component='Z')[0],fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))
    f,Dz3 = scipy.signal.welch(st3.select(component='Z')[0],fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))

    plt.figure(dpi=300,figsize=(8,6))
    plt.plot(f,10*np.log10(Dz),label="Raw")             
    plt.plot(f,10*np.log10(Dz1),label="EQ-Removed")             
    plt.plot(f,10*np.log10(Dz2),label="Rotated")             
    plt.plot(f,10*np.log10(Dz3),label="Compliance Removed")             

    plt.xscale('log')
    plt.title(st[0].stats.network+"."+st[0].stats.station+".BHZ  "+str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Pressure $(Pa^2/Hz)$[dB]\n or Acceleration$((m/s^2)^2/Hz)$ [dB] ")
    plt.grid(True)
    plt.ylim([0 ,80])
    # plt.ylim([-200 ,-80])
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.tight_layout()

    plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
    plt.plot(1/nlnm[0],nlnm[1],'--k')
    plt.legend(loc = "upper right",fontsize=8)
#%%
def psd_h(st,tw=6,nseg=2**11):
    
    ws = int(tw*60*60 * st[0].stats.sampling_rate)
    Z ,nd = sliding_window(st.select(component='Z')[0].data, ws = ws,hann=True)
    f,Dz = scipy.signal.welch(Z,fs=st[0].stats.sampling_rate,nperseg = nseg, noverlap=(nseg/2))

    plt.figure(dpi=300)
    for i in range(0,len(Z)):
        plt.plot(f,10*np.log10(Dz[i]),linewidth = 0.5,color="r")             
    
    plt.xscale('log')
    plt.title(st[0].stats.network+"."+st[0].stats.station+"."+st[3].stats.channel+"   " 
              +str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10] +" PSD for every " + str(tw) + " hours")    
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Acc [dB] ")
    plt.grid(True)
    plt.ylim([-200 ,-80])
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
    plt.plot(1/nlnm[0],nlnm[1],'--k')
    plt.legend(loc = "upper right",fontsize=8)
    
        
    plt.figure(dpi=300)
    # plt.subplot(211)
    plt.imshow(10*np.log10(Dz),aspect = 0.002,extent=[f[0],f[len(f)-1],0,nd*tw])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xscale('log')
    plt.xlim([((st[0].stats.sampling_rate*2)/nseg),1])
    plt.ylabel("Time (Hour)",fontsize=14)
    plt.xlabel('Frequency (Hz)',fontsize=14)
    plt.title("PSD of Z Channel Over Time1",fontsize=14)
    plt.colorbar()
    plt.tight_layout()
  
    # plt.subplot(212)
    # plt.imshow(Dz,aspect = 0.002,extent=[f[0],f[len(f)-1],0,nd*tw])
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    # plt.xscale('log')
    # plt.xlim(0.001,0.03)
    # plt.ylabel("Time (Hour)",fontsize=14)
    # plt.xlabel('Frequency (Hz)',fontsize=14)
    # plt.colorbar()
    
#%%

def coh(st,nseg=2**12,TP=5):
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
    
    f,C12 = scipy.signal.coherence(st.select(component='1')[0].data,
                               st.select(component='2')[0].data,
                               fs=st[0].stats.sampling_rate,
                               nperseg =nseg,noverlap=(nseg*0.9),
                               window=scipy.signal.windows.tukey(nseg,
                               (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,Cz1 = scipy.signal.coherence(st.select(component='Z')[0].data,
                               st.select(component='1')[0].data,
                               fs=st[0].stats.sampling_rate,
                               nperseg =nseg,noverlap=(nseg*0.9),
                               window=scipy.signal.windows.tukey(nseg,
                               (TP*60*st[0].stats.sampling_rate)/nseg))
    
    f,Cz2 = scipy.signal.coherence(st.select(component='Z')[0].data,
                               st.select(component='2')[0].data,
                               fs=st[0].stats.sampling_rate,
                               nperseg =nseg,noverlap=(nseg*0.9),
                               window=scipy.signal.windows.tukey(nseg,
                               (TP*60*st[0].stats.sampling_rate)/nseg))
    
    plt.figure(dpi=200,figsize=((8,8)))
    
    plt.subplot(221)
    plt.semilogx(f,Czp)
    plt.title("Coherence (P & Z)   "+
              st[0].stats.network+"."+
              st[0].stats.station+"  "+
              str(st[0].stats.starttime)[0:10]+"--"+
              str(st[0].stats.endtime)[0:10])    
    plt.ylabel("Coherence ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([(2/nseg),1])
    
    plt.subplot(222)
    plt.semilogx(f,C12)
    plt.title("Coherence (H1 & H2)")    
    plt.ylabel("Coherence ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([(2/nseg),1])
    
    plt.subplot(223)
    plt.semilogx(f,Cz1)
    plt.title("Coherence (H1 & Z) ")    
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Coherence ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([(2/nseg),1])
    
    
    plt.subplot(224)
    plt.semilogx(f,Cz2)
    plt.title("Coherence (H2 & Z)")    
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Coherence ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([(2/nseg),1])
    plt.tight_layout()

#%%  
def coh_h(st,tw=6,nseg=2**12,TP = 5):
    '''
     Coherence Over time
     Coherogram
     overlap = 0.9
     
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
    ws = int(tw*60*60 * st[0].stats.sampling_rate)
    Z ,nd = sliding_window(st.select(component='Z')[0].data, ws = ws,hann=True)
    P ,nd = sliding_window(st.select(component='H')[0].data, ws = ws,hann=True)

    plt.figure(dpi=300,figsize=(8,6))
    plt.subplot(211)
    f,Czp = scipy.signal.coherence(Z,P,fs=st[0].stats.sampling_rate,nperseg =nseg,
                                   noverlap=(nseg*0.9),
                                   window=scipy.signal.windows.tukey(nseg,
                                   (TP*60*st[0].stats.sampling_rate)/nseg))
    
    # loc = np.where((np.mean(Czp,axis=1)) == max(np.mean(Czp,axis=1)))
    
    for i in range(0,len(Czp)):
        plt.semilogx(f,Czp[i],linewidth = 0.25,color="r")
    # plt.semilogx(f,Czp[loc[0][0]],linewidth = 0.75,color="b")
    plt.semilogx(f,np.mean(Czp,axis=0),linewidth = 0.75,color="b")

    plt.title("Coherence (P&Z)"+"  for every " + str(tw) + " hours  "+st[0].stats.network+"."+st[0].stats.station+"  "+
              str(st[0].stats.starttime)[0:10]+"--"+str(st[0].stats.endtime)[0:10])    
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Coherence ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([0.001,1])

    plt.subplot(212)
    plt.imshow(Czp,aspect = 0.002,extent=[f[0],f[len(f)-1],0,nd*tw])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xscale('log')
    plt.xlim(0.001,1)
    plt.ylabel("Time (Hour)",fontsize=14)
    plt.xlabel('Frequency (Hz)',fontsize=14)
    plt.title("Coherence Between P and Z  As Function of Time",fontsize=14)
    plt.colorbar()

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
    
    plt.figure(dpi = 300,figsize=(8,8))
    plt.subplot(4,1,1)
    plt.pcolormesh(t, f, 10*np.log10(Sz))
    plt.ylabel('Frequency [Hz]')
    plt.yscale('log')
    plt.ylim([(st[0].stats.sampling_rate*2/nseg),1])
    plt.title("BHZ")
    # plt.colorbar()

    plt.subplot(4,1,2)
    plt.plot(t,Sz[10])
    plt.xlim([0,t[len(t)-1]])
    plt.ylabel('Frequency [Hz]')
    plt.yscale('log')
    # plt.ylim([1e-18,1e-14])
    plt.grid(True)

    plt.subplot(4,1,3)
    plt.pcolormesh(t, f, 10*np.log10(Sp))
    plt.ylabel('Frequency [Hz]')
    plt.yscale('log')
    plt.ylim([(st[0].stats.sampling_rate*2/nseg),1])
    plt.title("BDH")
    # plt.colorbar()
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    