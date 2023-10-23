#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 13:44:13 2023

@author: Mohammad-Amin
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
# import obspy
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
sta = "RR38"

inv = client.get_stations(
    network=net,
    station=sta,
    channel="BHZ",
    location="*",
    level="response")

A = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/')
#A = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/')

A.sort()

stream = read()
stream.clear()


for i in range(0,len(A)):
    stream = stream + read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/'+A[i])
    #stream = stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])
stream.merge(fill_value='interpolate')

#%% Reading rotated stream
A = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Rotated/')
A.sort()

rotated_stream = read()
rotated_stream.clear()

for i in range(0,len(A)):
    rotated_stream = rotated_stream + read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Rotated/'+A[i])
    # stream = stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])
rotated_stream.merge(fill_value='interpolate')

#%%
# DPG Calibration
import Pressure_calibration as DPG
gain_factor = DPG.calculate_spectral_ratio(stream,mag= 7 ,f_min=0.03,f_max=0.05)

#%%
#Compliance Funtion, Rotation's angles 
import compy 

rotated_stream,azimuth,angle = compy.Rotate(stream,time_window = 4)

High_Com_c,High_Czp,High_Com_Stream,f_c,f,uncertainty,uncertainty_theory = compy.Calculate_Compliance(rotated_stream,gain_factor=gain_factor)

High_Com1,High_Czp1,High_Com_Stream1 = compy.optimizer(High_Com_c,High_Czp,High_Com_Stream,f=f,alpha=0.97)
#%%
import pickle
import scipy
import compy
import numpy as np
import matplotlib.pyplot as plt 

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

Gzz = np.zeros([int(len(stream)/4) ,len(f_coh)])
Gpp = np.zeros([int(len(stream)/4) ,len(f_coh)])
for i in range(0,int(len(stream)/4)):
    F_G , Gzz[i] = scipy.signal.welch(stream[i].select(channel="BHZ"),fs = stream[0][0].stats.sampling_rate,nfft=2*len(f_coh)-1)
    F_G , Gpp[i] = scipy.signal.welch(stream[i].select(channel="BDH"),fs = stream[0][0].stats.sampling_rate,nfft=2*len(f_coh)-1)


# compliance_high,coh_high,stream = compy.optimizer(compliance_high, coh_high, stream, f_coh,alpha= 0.85,f_min_com=0.008,f_max_com=0.016)

print("calculating Median of Coherence And Compliance funtion")
Coherence_all = np.median(coh_high,axis=0)
Coherence_all_smoothed = scipy.signal.savgol_filter(Coherence_all,5,1)

Data = np.median(compliance_high,axis=0)

a1= 2
a2= -18
plt.rcParams.update({'font.size': 35})
plt.figure(dpi=300,figsize=(20,12))

# plt.subplot(211)
# plt.semilogx(f_coh,Coherence_all,'b')
plt.title("Coherence P/Z'' , YV.RR38")
plt.semilogx(f_coh,Coherence_all_smoothed,'b',linewidth = 3)
plt.vlines(x = f[a1], ymin=0, ymax=1,color='black',linestyles="dashed",linewidth=3)
plt.vlines(x = f[a2], ymin=0, ymax=1,color='black',linestyles="dashed",linewidth=3,label="Compliance Frequency Band")
plt.hlines(y = 0.9, xmin = f[a1] , xmax = f[a2],linestyles='dashed',linewidth=3,color = 'r',label="0.9 Treshhold")
plt.plot(F_G,np.mean(Gzz,axis=0)/max(np.mean(Gzz,axis=0)),label = "Normalized PSD (Z'')",linewidth=3)
plt.plot(F_G,np.mean(Gpp,axis=0)/max(np.mean(Gpp,axis=0)),label = "Normalized PSD (P)",linewidth=3)


plt.ylim([0,1])
plt.xlim([0.004,1])
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
plt.legend(loc='upper right',fontsize = 20)
plt.tight_layout()

# plt.subplot(212)
# plt.title(sta)
# plt.plot(f,Data)
# plt.vlines(x = f[a1], ymin=0, ymax=1e-10,color='black',linestyles="dashed")
# plt.vlines(x = f[a2], ymin=0, ymax=1e-10,color='black',linestyles="dashed")
# plt.ylabel('Compliance')
# plt.grid(True)



Data1 = np.median(compliance_high[:,a1:a2],axis=0)
Data = Data1
Data = scipy.signal.savgol_filter(Data1,10,5)
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
plt.hlines(y = 0.9, xmin = f[0] , xmax = f[-1],linestyles='dashed',linewidth=3,color = 'r')
plt.xlim(f[0],f[-1])
plt.xlabel('Frequency [Hz]')
plt.ylabel('Coherence')
plt.grid(True)
plt.tight_layout()
#%% Inversion
import inv_compy
from inv_compy import invert_compliace

sigma_v = 25
sigma_h = 15
alpha = 0.25
sediment_thickness = 10
n_sediment_layer = 6

iteration = 500000
n_layer= 6
starting_model,vs,mis_fit,ncompl,likeli_hood,accept_rate = inv_compy.invert_compliace(Data,
                                                                          f,
                                                                          depth_s =4560,
                                                                          starting_model = None,
                                                                          n_layer= n_layer,
                                                                          sediment_thickness = sediment_thickness,
                                                                          n_sediment_layer = n_sediment_layer,
                                                                          sigma_v = sigma_v,
                                                                          sigma_h = sigma_h,
                                                                          iteration = iteration,
                                                                          alpha = alpha)
burnin = 50
mis_fit_trsh = 3

inv_compy.plot_inversion(starting_model,
                         vs,mis_fit,
                         ncompl,
                         Data,
                         likelihood_data=likeli_hood,
                         freq=f,
                         sta=sta,
                         iteration=iteration,
                         s=np.std(Data)  ,
                         sigma_v = sigma_v,
                         sigma_h = sigma_h,
                         n_layer = starting_model.shape[0],
                         alpha = alpha, 
                         burnin = burnin,
                         mis_fit_trsh = mis_fit_trsh)

inv_compy.plot_hist(starting_model,burnin=burnin,mis_fit=mis_fit,mis_fit_trsh = mis_fit_trsh)

inv_compy.plot_hist2d(starting_model,burnin=burnin,mis_fit=mis_fit,mis_fit_trsh = mis_fit_trsh)

inv_compy.autocorreletion(starting_model, iteration)

#%%
# Alpha Parameterazation
import multiprocessing as mp
from joblib import Parallel , delayed
num_cores = mp.cpu_count()
from inv_compy import invert_compliace

alpha = np.linspace(0,1,20)
# alpha = np.logspace(np.log10(0.1),np.log10(0.4),10)
iteration = 10000
Results = Parallel(n_jobs=num_cores)(delayed(invert_compliace)(Data,
                                                               f,
                                                               depth_s =4560,
                                                               starting_model = None,
                                                               n_layer= n_layer,
                                                               sediment_thickness = sediment_thickness,
                                                               n_sediment_layer = n_sediment_layer,
                                                               sigma_v = 10,
                                                               sigma_h = 10,
                                                               iteration = iteration,
                                                               alpha=alpha[i])for i in range(0,len(alpha)))
MisFit = np.zeros([len(alpha)-1,len(Results[0][2][0])])
burnin = 1000
alpha_best = 10
min_misfit = np.zeros([len(alpha)-1])
max_misfit = np.zeros([len(alpha)-1])

for i in range(0,len(alpha)-1):
    MisFit[i] = Results[i][2][0]
    min_misfit[i] = np.min(MisFit[i][burnin:-1])
    max_misfit[i] = np.max(MisFit[i][burnin:-1])


plt.rcParams.update({'font.size': 25})
plt.figure(dpi=300,figsize=[20,15])
plt.subplot(211)
misfit_mean = np.zeros([len(alpha)])
for i, array in enumerate(MisFit):
    misfit_mean[i] = np.median(MisFit[i][burnin:-1])   
    plt.semilogx(MisFit[i][1:-1],linewidth=2,color='r')
plt.semilogx(MisFit[alpha_best][1:-1],linewidth=3,color='blue',label = "{:.2f}".format(alpha[alpha_best]))
plt.grid(True)
plt.vlines(burnin, np.min(MisFit), np.max(MisFit),linestyles="dashed",linewidth=5,colors='black',label="Burn-in")
plt.legend(loc='upper left')
plt.xlabel("Iteration")
plt.ylabel("$\\chi^2$")
# plt.xlim([1000,iteration])
plt.title("Missfit Trend for Different Alpha Values")

plt.subplot(212)
plt.plot(alpha[0:-1],misfit_mean[0:-1],linewidth=3,label="Median ",color = 'blue')

plt.plot(alpha[0:-1],misfit_mean[0:-1],linewidth=3,label="Best Alpha ", marker='o',color = 'blue', linestyle='-', markersize=18, 
           markevery=[alpha_best], markerfacecolor='red', markeredgecolor='black')

plt.plot(alpha[0:-1],min_misfit,'green',linewidth=3,label = "Minimum Value ")
# plt.plot(alpha[0:-1],max_misfit,'red',linewidth=3,label = "Maximum Value")
plt.xlabel("Alpha")
plt.ylabel("$\\chi^2$")
plt.title("Variation in Exploration Phase with Different Alpha Values")
plt.yscale('log')
plt.legend(loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.show()


#%%
# Number of Layers Parameterazation
import multiprocessing as mp
from joblib import Parallel , delayed
num_cores = mp.cpu_count()
from inv_compy import invert_compliace

n_layer = np.linspace(3,12,4)
n_layer = np.array([3,6,9,12])

# alpha = np.logspace(np.log10(0.1),np.log10(0.4),10)
iteration = 100000
Results = Parallel(n_jobs=num_cores)(delayed(invert_compliace)(Data,
                                                               f,
                                                               depth_s =4560,
                                                               starting_model = None,
                                                               n_layer= n_layer[i],
                                                               sediment_thickness = sediment_thickness,
                                                               n_sediment_layer = n_sediment_layer,
                                                               sigma_v = 25,
                                                               sigma_h = 15,
                                                               iteration = iteration,
                                                               alpha=0)for i in range(0,len(n_layer)))
MisFit = np.zeros([len(n_layer),len(Results[0][2][0])])
burnin = 1000
best_n_layer = 1
min_misfit = np.zeros([len(n_layer)])
max_misfit = np.zeros([len(n_layer)])

for i in range(0,len(n_layer)):
    MisFit[i] = Results[i][2][0]
    min_misfit[i] = np.min(MisFit[i][burnin:-1])
    max_misfit[i] = np.max(MisFit[i][burnin:-1])


plt.rcParams.update({'font.size': 25})
plt.figure(dpi=300,figsize=[20,15])
plt.subplot(211)
misfit_mean = np.zeros([len(n_layer)])
for i, array in enumerate(MisFit):
    misfit_mean[i] = np.median(MisFit[i][burnin:-1])   
    plt.semilogx(MisFit[i][1:-1],linewidth=2,color='r')
plt.semilogx(MisFit[best_n_layer][1:-1],linewidth=3,color='blue',label = "{:.2f}".format(n_layer[best_n_layer]))
plt.grid(True)
plt.vlines(burnin, np.min(MisFit), np.max(MisFit),linestyles="dashed",linewidth=5,colors='black',label="Burn-in")
plt.legend(loc='upper left')
plt.xlabel("Iteration")
plt.ylabel("$\\chi^2$")
# plt.xlim([1000,iteration])
plt.title("Missfit Trend for Different Models")

plt.subplot(212)
# plt.plot(n_layer,misfit_mean,linewidth=3,label="Median ",color = 'blue')

# plt.plot(n_layer,misfit_mean,linewidth=3,label="Best Number of Layers ", marker='o',color = 'blue', linestyle='-', markersize=18, 
#            markevery=[alpha_best], markerfacecolor='red', markeredgecolor='black')

plt.plot(n_layer,min_misfit,'green',linewidth=3,label = "Minimum Value ")
# plt.plot(n_layer,max_misfit,'red',linewidth=3,label = "Maximum Value ")
plt.xlabel("Number of Layers")
plt.ylabel("$\\chi^2$")
plt.title("Variation of $\\chi^2$ for Different Models[Various Number of Layers]")
plt.yscale('log')
plt.legend(loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.show()

#%%
# Vs Step Size
n_layer = 6
sigma_v = np.linspace(5,250,25)
Results = Parallel(n_jobs=num_cores)(delayed(invert_compliace)(Data,
                                                               f,
                                                               depth_s =4550,
                                                               starting_model = None,
                                                               n_layer = n_layer ,
                                                               sediment_thickness = sediment_thickness,
                                                               n_sediment_layer = 3,
                                                               sigma_v = sigma_v[i],
                                                               sigma_h = 0,
                                                               iteration = iteration,
                                                               alpha = 0)for i in range(0,len(sigma_v)))

MisFit = np.zeros([len(sigma_v)-1,len(Results[0][2][0])])
burnin = 1000
sigma_v_best = 9
min_misfit = np.zeros([len(sigma_v)-1])
max_misfit = np.zeros([len(sigma_v)-1])
accepet_rate = np.zeros([len(sigma_v)-1])

for i in range(0,len(sigma_v)-1):
    MisFit[i] = Results[i][2][0]
    accepet_rate[i] = Results[i][5]
    min_misfit[i] = np.min(MisFit[i][burnin:-1])
    max_misfit[i] = np.max(MisFit[i][burnin:-1])


plt.rcParams.update({'font.size': 35})
plt.figure(dpi=300,figsize=[25,30])
plt.subplot(211)
misfit_mean = np.zeros([len(sigma_v)])
for i, array in enumerate(MisFit):
    misfit_mean[i] = np.median(MisFit[i][burnin:-1])   
    plt.semilogx(MisFit[i][1:-1],linewidth=2,color='r')
plt.semilogx(MisFit[sigma_v_best][1:-1],linewidth=3,color='blue',label = "{:.2f}".format(sigma_v[sigma_v_best]))
plt.grid(True)
plt.vlines(burnin, np.min(MisFit), np.max(MisFit),linestyles="dashed",linewidth=5,colors='black',label="Burn-in")
plt.legend(loc='upper left')
plt.xlabel("Iteration")
plt.ylabel("$\\chi^2$")
# plt.xlim([1000,iteration])
plt.title("Missfit Trend for Different Simga Vs Values")

plt.subplot(212)
plt.plot(sigma_v[0:-1],misfit_mean[0:-1],linewidth=3,label="Median ",color = 'blue')

plt.plot(sigma_v[0:-1],misfit_mean[0:-1],linewidth=3,label="Best Simga Vs ", marker='o',color = 'blue', linestyle='-', markersize=18, 
           markevery=[sigma_v_best], markerfacecolor='red', markeredgecolor='black')

# plt.plot(accepet_rate,misfit_mean[0:-1],linewidth=3,label="Median ",color = 'blue')

# plt.plot(sigma_v[0:-1],accepet_rate*10,'k',linewidth=3,label = "Aceept Rate ")
# 
for x,y in zip(sigma_v,5*accepet_rate):

    label = "{:.2f}".format(y/5)

    plt.annotate(label, # this is the text
                 (x,y), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(0,-10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.plot(sigma_v[0:-1],min_misfit,'green',linewidth=3,label = "Minimum Value ")
plt.plot(sigma_v[0:-1],max_misfit,'red',linewidth=3,label = "Maximum Value")
plt.yscale('log')
plt.xlabel("sigma_v")
plt.ylabel("$\\chi^2$")
plt.title("Variation in Exploration Phase with Different Simga Vs Values")
plt.legend(loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.show()


#%%
n_layer = 3
iteration = 10000
sigma_h = np.linspace(5,100,20)
Results = Parallel(n_jobs=num_cores)(delayed(invert_compliace)(Data,
                                                               f,
                                                               depth_s =4550,
                                                               starting_model = None,
                                                               n_layer = n_layer,
                                                               sediment_thickness = sediment_thickness,
                                                               n_sediment_layer = 3,
                                                               sigma_v = 50,
                                                               sigma_h = sigma_h[i],
                                                               iteration = iteration,
                                                               alpha = 0.25)for i in range(0,len(sigma_h)))

MisFit = np.zeros([len(sigma_h)-1,len(Results[0][2][0])])
burnin = 1000
sigma_h_best = 2
min_misfit = np.zeros([len(sigma_h)-1])
max_misfit = np.zeros([len(sigma_h)-1])

for i in range(0,len(sigma_h)-1):
    MisFit[i] = Results[i][2][0]
    min_misfit[i] = np.min(MisFit[i][burnin:-1])
    max_misfit[i] = np.max(MisFit[i][burnin:-1])


plt.rcParams.update({'font.size': 25})
plt.figure(dpi=300,figsize=[20,15])
plt.subplot(211)
misfit_mean = np.zeros([len(sigma_h)])
for i, array in enumerate(MisFit):
    misfit_mean[i] = np.median(MisFit[i][burnin:-1])   
    plt.semilogx(MisFit[i][1:-1],linewidth=2,color='r')
plt.semilogx(MisFit[sigma_h_best][1:-1],linewidth=3,color='blue',label = "{:.2f}".format(sigma_h[sigma_h_best]))
plt.grid(True)
plt.vlines(burnin, np.min(MisFit), np.max(MisFit),linestyles="dashed",linewidth=5,colors='black',label="Burn-in")
plt.legend(loc='upper left')
plt.xlabel("Iteration")
plt.ylabel("$\\chi^2$")
# plt.xlim([1000,iteration])
plt.title("Missfit Trend for Different Sigma h Values")

plt.subplot(212)
plt.loglog(sigma_h[0:-1],misfit_mean[0:-1],linewidth=3,label="Median ",color = 'blue')

plt.loglog(sigma_h[0:-1],misfit_mean[0:-1],linewidth=3,label="Best Sigma h ", marker='o',color = 'blue', linestyle='-', markersize=18, 
           markevery=[sigma_h_best], markerfacecolor='red', markeredgecolor='black')

plt.loglog(sigma_h[0:-1],min_misfit,'green',linewidth=3,label = "Minimum Value ")
plt.loglog(sigma_h[0:-1],max_misfit,'red',linewidth=3,label = "Maximum Value")
plt.xlabel("Sigma h")
plt.ylabel("$\\chi^2$")
plt.title("Variation in Exploration Phase with Different Sigma h Values")
plt.legend(loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.show()
# Smoothing Curve

#%%
# Saving Rotated Stream

output_dir = "/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Rotated/'

os.makedirs(output_dir)

compy.split_and_save_stream(rotated_stream,30*24*60,output_dir)


#%% Saving Rotation Container  
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
#%%
# Loading Rotation Container 
import pickle
with open('/Users/mohammadamin/Desktop/Data/YV/'+sta+'/Compliance/rotation_container.pkl', 'rb') as file:
    loaded_rotation_container = pickle.load(file)
    
ang = loaded_rotation_container["angle"]

azi = loaded_rotation_container["azimuth"]

plt.rcParams.update({'font.size': 25})
plt.figure(dpi=300,figsize=[25,15])
plt.subplot(211)
plt.plot(ang,linewidth = 3,color="blue")
plt.ylabel('Angle')
plt.xlabel('Time')

plt.subplot(212)
plt.plot(azi,linewidth = 3,color='blue')
plt.ylabel('Azimuth')
plt.xlabel('Time')

#%% Saving compliance Container 
import pickle
com_container = {
    'Compliance': High_Com_c,
    'Frequency of Compliance': f_c,
    'Uncertainty':uncertainty,
    'Uncertainty Theory':uncertainty_theory,
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



#%% Saving Inversion Container
Inversion_container = {
    'Models': starting_model,
    'Shear Velocity': vs,
    'Misfit Fucntion': mis_fit,
    'Compliance Forward': ncompl,
    'compliance Measured': Data,
    'Liklihood': likeli_hood,
    'Compliance Frequency': f
    }

# Serialize and save the data container using pickle
file_path = "/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Compliance/'
os.makedirs(file_path)

with open(file_path+'Inversion_container.pkl', 'wb') as f1:
    pickle.dump(Inversion_container, f1)
    
print("Data container saved.")

#%%
# Load compliance Container
print("Loading Compliance Container...")
with open('/Users/mohammadamin/Desktop/Data/YV/'+sta+'/Compliance/Inversion_container.pkl', 'rb') as file:
    loaded = pickle.load(file)

# compliance_high = loaded["Compliance"]
# f = loaded["Frequency of Compliance"]
# coh_high = loaded["Coherence"]
# f_coh = loaded["Frequency of Coherence"]
# stream = loaded["Stream"]


# print("calculating Median of Coherence And Compliance funtion")


# Coherence_all = np.median(coh_high,axis=0)
# Data = np.median(compliance_high,axis=0)

# plt.subplot(211)
# plt.plot(f,Data)
# plt.subplot(212)
# plt.semilogx(f_coh,Coherence_all)


















