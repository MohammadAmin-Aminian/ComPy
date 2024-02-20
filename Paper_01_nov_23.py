#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 13:44:13 2023

@author: Mohammad-Amin
"""

from obspy import read
from obspy.clients.fdsn import Client
import os
import numpy as np
import matplotlib.pyplot as plt
# import ffplots as fp

client = Client("RESIF")
net = "YV"
sta = "RR50"

# client = Client("NCEDC")
# net = "BK"
# sta = "MOBB"
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
stream.merge(fill_value='interpolate')
#%% Rotating
import compy 
rotated_stream , azimuth,angle = compy.Rotate(stream,time_window = ((1/60)))
#%% Loading rotated stream
A = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Rotated/')
A.sort()

rotated_stream = read()
rotated_stream.clear()

for i in range(0,len(A)):
    rotated_stream = rotated_stream + read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Rotated/'+A[i])
    # stream = stream+ read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+A[i])
rotated_stream.merge(fill_value='interpolate')
# fp.coh_h(rotated_stream)
#%%
# DPG Calibration
import Pressure_calibration as DPG
gain_factor = DPG.calculate_spectral_ratio(stream,mag = 7 ,f_min=0.03,f_max=0.07)

#%%
#Compliance Funtion, Optimizing Compliance based on the Coherence Function

Good_Com_c,Good_Czp,Good_Com_Stream,f_c,f,uncertainty = compy.Calculate_Compliance_beta(rotated_stream,gain_factor=gain_factor)

Good_Com_c1,Good_Czp1,Good_Com_Stream1 = compy.optimizer(Good_Com_c,Good_Czp,Good_Com_Stream,f=f,alpha=0.97)


#%% Compliance Container
import pickle
import scipy
import compy

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


a1= 0
a2= -1


# Gzz = np.zeros([int(len(stream)/4) ,len(f_coh)])
# Gpp = np.zeros([int(len(stream)/4) ,len(f_coh)])
# for i in range(0,int(len(stream)/4)):
#     F_G , Gzz[i] = scipy.signal.welch(stream[i].select(channel="BHZ"),fs = stream[0][0].stats.sampling_rate,nfft=2*len(f_coh)-1)
#     F_G , Gpp[i] = scipy.signal.welch(stream[i].select(channel="BDH"),fs = stream[0][0].stats.sampling_rate,nfft=2*len(f_coh)-1)

# Gzz = Gzz*((2*np.pi*F_G)**4)

# compliance_high,coh_high,stream = compy.optimizer(compliance_high, 
#                                                   coh_high, 
#                                                   stream,
#                                                   f_coh,
#                                                   alpha= 0.8,
#                                                   beta = 0.6,
#                                                   f_min_com=0.008,
#                                                   f_max_com=0.016)

compliance_high,coh_high,stream = compy.optimizer_rms(compliance_high,coh_high,stream,f,a1,a2,
                                                      percentage = 50,
                                                      alpha=0.8,
                                                      beta=0.8)

print("calculating Median of Coherence And Compliance funtion")
Coherence_all = np.median(coh_high,axis=0)
Coherence_all_smoothed = scipy.signal.savgol_filter(Coherence_all,5,1)

Data = np.median(compliance_high,axis=0)[a1:a2]

treshhold = 0.8

s = np.std(compliance_high,axis=0)
# s = np.std(Data)
# Data = scipy.signal.savgol_filter(Data,3,1)
f = f[a1:a2]
s = s[a1:a2]

plt.rcParams.update({'font.size': 40})

plt.figure(dpi=300,figsize=(20,20))
plt.subplot(211)
plt.title(sta)
plt.plot(f,Data,'blue',linewidth=3,label='Compliance')
plt.errorbar(f, Data, yerr=s, fmt='o', ecolor='Black', capsize=15)
plt.ylabel('Compliance')
plt.legend(loc='upper left')
plt.grid(True)

plt.subplot(212)
plt.plot(f_coh,Coherence_all)
plt.hlines(y = treshhold, xmin = f[0] , xmax = f[-1],linestyles='dashed',linewidth=3,color = 'r')
plt.xlim(f[0],f[-1])
plt.xlabel('Frequency [Hz]')
plt.ylabel('Coherence')
plt.grid(True)
plt.tight_layout()
#%% Inversion
import compy
compy.plt_params()
import inv_compy
from inv_compy import invert_compliace
import ffplot as fp
fp.psd
sigma_v = 50
sigma_h = 25
alpha = 0
sediment_thickness = 10
n_sediment_layer = 3

alpha_pressure = 1

iteration = 100000
n_layer= 6

starting_model,vs,vs0,mis_fit,ncompl,likeli_hood,accept_rate = inv_compy.invert_compliace_beta(alpha_pressure * Data,
                                                                          f,
                                                                          depth_s =-inv[0][0][0].elevation,
                                                                          starting_model = None,
                                                                          n_layer= n_layer,
                                                                          sediment_thickness = sediment_thickness,
                                                                          n_sediment_layer = n_sediment_layer,
                                                                          sigma_v = sigma_v,
                                                                          sigma_h = sigma_h,
                                                                          iteration = iteration,
                                                                          s = s / 5 ,
                                                                          alpha = alpha,
                                                                          sta=sta)
d = mis_fit.copy()
d.sort()
d

mis_fit_trsh = 14

burnin =int( iteration * 0)

# inv_compy.plot_inversion_beta(starting_model,
#                             vs,
#                             vs0,
#                             mis_fit,
#                             ncompl,
#                             Data,
#                             likelihood_data=likeli_hood,
#                             freq=f,
#                             sta=sta,
#                             iteration=iteration,
#                             s=s,
#                             sigma_v = sigma_v,
#                             sigma_h = sigma_h,
#                             n_layer = starting_model.shape[0],
#                             alpha = alpha, 
#                             burnin = burnin,
#                             mis_fit_trsh = mis_fit_trsh)

compy.plt_params()
inv_compy.plot_inversion_density(vs,vs0,mis_fit,Data,s,freq=f,sta=sta,burnin=burnin,ncompl=(ncompl),iteration=iteration,mis_fit_trsh = mis_fit_trsh)


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
sigma_v = np.linspace(5,200,25)
Results = Parallel(n_jobs=num_cores)(delayed(invert_compliace)(Data,
                                                               f,
                                                               depth_s =4550,
                                                               starting_model = None,
                                                               n_layer = n_layer ,
                                                               sediment_thickness = 0.1,
                                                               n_sediment_layer = 1,
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



mis_fit_trsh = 4

burnin =int( iteration * 0.1)
for i in range(0,len(Results)):
    
    inv_compy.plot_inversion(Results[i][0],
                             Results[i][1],
                             Results[i][2],
                             Results[i][3],
                             Data,
                             likelihood_data=Results[i][4],
                             freq=f,
                             sta=sta,
                             iteration=iteration,
                             s=np.std(Data),
                             sigma_v = sigma_v[i],
                             sigma_h = 0,
                             n_layer = Results[i][0].shape[0],
                             alpha = 0, 
                             burnin = burnin,
                             mis_fit_trsh = mis_fit_trsh)



#%%
#H step Size Thickness
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
    'stream stats': stream.select(channel="*Z")[0].stats,
    'angle': angle,
    'azimuth': azimuth,
    'gain factor': gain_factor,
}

# Serialize and save the data container using pickle
file_path = "/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Compliance/'
os.makedirs(file_path)

with open(file_path+f'rotation_hourly_{sta}.pkl', 'wb') as f1:
    pickle.dump(rotation_container, f1)
    
print("Data container saved.")
#%%
# Loading Rotation Container 
import pickle
with open(file_path+f'rotation_hourly_{sta}.pkl', 'rb') as file:
    loaded_rotation_container = pickle.load(file)
    
angle = loaded_rotation_container["angle"]

azimuth = loaded_rotation_container["azimuth"]
stream_stats =loaded_rotation_container['stream stats']

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
# plt.title("YV." + str(stream[0].stats.station) +'  '+ str(stream[0].stats.starttime)[0:10]+'--'+str(stream[0].stats.endtime)[0:10]+" Tilt ")
# plt.plot(azimuth, '.', color='black')
# for i in range(0,len(eq_spans)):
#     plt.plot(int((eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)),
#               azimuth[int(eq_spans.start_times[i] - stream[0].stats.starttime) // (time_window*3600)],'o', color='red', markersize=10)

plt.xticks([])
plt.ylabel("Azimuth ["u"\u00b0]")
plt.legend(loc="lower right", facecolor='lightgreen')
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
plt.legend(loc="lower right", facecolor='lightgreen')
    
    # Set x-axis names with rotation and align ticks with the dates generated
plt.xticks(tick_positions, dates, rotation=65)
plt.ylim([-5, 5])
plt.tight_layout()

ii = 1

azimuth_reshaped = azimuth.reshape(24*ii,-1)
angle_reshaped = angle.reshape(24*ii,-1)

plt.figure(dpi=300, figsize=(30, 20))

plt.subplot(211)

for i in range(0,int(len(azimuth)//24*ii)):
    plt.plot(azimuth[i*24*ii:(i+1)*24*ii], '.', color='black')
plt.plot(np.mean(azimuth_reshaped,axis=1),linewidth=3)
         
plt.xticks([])
plt.ylabel("Azimuth ["u"\u00b0]")
plt.legend(loc="lower right", facecolor='lightgreen')
plt.grid(True)

plt.subplot(212)

for i in range(0,int(len(angle)//24*ii)):
    plt.plot(angle[i*24*ii:(i+1)*24*ii], '.', color='black')

plt.plot(np.mean(angle_reshaped,axis=1),linewidth=3)
plt.xlabel("Time [Hour]")
plt.ylabel("Incident Angle ["u"\u00b0]")
plt.grid(True)
plt.legend(loc="lower right", facecolor='lightgreen')
    
    # Set x-axis names with rotation and align ticks with the dates generated
# plt.xticks(tick_positions, dates, rotation=65)
plt.ylim([-2, 2])
plt.tight_layout()

#%% Saving compliance Container 
import pickle
com_container = {
    'Compliance': Good_Com_c,
    'Frequency of Compliance': f_c,
    'Uncertainty':uncertainty,
    'Coherence': Good_Czp,
    'Frequency of Coherence': f,
    'Stream': Good_Com_Stream
}

# Serialize and save the data container using pickle
file_path = "/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Compliance/'
os.makedirs(file_path)

with open(file_path+'com_container.pkl', 'wb') as f1:
    pickle.dump(com_container, f1)
    
print("Data container saved.")



#%% Saving Inversion Container
Inversion_container = {
    'Station': sta,
    'compliance Measured': Data,
    'compliance Frequency': f,
    'uncertainty': s,
    'compliance Forward': ncompl,
    'Models': starting_model,
    'Shear Velocity': vs,
    'Shear Velocity Starting': vs0,
    'Misfit Fucntion': mis_fit,
    'Liklihood': likeli_hood,
    'sigma_v': sigma_v,
    'sigma_h' : sigma_h,
    'alpha' : alpha,
    'alpha_pressure' : alpha_pressure,
    'iteration': iteration,
    'n_layer' : n_layer,
    'burnin' : burnin,
    'mis_fit_trsh' : mis_fit_trsh,
    'accept_rate' : accept_rate
    }

# Serialize and save the data container using pickle
file_path = "/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Compliance/'
# os.makedirs(file_path)

with open(file_path+'Inversion_container.pkl', 'wb') as f1:
    pickle.dump(Inversion_container, f1)
    
print("Data container saved.")

#%%
# Loading Inversion Container
import inv_compy

print("Loading Compliance Container...")
with open('/Users/mohammadamin/Desktop/Data/YV/'+sta+'/Compliance/Inversion_container.pkl', 'rb') as file:
    Inversion_container = pickle.load(file)

sta = Inversion_container["Station"]
sta = Inversion_container["Station"]
Data = Inversion_container["compliance Measured"]
f = Inversion_container["compliance Frequency"]
s = Inversion_container["uncertainty"]
ncompl = Inversion_container["compliance Forward"]
starting_model = Inversion_container["Models"]
vs = Inversion_container["Shear Velocity"]
vs0 = Inversion_container["Shear Velocity Starting"]
mis_fit = Inversion_container["Misfit Fucntion"]
likeli_hood = Inversion_container["Liklihood"]
sigma_v = Inversion_container["sigma_v"]
sigma_h = Inversion_container["sigma_h"]
alpha = Inversion_container["alpha"]
alpha_pressure = Inversion_container["alpha_pressure"]
iteration = Inversion_container["iteration"]
n_layer = Inversion_container["n_layer"]
burnin = Inversion_container["burnin"]
mis_fit_trsh = Inversion_container["mis_fit_trsh"]
accept_rate = Inversion_container["accept_rate"]


compy.plt_params()
inv_compy.plot_inversion_density(vs,vs0,mis_fit,Data,s,freq=f,sta=sta,burnin=burnin,ncompl=(ncompl),iteration=iteration,mis_fit_trsh = mis_fit_trsh)


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





