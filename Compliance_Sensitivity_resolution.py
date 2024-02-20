#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 16:29:03 2023

Compliance Resolution and sensitivity 

@author: Mohammad-Amin
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# import tiskit
# import obstools as obs


import numpy as np
import com_forward as cm
import matplotlib.pyplot as plt
import ffplot as fp

import inv_compy 
freq = np.arange(0.005, 0.02, 0.0005)

model0 = inv_compy.model_exp(1,100,8,2.0)[0]
model1 = inv_compy.model_exp(1,100,8,2.0)[0]
model2 = inv_compy.model_exp(1,100,8,2.0)[0]
model3 = inv_compy.model_exp(1,100,8,2.0)[0]
model4 = inv_compy.model_exp(1,100,8,2.0)[0]
model5 = inv_compy.model_exp(1,100,8,2.0)[0]
model6 = inv_compy.model_exp(1,100,8,2.0)[0]
model7 = inv_compy.model_exp(1,100,8,2.0)[0]

percentage = 1.2
# model2[:, 0][layer] = 5000 # Thicknes

# model2[:, 1][layer] = 100 # rho

# model2[:, 2][layer] = 6800 # Vp
model1[:, 3][0] = percentage * model1[:, 3][0] # Vs
model2[:, 3][1] = percentage * model2[:, 3][1] # Vs
model3[:, 3][2] = percentage * model3[:, 3][2] # Vs
model4[:, 3][3] = percentage * model4[:, 3][3] # Vs
model5[:, 3][4] = percentage * model5[:, 3][4] # Vs
model6[:, 3][5] = percentage * model6[:, 3][5] # Vs
model7[:, 3][6] = percentage * model7[:, 3][6] # Vs

ncomp0 = cm.calc_norm_compliance(4560, freq, model0)
ncomp1 = cm.calc_norm_compliance(4560, freq, model1)
ncomp2 = cm.calc_norm_compliance(4560, freq, model2)
ncomp3 = cm.calc_norm_compliance(4560, freq, model3)
ncomp4 = cm.calc_norm_compliance(4560, freq, model4)
ncomp5 = cm.calc_norm_compliance(4560, freq, model5)
ncomp6 = cm.calc_norm_compliance(4560, freq, model6)
ncomp7 = cm.calc_norm_compliance(4560, freq, model7)

vs0 = np.zeros([int(np.sum(model0[:, 0])), 1])
vs1 = np.zeros([int(np.sum(model1[:, 0])), 1])
vs2 = np.zeros([int(np.sum(model2[:, 0])), 1])
vs3 = np.zeros([int(np.sum(model2[:, 0])), 1])
vs4 = np.zeros([int(np.sum(model2[:, 0])), 1])
vs5 = np.zeros([int(np.sum(model2[:, 0])), 1])
vs6 = np.zeros([int(np.sum(model2[:, 0])), 1])
vs7 = np.zeros([int(np.sum(model2[:, 0])), 1])

for i in range(0, len(model0)):
    vs0[int(np.sum(model0[:, 0][0:i])):int(
        np.sum(model0[:, 0][0:i+1]))] = model0[:, 3][i]
    
for i in range(0, len(model1)):
    vs1[int(np.sum(model1[:, 0][0:i])):int(
        np.sum(model1[:, 0][0:i+1]))] = model1[:, 3][i]
        
for i in range(0, len(model2)):    
    vs2[int(np.sum(model2[:, 0][0:i])):int(
        np.sum(model2[:, 0][0:i+1]))] = model2[:, 3][i]

for i in range(0, len(model3)):    
    vs3[int(np.sum(model3[:, 0][0:i])):int(
        np.sum(model3[:, 0][0:i+1]))] = model3[:, 3][i]

for i in range(0, len(model4)):    
    vs4[int(np.sum(model4[:, 0][0:i])):int(
        np.sum(model4[:, 0][0:i+1]))] = model4[:, 3][i]

for i in range(0, len(model5)):    
    vs5[int(np.sum(model5[:, 0][0:i])):int(
        np.sum(model5[:, 0][0:i+1]))] = model5[:, 3][i]

for i in range(0, len(model6)):    
    vs6[int(np.sum(model6[:, 0][0:i])):int(
        np.sum(model6[:, 0][0:i+1]))] = model6[:, 3][i]

for i in range(0, len(model7)):    
    vs7[int(np.sum(model7[:, 0][0:i])):int(
        np.sum(model7[:, 0][0:i+1]))] = model7[:, 3][i]



# vp1 = np.zeros([int(np.sum(model1[:, 0])), 1])

# for i in range(0, len(model1)):
#     vp1[int(np.sum(model1[:, 0][0:i])):int(
#         np.sum(model1[:, 0][0:i+1]))] = model1[:, 2][i]
    
    
depth0 = np.arange(0, -int(np.sum(model1[:, 0, 0])), -1)
depth1 = np.arange(0, -int(np.sum(model1[:, 0, 0])), -1)
depth2 = np.arange(0, -int(np.sum(model2[:, 0, 0])), -1)
depth3 = np.arange(0, -int(np.sum(model2[:, 0, 0])), -1)
depth4 = np.arange(0, -int(np.sum(model2[:, 0, 0])), -1)
depth5 = np.arange(0, -int(np.sum(model2[:, 0, 0])), -1)
depth6 = np.arange(0, -int(np.sum(model2[:, 0, 0])), -1)
depth7 = np.arange(0, -int(np.sum(model2[:, 0, 0])), -1)

plt.rcParams.update({'font.size': 40})
plt.figure(dpi=300, figsize=(20, 20))
plt.subplot(211)
plt.plot(vs1, depth1, linestyle='dashed',  color='purple',linewidth=6)
plt.plot(vs2, depth2, linestyle='dotted' ,color='red',linewidth=6)
plt.plot(vs3, depth3, linestyle='dashdot',color='orange',linewidth=6)
plt.plot(vs4, depth4, linestyle='solid',color='green',linewidth=6)
plt.plot(vs5, depth5, linestyle='dashed',color='yellow',linewidth=6)
plt.plot(vs6, depth6, linestyle='dotted',color='pink',linewidth=6)
plt.plot(vs7, depth7, linestyle='dashdot',color='gray',linewidth=6)
plt.plot(vs0, depth0, linestyle='solid',  color='blue',linewidth=6)

plt.ylim([-12000,0])
plt.grid(True)
plt.legend(loc='upper left')
plt.title("Compliance Sensitivity")

plt.legend(loc='upper left')
plt.grid(True)
plt.xlabel("Shear Velocity [m/s]")
plt.ylabel("Depth [m]")
# plt.yscale('log')

plt.subplot(212)
plt.plot(freq, ncomp0, linestyle='solid',label="Original", color='blue',linewidth=6)
plt.plot(freq, ncomp1, linestyle='dashed',label="1st", color='purple',linewidth=6)
plt.plot(freq, ncomp2, linestyle='dotted',label="2nd", color='red',linewidth=6)
plt.plot(freq, ncomp3, linestyle='dashdot',label="3nd", color='orange',linewidth=6)
plt.plot(freq, ncomp4, linestyle='solid',label="4rd", color='green',linewidth=6)
plt.plot(freq, ncomp5, linestyle='dashed',label="5th", color='yellow',linewidth=6)
plt.plot(freq, ncomp6, linestyle='dotted',label="6th", color='pink',linewidth=6)
plt.plot(freq, ncomp7, linestyle='dashdot',label="7th", color='gray',linewidth=6)

# plt.plot(freq, ncomp0 - ncomp1, linestyle='dashed',label="1st", color='purple',linewidth=6)
# plt.plot(freq, ncomp0 - ncomp2, linestyle='dotted',label="2nd", color='red',linewidth=6)
# plt.plot(freq, ncomp0 - ncomp3, linestyle='dashdot',label="3nd", color='orange',linewidth=6)
# plt.plot(freq, ncomp0 - ncomp4, linestyle='solid',label="4rd", color='green',linewidth=6)
# plt.plot(freq, ncomp0 - ncomp5, linestyle='dashed',label="5th", color='yellow',linewidth=6)
# plt.plot(freq, ncomp0 - ncomp6, linestyle='dotted',label="6th", color='pink',linewidth=6)
# plt.plot(freq, ncomp0 - ncomp7, linestyle='dashdot',label="7th", color='gray',linewidth=6)

# plt.ylim([2e-11,3e-9])
plt.legend(loc='upper left',fontsize=25,facecolor='lightgreen')

plt.grid(True)
plt.title("Compliance")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Normalized Compliance")
plt.tight_layout()

(ncomp2 - ncomp1) / np.std(ncomp1)
np.mean((ncomp2 - ncomp1) / np.std(ncomp1))


#%%
import numpy as np
import com_forward as cm
import matplotlib.pyplot as plt
import ffplot as fp

import inv_compy 
freq = np.arange(0.007, 0.015, 0.0005)

model0 = inv_compy.model_exp(1,500,8,2.0)[0]
model1 = inv_compy.model_exp(1,500,8,2.0)[0]
model2 = inv_compy.model_exp(1,500,8,2.0)[0]


percentage = 1.05
percentage2 = 0.7
layer1 = 1
layer2 = 2
# model2[:, 0][layer] = 5000 # Thicknes

# model2[:, 1][layer] = 100 # rho

# model2[:, 2][layer] = 6800 # Vp
model1[:, 3][layer1] = 0.8 * model1[:, 3][layer1] # Vs
model2[:, 3][layer1] = percentage2 * model2[:, 3][layer1] # Vs

model1[:, 3][layer2] = 0.8 * model1[:, 3][layer2] # Vs
model2[:, 3][layer2] = percentage2 * model2[:, 3][layer2] # Vs

model2[:, 3][0] = percentage2 * model2[:, 3][0] # Vs
model1[:, 3][0] = percentage2 * model1[:, 3][0] # Vs

model1[:, 3][-6] = 1.05 * model1[:, 3][-6] # Vs
model1[:, 3][-5] = percentage * model1[:, 3][-5] # Vs
model2[:, 3][-5] = 0.9 * model2[:, 3][-5] # Vs

ncomp0 = cm.calc_norm_compliance(4265, freq, model0)
ncomp1 = cm.calc_norm_compliance(4265, freq, model1)
ncomp2 = cm.calc_norm_compliance(4265, freq, model2)


vs0 = np.zeros([int(np.sum(model0[:, 0])), 1])
vs1 = np.zeros([int(np.sum(model1[:, 0])), 1])
vs2 = np.zeros([int(np.sum(model2[:, 0])), 1])


for i in range(0, len(model0)):
    vs0[int(np.sum(model0[:, 0][0:i])):int(
        np.sum(model0[:, 0][0:i+1]))] = model0[:, 3][i]
    
for i in range(0, len(model1)):
    vs1[int(np.sum(model1[:, 0][0:i])):int(
        np.sum(model1[:, 0][0:i+1]))] = model1[:, 3][i]
        
for i in range(0, len(model2)):    
    vs2[int(np.sum(model2[:, 0][0:i])):int(
        np.sum(model2[:, 0][0:i+1]))] = model2[:, 3][i]



# vp1 = np.zeros([int(np.sum(model1[:, 0])), 1])

# for i in range(0, len(model1)):
#     vp1[int(np.sum(model1[:, 0][0:i])):int(
#         np.sum(model1[:, 0][0:i+1]))] = model1[:, 2][i]
    
    
depth0 = np.arange(0, -int(np.sum(model1[:, 0, 0])), -1)
depth1 = np.arange(0, -int(np.sum(model1[:, 0, 0])), -1)
depth2 = np.arange(0, -int(np.sum(model2[:, 0, 0])), -1)
depth3 = np.arange(0, -int(np.sum(model2[:, 0, 0])), -1)


plt.rcParams.update({'font.size': 40})
plt.figure(dpi=300, figsize=(20, 20))
plt.subplot(211)
plt.plot(vs1, depth1, linestyle='dashed',  color='green',linewidth=6)
plt.plot(vs2, depth2, linestyle='dotted' ,color='red',linewidth=6)
plt.plot(vs0, depth0, linestyle='solid' ,color='black',linewidth=6)

plt.ylim([-12000,0])
plt.grid(True)
plt.legend(loc='upper left')
plt.title("Compliance Sensitivity")

plt.legend(loc='upper left')
plt.grid(True)
plt.xlabel("Shear Velocity [m/s]")
plt.ylabel("Depth [m]")
# plt.yscale('log')

plt.subplot(212)
plt.plot(freq, ncomp0, linestyle='solid',label="Original", color='black',linewidth=6)
plt.plot(freq, ncomp1, linestyle='dashed',label="1st", color='green',linewidth=6)
plt.plot(freq, ncomp2, linestyle='dotted',label="2nd", color='red',linewidth=6)


# plt.plot(freq, ncomp0 - ncomp1, linestyle='dashed',label="1st", color='purple',linewidth=6)
# plt.plot(freq, ncomp0 - ncomp2, linestyle='dotted',label="2nd", color='red',linewidth=6)
# plt.plot(freq, ncomp0 - ncomp3, linestyle='dashdot',label="3nd", color='orange',linewidth=6)
# plt.plot(freq, ncomp0 - ncomp4, linestyle='solid',label="4rd", color='green',linewidth=6)
# plt.plot(freq, ncomp0 - ncomp5, linestyle='dashed',label="5th", color='yellow',linewidth=6)
# plt.plot(freq, ncomp0 - ncomp6, linestyle='dotted',label="6th", color='pink',linewidth=6)
# plt.plot(freq, ncomp0 - ncomp7, linestyle='dashdot',label="7th", color='gray',linewidth=6)

# plt.ylim([2e-11,3e-9])
# plt.legend(loc='upper left',fontsize=25,facecolor='lightgreen')

plt.grid(True)
plt.title("Compliance")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Normalized Compliance")
plt.tight_layout()

(ncomp2 - ncomp1) / np.std(ncomp1)
np.mean((ncomp2 - ncomp1) / np.std(ncomp1))


