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
freq = np.arange(0.008, 0.015, 0.0005)
plt.rcParams.update({'font.size': 20})

model1 = inv_compy.model_exp(1,50,50,1.05)[0]
model2 = inv_compy.model_exp(1,50,50,1.05)[0]

layer = 49

model2[:, 0][layer] = 6000*10e3 # Thicknes


model2[:, 1][layer] = 4200 # Vs

model2[:, 2][layer] = 6800 # Vp

model2[:, 3][layer] = 4200 #Density

ncomp1 = cm.calc_norm_compliance(4560, freq, model1)
ncomp2 = cm.calc_norm_compliance(4560, freq, model2)

vs1 = np.zeros([int(np.sum(model1[:, 0])), 1])
vs2 = np.zeros([int(np.sum(model2[:, 0])), 1])

for i in range(0, len(model1)):
    vs1[int(np.sum(model1[:, 0][0:i])):int(
        np.sum(model1[:, 0][0:i+1]))] = model1[:, 3][i]
    vs2[int(np.sum(model2[:, 0][0:i])):int(
        np.sum(model2[:, 0][0:i+1]))] = model2[:, 3][i]


# vp1 = np.zeros([int(np.sum(model1[:, 0])), 1])

# for i in range(0, len(model1)):
#     vp1[int(np.sum(model1[:, 0][0:i])):int(
#         np.sum(model1[:, 0][0:i+1]))] = model1[:, 2][i]
    
    

depth = np.arange(0, -int(np.sum(model1[:, 0, 0])), -1)
depth2 = np.arange(0, -int(np.sum(model2[:, 0, 0])), -1)


plt.figure(dpi=300, figsize=(12, 12))
plt.subplot(211)
plt.plot(vs1, depth, '-', label="VS ", color='blue')
plt.plot(vs2, depth2, ':', label="VS ", color='red')
# plt.plot(vp1, depth, ':', label="VP ", color='blue')
# plt.ylim([-20000,0])
plt.grid(True)
plt.legend(loc='upper left')
plt.title("Earth Model- CRUST LASKE MODEL(changed)")

plt.legend(loc='upper left')
plt.grid(True)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Normalized Compliance")
# plt.yscale('log')

plt.subplot(212)
plt.plot(freq, ncomp1, label="Synthetic", color='blue')
plt.plot(freq, ncomp2, ":",label="Synthetic", color='red')

# plt.ylim([2e-11,3e-9])
plt.legend(loc='upper left')

plt.grid(True)
plt.title("Compliance")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Normalized Compliance")
plt.tight_layout()

(ncomp2 - ncomp1) / np.sqrt(np.cov(ncomp1))
