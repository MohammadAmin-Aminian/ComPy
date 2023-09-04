#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 17:26:47 2023

@author: mohammadamin
"""


a1 = 0
a2 = 20



plt.figure(dpi=300,figsize=(12,16))
plt.subplot(311)
plt.errorbar(freq[a1:a2], Data[a1:a2], yerr=np.std(Data[a1:a2]), ecolor=('r'), color='black',
                linewidth=3, label='Measured Compliance')

plt.subplot(312)
plt.errorbar(freq[a1:a2], Data[a1:a2], yerr=uncertainty_theory[a1:a2], ecolor=('r'), color='black',
                linewidth=3, label='Measured Compliance')

plt.subplot(313)
plt.errorbar(freq[a1:a2], Data[a1:a2], yerr=uncertainty[a1:a2], ecolor=('r'), color='black',
                linewidth=3, label='Measured Compliance')
plt.tight_layout()
