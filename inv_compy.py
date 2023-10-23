#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 14:04:23 2023

@author: mohammadamin

"""
import random
import numpy as np
import matplotlib.pyplot as plt
def invert_compliace(Data,f,depth_s,starting_model = None,n_layer = 3,sediment_thickness = 80 ,n_sediment_layer=3,sigma_v = 5,sigma_h = 5,iteration = 100000,alpha = 0,sta="RR29"):
# With Liklihood thickness changing + Constrain for thickness
# Accept Probablity = delta s / s ^ 2 (L[i]/L[i+1])
    np.random.seed(0)
    #Data uncertainty 
    # s = np.sqrt(np.var(Data))
    s = np.std(Data) / 2
    # s = uncertainty
    # starting_model,vs0,vsi = model_exp(iteration,first_layer=50,n_layer = 10,power_factor=1.7)
    starting_model,vs0,vsi = Model(iteration, first_layer = 200, n_layer = 13,
                                   power_factor = 1.17,sediment_thickness = sediment_thickness ,n_sediment_layer = n_sediment_layer)
    
    starting_model,vs0,vsi = Model_V2(iteration,n_layer = n_layer)
    
    #starting_model,vs0,vsi = model_crust2(iteration)
    #Constrains
    likelihood_Model = np.zeros([1, iteration])
    likeli_hood = np.zeros([1, iteration])
    likelihood_data = np.zeros([1, iteration])
    mis_fit = np.zeros([1, iteration])
    accept = 0

    ncompl = np.zeros([iteration, np.shape(f)[0]])

            
    for i in range(1, iteration):
        print(iteration - i)
        # r = np.random.uniform(-1, 1)
        r = np.random.randn()
        
        # if mis_fit[0, i] < 3:
        #     # adaptive_step = 1 - (i / iteration)
        #     adaptive_step = 0.75
        # elif mis_fit[0, i] < 2:
        #     # adaptive_step = 1 - (i / iteration)
        #     adaptive_step = 0.5
        # elif mis_fit[0, i] < 2:
        #     # adaptive_step = 1 - (i / iteration)
        #     adaptive_step = 0.25
                
        # else :
        #     adaptive_step = 1
        adaptive_step = 1

        starting_model[:, :, i] = starting_model[:, :, i-1]
        layer = np.random.randint(0, starting_model.shape[0] -1)

        if np.random.randint(0,2) == 0: # Choosing between layer thickness or layer VS
            step_vs = r * sigma_v * adaptive_step  
            
            starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs
            
            #These lines are for keeping inversion in a band
            #it works line model misfit ...
            # if np.sum(starting_model[0:layer,0,i]) < 2000:
            #         if starting_model[layer, 3, i] > 4000:
            #             step_vs = - r * sigma_v * adaptive_step       
            #             starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 
                        
            # if np.sum(starting_model[0:layer,0,i]) < 6000:
            #         if starting_model[layer, 3, i] > 4300:
            #             step_vs = - r * sigma_v * adaptive_step 
            #             starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 

            # if np.sum(starting_model[0:layer,0,i]) > 8000:
            #         if starting_model[layer, 3, i] < 3800:
            #             step_vs = np.abs(r) * sigma_v * adaptive_step 
            #             starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 
                
            # Vp changes based on Vs perturbation (poisson ratio)
            starting_model[layer, 2, i] = velp(starting_model[layer, 3, i])
            
            # Density Changes based on Vs perturbation(Gardner's Eq)
            starting_model[layer, 1, i] = density(starting_model[layer, 3, i])
            
            #updating the S with every iteration
        
        else: 
                
            layer = np.random.randint(0, starting_model.shape[0])
            step_th = r * sigma_h * adaptive_step  
            # layer = np.random.randint(0, 2) # Just sediment layer can change in thickness
            starting_model[layer, 0, i] = starting_model[layer, 0, i-1] + step_th # added to rectified the bug cause by model issue

            if starting_model[layer, 0, i] < 0 :
                step_th =  np.abs(r) * sigma_h * adaptive_step 
                starting_model[layer, 0, i] = starting_model[layer, 0, i-1] + step_th 

            # layer1 = starting_model.shape[0]-1

            starting_model[starting_model.shape[0] -2, 0, i] = starting_model[starting_model.shape[0] -2, 0, i] + step_th
            

        ncompl[i, :] = calc_norm_compliance(depth_s, f, starting_model[:, :, i])

        for ii in range(0, len(starting_model)):
            vsi[0, int(np.sum(starting_model[:, :, i][:, 0][0:ii])):int(np.sum(
                starting_model[:, :, i][:, 0][0:ii+1])), 0] = starting_model[:, 3][ii][i]

        # if i > int(iteration * 0.5):
        #         s_variable = 0.75
        # elif i > int(iteration * 0.75 ):
        #         s_variable = 0.5
        # elif i > int(iteration * 0.90 ):
        #         s_variable = 0.25
        # elif i > int(iteration * 0.95 ):
        #         s_variable = 1 - (i / iteration)
        # else :
        s_variable = 1
            
        
        likelihood_data[0, i] = liklihood(Data, ncompl[i, :], s = s )
    
        # likeli_hood[0, i] = likelihood_data[0, i] + (likelihood_Model[0, i])  # Likihood Data + Model
        
        # likeli_hood[0, i] = liklihood(Data, ncompl[i, :], s=s)  # Liklihood Data
        
        likeli_hood[0, i] = liklihood_roughness(Data, ncompl[i, :],vsi, s = (s * s_variable),alpha= alpha,order=2)    # Liklihood Data + Roughness
        # likeli_hood[0, i] = liklihood_all(Data, ncompl[i, :],vs = vsi,vs_prior=vs0 ,s=s,alpha= 0,beta=0.01,lamda=0,order=2)    # Liklihood Data + Roughness
    
        mis_fit[0, i] = misfit(Data,ncompl[i, :],l=2,s=s)  # Misfit Data
        


        print(mis_fit[0, i])
    
        if likeli_hood[0, i] >= likeli_hood[0, i-1]:
            # starting_model[:, :, i] = starting_model[:, :, i]
            accept+=1

        else:
            # Probablity of Acceptance
            p_candidate = np.random.rand(1)[0]
            # print((likeli_hood[0,i]/likeli_hood[0,i-1]))

            if p_candidate < (likeli_hood[0,i]/likeli_hood[0,i-1]):
                accept+=1
                # starting_model[:, :, i] = starting_model[:, :, i]
            else:                                   #New Line
                  starting_model[:, :, i] = starting_model[:, :, i-1]
                  likeli_hood[0, i] = likeli_hood[0, i-1]
                  mis_fit[0, i] = mis_fit[0, i-1]
                  likelihood_Model[0, i] = likelihood_Model[0, i-1]
                  likelihood_data[0, i] = likelihood_data[0, i-1]
                  ncompl[i, :] = ncompl[i-1, :] 
                  
    vs = np.zeros([iteration, int(np.sum(starting_model[0:-1, 0, 0])), 1])

    for j in range(0, iteration):
        print(iteration - j)
        for i in range(0, len(starting_model)):
            vs[j, int(np.sum(starting_model[:, :, j][:, 0][0:i])):int(np.sum(
                starting_model[:, :, j][:, 0][0:i+1])), 0] = starting_model[:, 3,j][i]
            
    print("Accept Rate Is " + str(accept/iteration) +"%")
    accept_rate = accept/iteration
    # burnin = int(0.2 * iteration)
    # # burnin = 10000
    # plot_inversion(starting_model,vs,mis_fit,ncompl,Data,likelihood_data=likeli_hood,freq=f,sta=sta,
    #                iteration=iteration,s=s,burnin = burnin,mis_fit_trsh = 5)
    
    # plot_hist(starting_model,burnin=burnin,mis_fit=mis_fit,mis_fit_trsh = 1)
    
    # autocorreletion(starting_model,iteration)
    
    return(starting_model,vs,mis_fit,ncompl,likelihood_data,accept_rate)
#%%
def Lcurve(Data,f,depth_s,starting_model = None,sigma_v = 5,sigma_h = 1,iteration = 100000):
    
    burnin = int(0.8 * iteration)
    
    # Define the parameters for the logarithmic array
    start = 0.01  
    stop = 10  # Stop value
    num_points = 10  # Number of points in the array
    base = 10  # Logarithm base

    # Create a logarithmic array using numpy's logspace function
    alpha = np.logspace(np.log10(start), np.log10(stop), num=num_points, base=base)
    
    # alpha = np.linspace(start, stop, num=num_points)
    
    l_curve = []
    for i in range(0,len(alpha)):
        
        starting_model,vs,mis_fit,ncompl,likeli_hood = invert_compliace(Data
                                                                   ,f
                                                                   ,depth_s
                                                                   ,starting_model = None
                                                                   ,sigma_v = sigma_v
                                                                   ,sigma_h = sigma_h
                                                                   ,iteration = iteration
                                                                   ,alpha = alpha[i])
        l_curve.append(mis_fit)
        
    l_curve_burnin = np.zeros([len(l_curve),1])
    
    for i in range(0,len(l_curve)):
        l_curve_burnin[i] = np.median(l_curve[i][0][burnin:iteration])
    
    
    plt.loglog(alpha,l_curve_burnin)
    
    return(l_curve,alpha)

#%%
def model_exp(iteration, first_layer = 200, n_layer = 15,power_factor = 1.15):
    starting_model = np.zeros([n_layer+1, 4, iteration])
    dep = 0
    for i in range(0,len(starting_model)):
        starting_model[i][0][0] = np.float16(first_layer*(power_factor)**i) # Thickness
    
        # starting_model[0][3][0] = 150 # VS First Sediment layer
        # starting_model[1][3][0] = 250 # VS Second Sediment layer
        # starting_model[2][3][0] = 350 # VS Third Sediment layer    
        # starting_model[3][3][0] = 500 # VS Forth Sediment layer    
        
        # starting_model[1][3][0] = 350 # VS
        starting_model[len(starting_model)-1][0][0] = 100000 # Thickness
        starting_model[len(starting_model)-1][1][0] = 3340
        starting_model[len(starting_model)-1][2][0] = 8120
        starting_model[len(starting_model)-1][3][0] =  4510
        
        print("Thickness of layer " + str(i+1) + " is " + str(starting_model[i][0][0]) + " m")
    
        dep = starting_model[i][0][0] + dep
    
        starting_model[i][3][0] = 1200 + ((i/len(starting_model))*3000) # VS
    
        # starting_model[i][3][0] = 2000

        starting_model[i][2][0] = velp(starting_model[i][3][0]) # Vp
    
        starting_model[i][1][0] = density(starting_model[i][3][0]) # Density
        #RR38
    
    #Half space


    
    #     model4 =   np.array([[700, 2550, 5000, 2700]
    #               ,[1540, 2850, 6500, 3700]
    #               ,[4750, 3050, 7100, 4050]
    #               ,[3010, 3100, 7530, 4190]])

    # # starting_model[:,:,0] = DD
    # starting_model[:,:,0] = model4

    # starting_model[0][3][0] = 2700
    # starting_model[1][3][0] = 2700
    # starting_model[2][3][0] = 2700


    # starting_model[3][3][0] = 3700
    # starting_model[4][3][0] = 3700
    # starting_model[5][3][0] = 3700


    # starting_model[6][3][0] = 4050
    # starting_model[7][3][0] = 4050
    # starting_model[8][3][0] = 4050
    # starting_model[9][3][0] = 4050
    # starting_model[10][3][0] = 4050

    # starting_model[11][3][0] = 4370
    # # starting_model[12][3][0] = 4370
    # # starting_model[13][3][0] = 4370
    # # starting_model[24][3][0] = 4190

    # # starting_model[9][0][0] = 10000
    
    vs0 = np.zeros([1, int(np.sum(starting_model[:, 0, 0])), 1])
    for i in range(0, len(starting_model)):
        vs0[0, int(np.sum(starting_model[:, :, 0][:, 0][0:i])):int(np.sum(
            starting_model[:, :, 0][:, 0][0:i+1])), 0] = starting_model[:, 3][i][0]

    vsi = np.zeros([1, int(np.sum(starting_model[:, 0, 0])), 1])
    print("Depth above half-space is "+ str(dep - starting_model[len(starting_model)-1][0][0]) + " m")
    return(starting_model,vs0,vsi)
#%%

def Model(iteration, first_layer = 200, n_layer = 10,power_factor = 1.17,sediment_thickness = 80,n_sediment_layer = 1):
    
    # Calculate the logarithmic series that sums up to the total distance
    log_sum_series = np.geomspace(1, sediment_thickness, n_sediment_layer)
    # Normalize the series to make the sum exactly equal to the total distance
    log_sum_series /= np.sum(log_sum_series)  # Normalize
    # Calculate the distances for each layer
    log_spaced_distances = sediment_thickness * log_sum_series
    
    vs_sediments = np.geomspace(340, 600, n_sediment_layer)
    
    starting_model = np.zeros([n_layer+1+n_sediment_layer, 4, iteration])
    dep = 0
    
      
    #Sediments
    for i in range(0,n_sediment_layer):
        starting_model[i][0][0] = log_spaced_distances[i]# Thickness
  
    
        starting_model[i][3][0] = vs_sediments[i] # VS
    
        # starting_model[i][3][0] = 2000

        starting_model[i][2][0] = velp(starting_model[i][3][0]) # Vp
    
        starting_model[i][1][0] = density(starting_model[i][3][0]) # Density
        print("Thickness of Sediment-layer " + str(i+1) + " is " + str(starting_model[i][0][0]) + " m")
        dep = starting_model[i][0][0] + dep
        
    for i in range(n_sediment_layer,len(starting_model)):
        starting_model[i][0][0] = np.float16(first_layer*(power_factor)**(i - n_sediment_layer + 1)) # Thickness

        
        
        print("Thickness of layer " + str(i+1) + " is " + str(starting_model[i][0][0]) + " m")
    
        dep = starting_model[i][0][0] + dep
    
        starting_model[i][3][0] = 1200 + (((i - n_sediment_layer + 1)/len(starting_model))*3000) # VS
    
        # starting_model[i][3][0] = 2000

        starting_model[i][2][0] = velp(starting_model[i][3][0]) # Vp
    
        starting_model[i][1][0] = density(starting_model[i][3][0]) # Density
    
    #RR29
    starting_model[1+n_sediment_layer][3][0] =  2700# VS
    starting_model[1+n_sediment_layer][2][0] =  5000# Vp
    starting_model[1+n_sediment_layer][1][0] =  2550#density
    
    starting_model[2+n_sediment_layer][3][0] =  2700# VS
    starting_model[2+n_sediment_layer][2][0] =  5000# Vp
    starting_model[2+n_sediment_layer][1][0] =  2550#density
    
    starting_model[3+n_sediment_layer][3][0] =  2700# VS
    starting_model[3+n_sediment_layer][2][0] =  5000# Vp
    starting_model[3+n_sediment_layer][1][0] =  2550#density    
    
    
    
    starting_model[4+n_sediment_layer][3][0] =  3700# VS
    starting_model[4+n_sediment_layer][2][0] =  6500# Vp
    starting_model[4+n_sediment_layer][1][0] =  2850#density
    
    starting_model[5+n_sediment_layer][3][0] =  3700# VS
    starting_model[5+n_sediment_layer][2][0] =  6500# Vp
    starting_model[5+n_sediment_layer][1][0] =  2850#density
    
    starting_model[6+n_sediment_layer][3][0] =  3700# VS
    starting_model[6+n_sediment_layer][2][0] =  6500# Vp
    starting_model[6+n_sediment_layer][1][0] =  2850#density
    

    starting_model[7+n_sediment_layer][3][0] =  4050# VS
    starting_model[7+n_sediment_layer][2][0] =  7100# Vp
    starting_model[7+n_sediment_layer][1][0] =  3050#density
    
    starting_model[8+n_sediment_layer][3][0] =  4050# VS
    starting_model[8+n_sediment_layer][2][0] =  7100# Vp
    starting_model[8+n_sediment_layer][1][0] =  3050#density
    
    starting_model[9+n_sediment_layer][3][0] =  4050# VS
    starting_model[9+n_sediment_layer][2][0] =  7100# Vp
    starting_model[9+n_sediment_layer][1][0] =  3050#density
    
    starting_model[10+n_sediment_layer][3][0] =  4050# VS
    starting_model[10+n_sediment_layer][2][0] =  7100# Vp
    starting_model[10+n_sediment_layer][1][0] =  3050#density
    
    starting_model[11+n_sediment_layer][0][0] =  1638# Thickness
    starting_model[11+n_sediment_layer][3][0] =  4050# VS
    starting_model[11+n_sediment_layer][2][0] =  7100# Vp
    starting_model[11+n_sediment_layer][1][0] =  3050#density
    
    dep = dep + 513
    
    starting_model[12+n_sediment_layer][0][0] =  50000# Thickness
    starting_model[12+n_sediment_layer][3][0] =  4510# VS
    starting_model[12+n_sediment_layer][2][0] =  8120# Vp
    starting_model[12+n_sediment_layer][1][0] =  3340#density
    
    
    #Half-Space   
    starting_model[len(starting_model)-1][0][0] = 100000 # Thickness
    starting_model[len(starting_model)-1][1][0] = 3340
    starting_model[len(starting_model)-1][2][0] = 8120
    starting_model[len(starting_model)-1][3][0] =  4510
        
    vs0 = np.zeros([1, int(np.sum(starting_model[:, 0, 0])), 1])
    for i in range(0, len(starting_model)):
        vs0[0, int(np.sum(starting_model[:, :, 0][:, 0][0:i])):int(np.sum(
            starting_model[:, :, 0][:, 0][0:i+1])), 0] = starting_model[:, 3][i][0]

    vsi = np.zeros([1, int(np.sum(starting_model[:, 0, 0])), 1])
    print("Depth above half-space is "+ str(dep) + " m")
    # print("Depth above half-space is "+ str(dep - starting_model[len(starting_model)-1][0][0]) + " m")
    return(starting_model,vs0,vsi)
#%%

def Model_V2(iteration,n_layer=4):
    dep = 0

    if n_layer == 3:
        starting_model = np.zeros([n_layer+3, 4, iteration])
    
        #RR38
        starting_model[0][3][0] =  340 # VS
        starting_model[0][2][0] =  1500 # Vp
        starting_model[0][1][0] =  1820 # Density
        starting_model[0][0][0] =  10 #Thickness
    
        starting_model[1][3][0] =  2700 # VS
        starting_model[1][2][0] =  5000 # Vp
        starting_model[1][1][0] =  2550 # Density
        starting_model[1][0][0] =  700 #Thickness
        
        starting_model[2][3][0] =  3700 # VS
        starting_model[2][2][0] =  6500 # Vp
        starting_model[2][1][0] =  2850 # Density
        starting_model[2][0][0] =  1540 #Thickness
        
        starting_model[3][3][0] =  4050 # VS
        starting_model[3][2][0] =  7100 # Vp
        starting_model[3][1][0] =  3050 # Density
        starting_model[3][0][0] =  4750 #Thickness
        
        starting_model[4][3][0] =  4190 # VS
        starting_model[4][2][0] =  7530 # Vp
        starting_model[4][1][0] =  3100 # Density
        starting_model[4][0][0] =  10000 #Thickness
        
        starting_model[5][3][0] =  4190 # VS
        starting_model[5][2][0] =  7530 # Vp
        starting_model[5][1][0] =  3100 # Density
        starting_model[5][0][0] =  10000 #Thickness
        
    elif n_layer == 6:
        starting_model = np.zeros([n_layer+3, 4, iteration])
    
        #RR38
        starting_model[0][3][0] =  340 # VS
        starting_model[0][2][0] =  1500 # Vp
        starting_model[0][1][0] =  1820 # Density
        starting_model[0][0][0] =  10 #Thickness
    
        starting_model[1][3][0] =  2700 # VS
        starting_model[1][2][0] =  5000 # Vp
        starting_model[1][1][0] =  2550 # Density
        starting_model[1][0][0] =  350 #Thickness
        
        starting_model[2][3][0] =  2700 # VS
        starting_model[2][2][0] =  5000 # Vp
        starting_model[2][1][0] =  2550 # Density
        starting_model[2][0][0] =  350 #Thickness
        
        starting_model[3][3][0] =  3700 # VS
        starting_model[3][2][0] =  6500 # Vp
        starting_model[3][1][0] =  2850 # Density
        starting_model[3][0][0] =  770 #Thickness
        
        starting_model[4][3][0] =  3700 # VS
        starting_model[4][2][0] =  6500 # Vp
        starting_model[4][1][0] =  2850 # Density
        starting_model[4][0][0] =  770 #Thickness
        
        
        starting_model[5][3][0] =  4050 # VS
        starting_model[5][2][0] =  7100 # Vp
        starting_model[5][1][0] =  3050 # Density
        starting_model[5][0][0] =  2375 #Thickness
        
        starting_model[6][3][0] =  4050 # VS
        starting_model[6][2][0] =  7100 # Vp
        starting_model[6][1][0] =  3050 # Density
        starting_model[6][0][0] =  2375 #Thickness
        
        starting_model[7][3][0] =  4190 # VS
        starting_model[7][2][0] =  7530 # Vp
        starting_model[7][1][0] =  3100 # Density
        starting_model[7][0][0] =  10000 #Thickness    
    
        starting_model[8][3][0] =  4190 # VS
        starting_model[8][2][0] =  7530 # Vp
        starting_model[8][1][0] =  3100 # Density
        starting_model[8][0][0] =  10000 #Thickness   
        
    elif n_layer == 9:
        starting_model = np.zeros([n_layer+3, 4, iteration])
    
        #RR38
        starting_model[0][3][0] =  340 # VS
        starting_model[0][2][0] =  1500 # Vp
        starting_model[0][1][0] =  1820 # Density
        starting_model[0][0][0] =  10 #Thickness
    
        starting_model[1][3][0] =  2700 # VS
        starting_model[1][2][0] =  5000 # Vp
        starting_model[1][1][0] =  2550 # Density
        starting_model[1][0][0] =  233 #Thickness
        
        starting_model[2][3][0] =  2700 # VS
        starting_model[2][2][0] =  5000 # Vp
        starting_model[2][1][0] =  2550 # Density
        starting_model[2][0][0] =  233 #Thickness
        
        starting_model[3][3][0] =  2700 # VS
        starting_model[3][2][0] =  5000 # Vp
        starting_model[3][1][0] =  2550 # Density
        starting_model[3][0][0] =  233 #Thickness
        
        
        starting_model[4][3][0] =  3700 # VS
        starting_model[4][2][0] =  6500 # Vp
        starting_model[4][1][0] =  2850 # Density
        starting_model[4][0][0] =  513 #Thickness
        
        starting_model[5][3][0] =  3700 # VS
        starting_model[5][2][0] =  6500 # Vp
        starting_model[5][1][0] =  2850 # Density
        starting_model[5][0][0] =  513 #Thickness
        
        starting_model[6][3][0] =  3700 # VS
        starting_model[6][2][0] =  6500 # Vp
        starting_model[6][1][0] =  2850 # Density
        starting_model[6][0][0] =  513 #Thickness
        
        
        starting_model[7][3][0] =  4050 # VS
        starting_model[7][2][0] =  7100 # Vp
        starting_model[7][1][0] =  3050 # Density
        starting_model[7][0][0] =  1583 #Thickness
        
        starting_model[8][3][0] =  4050 # VS
        starting_model[8][2][0] =  7100 # Vp
        starting_model[8][1][0] =  3050 # Density
        starting_model[8][0][0] =  1583 #Thickness
        
        starting_model[9][3][0] =  4050 # VS
        starting_model[9][2][0] =  7100 # Vp
        starting_model[9][1][0] =  3050 # Density
        starting_model[9][0][0] =  1583 #Thickness
        
        starting_model[10][3][0] =  4190 # VS
        starting_model[10][2][0] =  7530 # Vp
        starting_model[10][1][0] =  3100 # Density
        starting_model[10][0][0] =  10000 #Thickness    
        
        starting_model[11][3][0] =  4190 # VS
        starting_model[11][2][0] =  7530 # Vp
        starting_model[11][1][0] =  3100 # Density
        starting_model[11][0][0] =  10000 #Thickness    
        
        
    else :
        n_layer = 12
        starting_model = np.zeros([n_layer+3, 4, iteration])
    
        #RR38
        starting_model[0][3][0] =  340 # VS
        starting_model[0][2][0] =  1500 # Vp
        starting_model[0][1][0] =  1820 # Density
        starting_model[0][0][0] =  10 #Thickness
        
    
    
        starting_model[1][3][0] =  2700 # VS
        starting_model[1][2][0] =  5000 # Vp
        starting_model[1][1][0] =  2550 # Density
        starting_model[1][0][0] =  175 #Thickness
        
        starting_model[2][3][0] =  2700 # VS
        starting_model[2][2][0] =  5000 # Vp
        starting_model[2][1][0] =  2550 # Density
        starting_model[2][0][0] =  175 #Thickness
        
        starting_model[3][3][0] =  2700 # VS
        starting_model[3][2][0] =  5000 # Vp
        starting_model[3][1][0] =  2550 # Density
        starting_model[3][0][0] =  175 #Thickness
        
        starting_model[4][3][0] =  2700 # VS
        starting_model[4][2][0] =  5000 # Vp
        starting_model[4][1][0] =  2550 # Density
        starting_model[4][0][0] =  175 #Thickness
        
        
        starting_model[5][3][0] =  3700 # VS
        starting_model[5][2][0] =  6500 # Vp
        starting_model[5][1][0] =  2850 # Density
        starting_model[5][0][0] =  385 #Thickness
        
        starting_model[6][3][0] =  3700 # VS
        starting_model[6][2][0] =  6500 # Vp
        starting_model[6][1][0] =  2850 # Density
        starting_model[6][0][0] =  385 #Thickness
        
        starting_model[7][3][0] =  3700 # VS
        starting_model[7][2][0] =  6500 # Vp
        starting_model[7][1][0] =  2850 # Density
        starting_model[7][0][0] =  385 #Thickness
                
        starting_model[8][3][0] =  3700 # VS
        starting_model[8][2][0] =  6500 # Vp
        starting_model[8][1][0] =  2850 # Density
        starting_model[8][0][0] =  385 #Thickness
        
        
        
        starting_model[9][3][0] =  4050 # VS
        starting_model[9][2][0] =  7100 # Vp
        starting_model[9][1][0] =  3050 # Density
        starting_model[9][0][0] =  1187 #Thickness
        
        starting_model[10][3][0] =  4050 # VS
        starting_model[10][2][0] =  7100 # Vp
        starting_model[10][1][0] =  3050 # Density
        starting_model[10][0][0] =  1187 #Thickness
        
        starting_model[11][3][0] =  4050 # VS
        starting_model[11][2][0] =  7100 # Vp
        starting_model[11][1][0] =  3050 # Density
        starting_model[11][0][0] =  1187 #Thickness
        
        starting_model[12][3][0] =  4050 # VS
        starting_model[12][2][0] =  7100 # Vp
        starting_model[12][1][0] =  3050 # Density
        starting_model[12][0][0] =  1187 #Thickness
        
        
        starting_model[13][3][0] =  4190 # VS
        starting_model[13][2][0] =  7530 # Vp
        starting_model[13][1][0] =  3100 # Density
        starting_model[13][0][0] =  10000 #Thickness         
        
        starting_model[14][3][0] =  4190 # VS
        starting_model[14][2][0] =  7530 # Vp
        starting_model[14][1][0] =  3100 # Density
        starting_model[14][0][0] =  10000 #Thickness    
        
        
    vs0 = np.zeros([1, int(np.sum(starting_model[:, 0, 0])), 1])
    for i in range(0, len(starting_model)):
        vs0[0, int(np.sum(starting_model[:, :, 0][:, 0][0:i])):int(np.sum(
            starting_model[:, :, 0][:, 0][0:i+1])), 0] = starting_model[:, 3][i][0]

    vsi = np.zeros([1, int(np.sum(starting_model[:, 0, 0])), 1])
    
    # print("Depth above half-space is "+ str(dep) + " m")
    # print("Depth above half-space is "+ str(dep - starting_model[len(starting_model)-1][0][0]) + " m")
    return(starting_model,vs0,vsi)
    
    # return(starting_model)
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
     # misfit = np.sqrt(np.sum(((d-m)**l)/(s**l)))
     misfit = np.linalg.norm((d-m)/s,ord = l)
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
     # L =   k * np.exp(-0.5*np.sum((d-m)**2/(s**2)))
     # L =   k * np.exp(-0.5*np.linalg.norm(d-m)**2/(s**2))
     
     L =   k * np.exp(-0.5*np.linalg.norm((d-m)/s,ord = 2))

     return(L)
#%%
def liklihood_all(d,m,vs,vs_prior,k=1,s=1,sm=1,alpha=1,beta=1,lamda=1,order=2):

      R_m = alpha*Roughness(vs.flatten(),order)
      L2_data = np.linalg.norm((d - m)/s,ord=2)
      # L2_Model = beta*np.linalg.norm((vs_prior - vs)/sm)
      
      L2_Model = beta*np.sqrt(np.sum(((vs_prior - vs)**2)/(sm**2)))
      
      Damping = lamda * np.sqrt(np.linalg.norm(m/s)**2)
      L =   k * np.exp ( -0.5 * ( ( L2_data + L2_Model + R_m + Damping) ))

      return(L)
 #%%
def liklihood_roughness(d,m,vs,k=1,s=1,alpha=1,order=2):
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
      # L =   k * np.exp(-0.5*np.sum((d-m)**2/(s**2)))
      R_m = alpha*Roughness(vs.flatten(),order)
      # L =   k * np.exp ( -0.5 * ( ( np.linalg.norm(d - m)**2 / (s**2)) + R_m) )
      
      L =   k * np.exp ( -0.5 * ( ( np.linalg.norm((d - m)/s,ord=2)  + R_m) ))

      return(L)
  
#%%
def Roughness(vs,order):
    R_m = np.sqrt(np.sum(np.gradient(vs,order)))
    return(R_m)
#%%
# stable hyperbolic tangent
def dtanh(x):

    a=np.exp(x*(x<=50))
    one=np.ones(np.shape(x))

    y= (abs(x) > 50) * (abs(x)/x) + (abs(x)<=50)*((a-one/a)/(a+one/a))
    return y
#%%

def gravd(W,h):
    # GRAVD Gravity wave wavenumber determination
    # K = rad/meter
    # W = ang freq (can be a vector)
    # h = water depth (meters)
    
    
    if type(W)!=np.ndarray:
        W=np.array([W])
    G=9.79329
    N=len(W)
    W2=W*W
    kDEEP=W2/G
    kSHAL=W/(np.sqrt(G*h))
    erDEEP=np.ones(np.shape(W)) - G*kDEEP*dtanh(kDEEP*h)/W2
    one = np.ones(np.shape(W))
    d=np.copy(one)
    done=np.zeros(np.shape(W))
    nd=np.where(done==0)
    
    k1=np.copy(kDEEP)
    k2=np.copy(kSHAL)
    e1=np.copy(erDEEP)
    ktemp=np.copy(done)
    e2=np.copy(done)
        
    while True:
        e2[nd] = one[nd] - G*k2[nd] * dtanh(k2[nd]*h)/W2[nd]
        d=e2*e2
        done=d<1e-20
        if done.all():
            K=k2
            break
        
        nd=np.where(done==0)
        ktemp[nd]=k1[nd]-e1[nd]*(k2[nd]-k1[nd])/(e2[nd]-e1[nd])
        k1[nd]=k2[nd]
        k2[nd]=ktemp[nd]
        e1[nd]=e2[nd]
   
    return K

#%%

def argdtray(wd,h):
    
    hh = np.sqrt(abs(h)); #% magnitude of wavenumber/freq
    th=wd*hh;			#% number of waves (or e-foldings) in layer in radians
    if th >= 1.5e-14:
        if h <= 0:   #% propagating wave
            c =  np.cos(th);
            s = -np.sin(th)/hh;
        else:			#% evenescent wave
           	d=np.exp(th);
           	c =  0.5*(d + 1/d);
           	s = -0.5*(d - 1/d)/hh;
    else:
        c = 1;
        s = -wd;
    
    return c,s
#%%
def gravd(W,h):
     # GRAVD Gravity wave wavenumber determination
     # K = rad/meter
     # W = ang freq (can be a vector)
     # h = water depth (meters)
     
     
     if type(W)!=np.ndarray:
         W=np.array([W])
     G=9.79329
     N=len(W)
     W2=W*W
     kDEEP=W2/G
     kSHAL=W/(np.sqrt(G*h))
     erDEEP=np.ones(np.shape(W)) - G*kDEEP*dtanh(kDEEP*h)/W2
     one = np.ones(np.shape(W))
     d=np.copy(one)
     done=np.zeros(np.shape(W))
     nd=np.where(done==0)
     
     k1=np.copy(kDEEP)
     k2=np.copy(kSHAL)
     e1=np.copy(erDEEP)
     ktemp=np.copy(done)
     e2=np.copy(done)
         
     while True:
         e2[nd] = one[nd] - G*k2[nd] * dtanh(k2[nd]*h)/W2[nd]
         d=e2*e2
         done=d<1e-20
         if done.all():
             K=k2
             break
         
         nd=np.where(done==0)
         ktemp[nd]=k1[nd]-e1[nd]*(k2[nd]-k1[nd])/(e2[nd]-e1[nd])
         k1[nd]=k2[nd]
         k2[nd]=ktemp[nd]
         e1[nd]=e2[nd]
    
     return K
#%%
def     raydep(P,om,d,ro,vp2,vs2):
    
    """
        %RAYDEP	Propagator matrix sol'n for P-SV waves, minor vector method
    %   u   =  horizontal velocity AT TOP OF EACH LAYER
    %   v   =  vertical velocity AT TOP OF EACH LAYER
    %   sigzx = horizontal stress AT TOP OF EACH LAYER
    %   sigzz = vertical stress AT TOP OF EACH LAYER
    %   {  Normalized compliance = -k*v/(omega*sigzz)  }
    % [v u sigzz sigzx] = raydep(p,omega,d,ro,vp2,vs2)
    %   p     = slowness (s/m) of surface wave
    %   omega = angular frequency (radians/sec) of surface wave
    %   d     = thicknesses of the model layers (meters?)
    %   rho   = density of the layer (kg/m^3) (= gm/cc * 1000)
    %   vp2   = compressional velocity squared (m/s)^2
    %   vs2   = shear velocity squared (m/s)^2
    
    % W.C. Crawford 1-5-89
    %   x = stress-displacement vector:
    %		( vertical velocity,  horiz velocity, vert stress, horiz stress )
    %	y = minor vector matrix
    """


    mu = ro * vs2;
    n = len(d);
    ist = n-1;
    ysav = 0;
    psq = P*P;
    r2 = 2*mu[ist]*P;
    #% R and S are the "Wavenumbers" of compress and shear waves in botlayer
    #% RoW and SoW are divided by ang freq
    RoW = np.sqrt(psq- 1/vp2[ist]);
    SoW = np.sqrt(psq- 1/vs2[ist]);
    ym = np.zeros((ist+1,5));
    i = ist;
    y=np.zeros((5,))
    x=np.zeros((i+1,4))
    y[3-1] =  RoW;
    y[4-1] = -SoW;
    y[1-1] = (RoW*SoW - psq)/ro[i];
    y[2-1] = r2*y[1-1] + P;
    y[5-1] = ro[i] - r2*(P + y[2-1]);
    ym[i,:] = y;
    #%*****PROPAGATE UP LAYERS*********
    while i > 0:
        i = i-1;
        ha = psq-1/vp2[i];
        ca,sa = argdtray(om*d[i],ha);
        hb = psq - 1/vs2[i];
        cb,sb = argdtray(om*d[i],hb);
        hbs = hb*sb;
        has = ha*sa;
        r1 = 1/ro[i];
        r2 = 2*mu[i]*P;
        b1 = r2*y[1-1] - y[2-1];
        g3 = ( y[5-1] + r2*(y[2-1]-b1) ) * r1;
        g1 = b1 + P*g3;
        g2 = ro[i]*y[1-1] - P*(g1+b1);
        e1 = cb*g2 - hbs*y[3-1];
        e2 = -sb*g2 + cb*y[3-1];
        e3 = cb*y[4-1] + hbs*g3;
        e4 = sb*y[4-1] + cb*g3;
        y[3-1] = ca*e2 - has*e4;
        y[4-1] = sa*e1 + ca*e3;
        g3 = ca*e4 - sa*e2;
        b1 = g1 - P*g3;
        y[1-1] = (ca*e1 + has*e3 + P*(g1+b1))*r1;
        y[2-1] = r2*y[1-1] - b1;
        y[5-1] = ro[i]*g3 - r2*(y[2-1] - b1);
        ym[i,:] = y;
    
    de = y[5-1]/np.sqrt(y[1-1]*y[1-1] + y[2-1]*y[2-1]);
    ynorm = 1/y[3-1];
    y[1-1:4] = np.array([0 , -ynorm,  0,  0]);
    #%*****PROPAGATE BACK DOWN LAYERS*********
    while i <= ist:
        x[i,1-1] = -ym[i,2-1]*y[1-1] - ym[i,3-1]*y[2-1] + ym[i,1-1]*y[4-1];
        x[i,2-1] = -ym[i,4-1]*y[1-1] + ym[i,2-1]*y[2-1] - ym[i,1-1]*y[3-1];
        x[i,3-1] = -ym[i,5-1]*y[2-1] - ym[i,2-1]*y[3-1] - ym[i,4-1]*y[4-1]
        x[i,4-1] =  ym[i,5-1]*y[1-1] - ym[i,3-1]*y[3-1] + ym[i,2-1]*y[4-1];
        ls = i;
        if i >=2-1:
            sum = abs( x[i,1-1] + i*x[i,2-1]);
            pbsq = 1/vs2[i];
            if sum < 1e-4:
                break
                
        ha = psq - 1/vp2[i];
        ca, sa = argdtray(om*d[i],ha);
        hb = psq-1/vs2[i];
        cb,sb = argdtray(om*d[i],hb);
        hbs = hb*sb;
        has = ha*sa;
        r2 = 2*P*mu[i];
        e2 = r2*y[2-1] - y[3-1];
        e3 = ro[i]*y[2-1] - P*e2;
        e4 = r2*y[1-1] - y[4-1];
        e1 = ro[i]*y[1-1] - P*e4;
        e6 = ca*e2 - sa*e1;
        e8 = cb*e4 - sb*e3;
        y[1-1] = (ca*e1 - has*e2+P*e8) / ro[i];
        y[2-1] = (cb*e3 - hbs*e4+P*e6) / ro[i];
        y[3-1] = r2*y[2-1] - e6;
        y[4-1] = r2*y[1-1] - e8;
        i = i+1;
    #
    #if x(1,3) == 0
    #  error('vertical surface stress = 0 in DETRAY');
    #end
    ist = ls;
    v =	 x[:,1-1];
    u =  x[:,2-1];
    zz = x[:,3-1];
    zx = x[:,4-1];
    
    return v, u, zz, zx

#%%
def calc_norm_compliance(depth,freq,model):
    """ calculate normalized compliance for a model and water depth
    
    model=[thick(m)  rho(g/cc) vp(m/s) vs(m/s)]
    """

    thick=model[:,0];
    rho=model[:,1];
    vpsq=model[:,2]*model[:,2];
    vssq=model[:,3]*model[:,3];
    omega = 2*np.pi*freq;
    k = gravd(omega,depth);
    p = k / omega;
    
    ncomp=np.zeros((len(p)))
    
    for i in np.arange((len(p))):
        v, u, sigzz, sigzx = raydep( p[i],omega[i],thick,rho,vpsq,vssq)
        ncomp[i] = -k[i]*v[1-1]/(omega[i]*sigzz[1-1]);
        
        #   u   =  horizontal velocity AT TOP OF EACH LAYER
        #   v   =  vertical velocity AT TOP OF EACH LAYER
        #   sigzx = horizontal stress AT TOP OF EACH LAYER
        #   sigzz = vertical stress AT TOP OF EACH LAYER
        
    return ncomp

#%%
def plot_inversion(starting_model,vs,mis_fit,ncompl,Data,likelihood_data,freq,sta,iteration,s,sigma_v,sigma_h,n_layer,alpha,burnin = 50000,mis_fit_trsh = 1):

    depth = np.arange(0, -int(np.sum(starting_model[0:-1, 0, 0])), -1)
    # var_vs = np.zeros([1,vs.shape[0]])

    # for i in range(0, vs.shape[0]):
    #     var_vs[0,i] = np.var(vs[i])
    layer_number  = -2
    Vs_final = np.zeros(vs[0].shape)
    N = 0
    for i in range(burnin,iteration):
        if mis_fit[0][i] < mis_fit_trsh:
                if starting_model[layer_number][3][i] < 5000:
                    Vs_final = Vs_final + vs[i]
                    N = N + 1
    
    Vs_final = Vs_final/N
    print(N)

    plt.rcParams.update({'font.size': 40})
    nn = int((iteration - burnin)/1000) # I want to see just 1000 points of the data
    # nn = 1
    plt.figure(dpi=300, figsize=(30, 25))
    plt.suptitle("Sigma V = " + str(sigma_v) + ", Sigma h = " + str(sigma_h) + ", Alpha = " + str(alpha) + ",N of Layer =" +str(n_layer))
    plt.subplot(121)

    for i in range(burnin, ncompl.shape[0], nn):
        if mis_fit[0][i] < mis_fit_trsh:
            plt.plot(freq, ncompl[i], color='green', linewidth=0.25)
        
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Normalized Compliance')
        
        # plt.plot(freq, Data, color='black', label='Measured Compliance')
    plt.errorbar(freq, Data, yerr=s, ecolor=('black'), color='black',fmt='none',
                    linewidth= 5, label='Measured Compliance',capsize=15)
        
    # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue', 
    #             label='Median of Brun-in')
    plt.plot(freq, ncompl[1],color='black', label='Start Compliance',linewidth= 5,linestyle='dashed')
    # plt.ylim([10e-14,10e-10])
    plt.grid(True)
    plt.legend(loc='upper left',fontsize=30)
    # plt.ylim([1e-11,6e-11])
    # plt.yscale('log')
    
    plt.subplot(122)
    plt.title(sta)
    for i in range(burnin, vs.shape[0]):
        if mis_fit[0][i] < mis_fit_trsh:
                if starting_model[layer_number][3][i] < 5000:

                    plt.plot(vs[i], depth, color = 'green', linewidth = 0.01)
    plt.plot(vs[0], depth, color='black', label='Start Model',linewidth= 5 ,linestyle='dashed')

    # Add labels, legend, and colorbar as in your code
    plt.xlabel('Shear Velocity [m/s]')
    plt.ylabel('Depth [m Below Seafloor]')
    plt.ylim([-12000, 0])
    plt.grid(True)
    plt.legend(loc='lower left')
    plt.tight_layout()

    # plt.show()  # Don't forget to show the plot

    # plt.rcParams.update({'font.size': 30})
    # plt.figure(dpi=300, figsize=(10, 10))
    # for i in range(burnin, vs.shape[0], nn):
    #     if mis_fit[0][i] < mis_fit_trsh:
    #         plt.plot(vs[i], depth, color='grey', linewidth= 1)
            
    # plt.plot(Vs_final, depth, color='black', label='Final Result',linewidth=3)
    # plt.plot(vs[0], depth, color='green', label='Start Model',linewidth=5,linestyle='dashed')
    
    # plt.grid(True)
    # plt.xlabel('Shear Velocity [m/s]')
    # plt.ylabel('Depth [m]')
    # plt.ylim([-int(np.sum(starting_model[:, 0, 0]))+2000, 0])
    
    # plt.ylim([-10000, 0])
    # plt.xlim([0,4500])
    plt.figure(dpi=300,figsize=(16,24))
    plt.subplot(211)
    # plt.plot(likeli_hood[0, 1:iteration-1],color='Blue',label="Total Likelihood")
    # plt.plot(likelihood_Model[0, 1:iteration-1],color='red',label="Model Likelihood")
    plt.plot(likelihood_data[0, 1:iteration-1],color='green',label="Data Likelihood")
    plt.xscale('log')
    plt.vlines(x=burnin, ymin=0, ymax=1, color='r',
                label='Burn-in Region', linestyles='dashed')
    plt.legend(loc='upper left')
    
    plt.title('Likelihood')
    plt.xlabel('Iteration')
    plt.ylabel('Liklihood')
    plt.grid(True)
    
    plt.subplot(212)
    plt.plot(mis_fit[0,1:iteration-1])
    plt.xscale('log')
    # plt.yscale('log')
    
    plt.vlines(x=burnin, ymin=0, ymax=np.max(mis_fit[0,1:iteration-1]), color='r',
               label='Burn-in Region', linestyles='dashed',linewidth=3)
    
    plt.title('Data Misfit')
    plt.xlabel('Iteration')
    plt.ylabel("$\\chi^2$")
    plt.grid(True)
    plt.minorticks_on()

    plt.tight_layout()

    vs_burnin = np.zeros([iteration, int(np.sum(starting_model[:, 0, 0])), 1])
    plt.legend(loc='lower left')
    plt.tight_layout()
    print(N)
#%%    

def plot_hist(starting_model,burnin,mis_fit,mis_fit_trsh): 
    
    starting_model_opt =[]
    for i in range(burnin,len(starting_model[0][0][:])-1):
        if mis_fit[0][i] < mis_fit_trsh:
            starting_model_opt.append(starting_model[:,:,i])    
    starting_model_opt = np.array(starting_model_opt)
    plt.rcParams.update({'font.size': 40})

    plt.figure(dpi=300, figsize=(35, 70))
    plt.suptitle('Shear Velocity Distribution of Layers',y= 0.99)
    for i in range(0, (len(starting_model_opt[0]) -2)):
        
            plt.subplot((len(starting_model_opt[0]) - 1),4,i+1)

            # plt.hist(starting_model_opt[i-1,3,burnin:-1], 100, density=True, facecolor='k', alpha=1)
            plt.hist(starting_model_opt[:,i,3], bins = 20, density = True, histtype = "barstacked",facecolor='k', alpha=1,stacked = True,label='Layer ' + str(i+1))
            plt.xlim([0,4500])
            plt.legend(loc='upper left')
            # plt.title('Layer ' + str(i+1))
            plt.tight_layout()
            
        
    plt.figure(dpi=300, figsize=(35, 70))
    plt.suptitle('Thickness Distribution of Layers',y= 0.99)
    for i in range(0, (len(starting_model_opt[0]) -2)):
        
            plt.subplot((len(starting_model_opt[0]) - 1),4,i+1)

            # plt.hist(starting_model_opt[i-1,3,burnin:-1], 100, density=True, facecolor='k', alpha=1)
            plt.hist(starting_model_opt[:,i,0], bins = 20, histtype = "barstacked",facecolor='k', alpha=1,stacked = True,label='Layer ' + str(i+1))
            plt.xlim([0,2000])
            plt.legend(loc='upper left')
            # plt.title('Layer ' + str(i+1))
            plt.tight_layout()
    

#%%
def plot_hist2d(starting_model,burnin,mis_fit,mis_fit_trsh): 
    
    starting_model_opt =[]
    for i in range(burnin,len(starting_model[0][0][:])-1):
        if mis_fit[0][i] < mis_fit_trsh:
            starting_model_opt.append(starting_model[:,:,i])    
    starting_model_opt = np.array(starting_model_opt)
    
    plt.rcParams.update({'font.size': 40})
    plt.figure(dpi=300, figsize=(50, 100))
    plt.suptitle('Shear Velocity and Thickness Distribution of Layers',y= 0.99)
    for i in range(0, (len(starting_model_opt[0]) -2)):
        
            plt.subplot((len(starting_model_opt[0]) - 1),4,i+1)

            # plt.hist(starting_model_opt[i-1,3,burnin:-1], 100, density=True, facecolor='k', alpha=1)
            plt.hist2d(starting_model_opt[:,i,0],starting_model_opt[:,i,3],density = True,bins=10, cmap='jet')
            # bins = 50, density = True, histtype = "barstacked",facecolor='k', alpha=1,stacked = True,label='Layer ' + str(i+1))
            plt.xlabel("Thickness[m]")
            plt.ylabel('Vs [m/s]')
            # plt.legend(loc='upper left')
            # plt.title('Layer ' + str(i+1))
            plt.colorbar()
            plt.tight_layout()

#%%
def autocorreletion(starting_model,N):
    
    cc = np.zeros([len(starting_model),len(starting_model[0,3,:])])
    
    for i in range(0,len(cc)-1):
        s = starting_model[i,3,:]
        dd = np.correlate(s-np.mean(s),s-np.mean(s),'full')/np.sum((s-np.mean(s))**2)
        # Auto-correlations.
        cc[i] = dd[N-1:]
        

# Estimate of the effective sample size (Gelman et al., 2013).
    Neff = np.zeros([len(starting_model)])
    for i in range(0,len(cc)):
        for j in range(N-1):
            if (cc[i][j]+cc[i][j+1]>0.0):
                Neff[i]+=cc[i][j]
        
    Neff=N/(1.0+2.0*Neff)
    # for i in range(0,len(cc)):
    #     print('Effective Sample Size (parameter 1): %f' % Neff[i])

# Plot autocorrelation function.
    # plt.figure(dpi = 300,figsize=(40,40))
    # plt.rcParams.update({'font.size': 30})
    # for i in range(0,len(cc)-1):
    #     plt.plot(cc[i][0:N],'k',linewidth=3,label='Neff = : %f' % Neff[i])
    #     plt.xlabel('Iteration',labelpad=15)
    #     plt.xlim([0,N])
    #     plt.title('Layer ' + str(i+1))
    #     plt.legend(loc ="upper right",fontsize = 15)
    #     plt.grid()
    #     plt.tight_layout()
        
        
    plt.figure(dpi = 300,figsize=(30,40))
    plt.rcParams.update({'font.size': 30})
    for i in range(0,len(cc)-1):
        plt.subplot(len(cc)-1,3,i+1)
        plt.plot(cc[i][0:N],'k',linewidth=3,label='Neff = : %f' % Neff[i])
        plt.xlabel('Iteration',labelpad=15)
        plt.xlim([0,N])
        plt.title('Layer ' + str(i+1))
        plt.legend(loc ="upper right",fontsize = 15)
        plt.grid()
        plt.tight_layout()

    plt.figure(dpi = 300,figsize=(30,40))
    plt.rcParams.update({'font.size': 30})
    plt.suptitle("Vs")
    for i in range(0,len(cc)-1):
        plt.subplot(len(cc)-1,3,i+1)
        plt.plot(starting_model[i,3,:],'k',linewidth=2)
        plt.xlabel('Iteration',labelpad=15)
        plt.xlim([0,N])
        plt.title('Layer ' + str(i+1))
        plt.grid()
        plt.tight_layout()
        
    plt.figure(dpi = 300,figsize=(30,40))
    plt.rcParams.update({'font.size': 30})
    plt.suptitle("Thickness")
    for i in range(0,len(cc)-1):
        plt.subplot(len(cc)-1,3,i+1)
        plt.plot(starting_model[i,0,:],'k',linewidth=2)
        plt.xlabel('Iteration',labelpad=15)
        plt.xlim([0,N])
        plt.title('Layer ' + str(i+1))
        plt.grid()
        plt.tight_layout()        
    
    return(Neff)
