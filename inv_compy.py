#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 14:04:23 2023

@author: mohammadamin

"""
import random
import numpy as np
import matplotlib.pyplot as plt
def invert_compliace(Data,f,depth_s,starting_model = None,s=None,n_layer = 3,sediment_thickness = 80 ,n_sediment_layer=3,sigma_v = 5,sigma_h = 5,iteration = 100000,alpha = 0.25,sta="RR38"):
# With Liklihood thickness changing + Constrain for thickness
# Accept Probablity = delta s / s ^ 2 (L[i]/L[i+1])
    borders = np.array([0.5,0.1])
    
    np.random.seed(0)
    #Data uncertainty 
    if s == None:
        s = np.sqrt(np.var(Data))/20
    # s = uncertainty
    # starting_model,vs0,vsi = model_exp(iteration,first_layer=50,n_layer = 10,power_factor=1.7)
    starting_model,vs0,vsi = Model(iteration, first_layer = 200, n_layer = 13,
                                   power_factor = 1.17,sediment_thickness = sediment_thickness ,n_sediment_layer = n_sediment_layer)
    
    starting_model,vs0,vsi = Model_V2(iteration,n_layer = n_layer,sta=sta)
    
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


        if mis_fit[0, i]/(np.sqrt(np.var(Data)) / s) < 1:
            # adaptive_step = 1 - (i / iteration)
            adaptive_step = 0.5 
        else :
            adaptive_step = 1
        adaptive_step = 1

        starting_model[:, :, i] = starting_model[:, :, i-1]
        layer = np.random.randint(0, starting_model.shape[0] -1)

        if np.random.randint(0,2) == 0: # Choosing between layer thickness or layer VS
            step_vs = r * sigma_v * adaptive_step  
            
            starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs
            
            
            
            if starting_model[layer, 3, i] < starting_model[layer, 3, 0] * (1 - borders[0]) or starting_model[layer, 3, i] > starting_model[layer, 3, 0] * (1 + borders[1]) :
                step_vs = - r * sigma_v * adaptive_step       
                starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 
                           
            
            # These lines are for keeping inversion in a band
            # it works like model misfit ...
            # if 0 < np.sum(starting_model[0:layer,0,i]) < 500:
            #         if starting_model[layer, 3, i] < 300 or starting_model[layer, 3, i] > 3000 :
            #             step_vs = - r * sigma_v * adaptive_step       
            #             starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 
                        
                        
            # if 500 < np.sum(starting_model[0:layer,0,i]) < 2200:
            #         if starting_model[layer, 3, i] < 3300 or starting_model[layer, 3, i] > 4070:
            #             step_vs = - r * sigma_v * adaptive_step       
            #             starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 
                                                
            # if 2200 < np.sum(starting_model[0:layer,0,i]) < 7000:
            #         if starting_model[layer, 3, i] < 3600 or starting_model[layer, 3, i] > 4400:
            #             step_vs = - r * sigma_v * adaptive_step       
            #             starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 
                        
            # if 7000 < np.sum(starting_model[0:layer,0,i]) < 17500:
            #         if starting_model[layer, 3, i] < 4000 or starting_model[layer, 3, i] > 4400:
            #             step_vs = - r * sigma_v * adaptive_step       
            #             starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 

    
                        
            # if 0 < np.sum(starting_model[0:layer,0,i]) < 710:
            #         while starting_model[layer, 3, i] < 300 or starting_model[layer, 3, i] > 3000 :
            #             step_vs = np.random.randn() * sigma_v * adaptive_step       
            #             starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 
            #             print(1)
            #             print(starting_model[layer, 3, i])
                        
            # if 710 < np.sum(starting_model[0:layer,0,i]) < 2250:
            #         while starting_model[layer, 3, i] < 3300 or starting_model[layer, 3, i] > 4070:
            #             step_vs = np.random.randn() * sigma_v * adaptive_step       
            #             starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 
            #             print(2)
            #             print(starting_model[layer, 3, i])
                                                
            # if 2250 < np.sum(starting_model[0:layer,0,i]) < 7000:
            #         while starting_model[layer, 3, i] < 3600 or starting_model[layer, 3, i] > 4400:
            #             step_vs = np.random.randn() * sigma_v * adaptive_step       
            #             starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 
            #             print(3)
            #             print(starting_model[layer, 3, i])
                      
            # if 7000 < np.sum(starting_model[0:layer,0,i]) < 15000:
            #         while starting_model[layer, 3, i] < 4000 or starting_model[layer, 3, i] > 4400:
            #             step_vs = np.random.randn() * sigma_v * adaptive_step       
            #             starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs
            #             print(4)
            #             print(starting_model[layer, 3, i])
                        

                                         
                        
                                
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
        


        print(mis_fit[0, i] / (np.sqrt(np.var(Data)) / s))
    
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
def invert_compliace_beta(Data,f,depth_s,starting_model = None,s=None,n_layer = 3,sediment_thickness = 80 ,n_sediment_layer=3,sigma_v = 5,sigma_h = 5,iteration = 100000,alpha = 0.25,sta="RR38"):
# With Liklihood thickness changing + Constrain for thickness
# Accept Probablity = delta s / s ^ 2 (L[i]/L[i+1])
    borders = np.array([0.85,0.20])
    borders_deep = np.array([0.1,0.1])
    
    # borders = np.array([5,5])
    # borders_deep = np.array([5,5])
    
    np.random.seed(0)
    #Data uncertainty 
    # if s == None:
    #     s = np.sqrt(np.var(Data))/20
    # s = uncertainty
    # starting_model,vs0,vsi = model_exp(iteration,first_layer=50,n_layer = 10,power_factor=1.7)
    # starting_model,vs0,vsi = Model(iteration, first_layer = 200, n_layer = 13,
    #                                power_factor = 1.17,sediment_thickness = sediment_thickness ,n_sediment_layer = n_sediment_layer)
    
    starting_model,vs0,vsi = Model_V2(iteration,n_layer = n_layer,sta=sta)
    
    # starting_model,vs0,vsi = model_crust2(iteration)
    
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

        adaptive_step = 1

        starting_model[:, :, i] = starting_model[:, :, i-1]
        layer = np.random.randint(0, starting_model.shape[0] -1)

        if np.random.randint(0,2) == 0: # Choosing between layer thickness or layer VS
            step_vs = r * sigma_v * adaptive_step  
            
            starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs
            
            
            if layer < int(starting_model.shape[0]-2) :    
                
                # if starting_model[layer, 3, i] < starting_model[layer, 3, 0] * (1 - borders[0]) or starting_model[layer, 3, i] > starting_model[layer, 3, 0] * (1 + borders[1]) :
                #     step_vs = - r * sigma_v * adaptive_step       
                #     starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 
                
                while starting_model[layer, 3, i] < starting_model[layer, 3, 0] * (1 - borders[0]) or starting_model[layer, 3, i] > starting_model[layer, 3, 0] * (1 + borders[1]) :
                    r = np.random.randn()
                    step_vs =  r * sigma_v * adaptive_step   
                    
                    print("Shallow Borders")
                    
                    starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 


            if layer == int(starting_model.shape[0]-2) :    
                if starting_model[layer, 3, i] < starting_model[layer, 3, 0] * (1 - borders_deep[0]) or starting_model[layer, 3, i] > starting_model[layer, 3, 0] * (1 + borders_deep[1]) :
                    step_vs = - r * sigma_v * adaptive_step       
                    starting_model[layer, 3, i] = starting_model[layer, 3, i-1] + step_vs 
                    print("Deep Borders")
                        
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

        s_variable = 1
        
        likelihood_data[0, i] = liklihood(Data, ncompl[i, :], s = s )

        likeli_hood[0, i] = liklihood_roughness(Data, ncompl[i, :],vsi, s = (s * s_variable),alpha= alpha,order=2)    # Liklihood Data + Roughness
    
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
    
    return(starting_model,vs,vs0,mis_fit,ncompl,likelihood_data,accept_rate)

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
    
        starting_model[i][3][0] = 3000 + ((i/len(starting_model))*3000) # VS
    
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
def Model_V2(iteration,n_layer=3,sta = "RR38"):
    dep = 0
    if sta == "RR28":
        if n_layer == 3:
            starting_model = np.zeros([n_layer+3, 4, iteration])
        # Thickness, Density, Vp, Vs 
            starting_model[:,:,0]= np.array([[   260.,  1820.,  1750.,   340.],
                                             [  700.,  2550.,  5000.,  2700.],
                                             [ 1540.,  2850.,  6500.,  3700.],
                                             [ 4750.,  3050.,  7100.,  4050.],
                                             [10000.,  3100.,  7530.,  4190.],
                                             [10000.,  3100.,  7530.,  4190.]])
       
        
        elif n_layer == 6:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            
            starting_model[:,:,0] = np.array([[   260.,  1820.,  1750.,   340.],
                                              [  350.,  2550.,  5000.,  2700.],
                                              [  350.,  2550.,  5000.,  2700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.]])

        
        elif n_layer == 9:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[   260.,  1820.,  1750.,   340.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.]])
        
        else :
            n_layer = 12
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[   260.,  1820.,  1750.,   340.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.]])

    elif sta == "RR29":
        if n_layer == 3:
            starting_model = np.zeros([n_layer+3, 4, iteration])
        
            starting_model[:,:,0]= np.array([[   10.,  1820.,  1750.,   340.],
                                             [  700.,  2550.,  5000.,  2700.],
                                             [ 1540.,  2850.,  6500.,  3700.],
                                             [ 4750.,  3050.,  7100.,  4050.],
                                             [10000.,  3100.,  7530.,  4190.],
                                             [10000.,  3100.,  7530.,  4190.]])
       
        
        elif n_layer == 6:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            
            starting_model[:,:,0] = np.array([[   10.,  1820.,  1750.,   340.],
                                              [  350.,  2550.,  5000.,  2700.],
                                              [  350.,  2550.,  5000.,  2700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.]])

        
        elif n_layer == 9:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[   10.,  1820.,  1750.,   340.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.]])
        
        else :
            n_layer = 12
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[   10.,  1820.,  1750.,   340.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.]])
            
    elif sta == "RR34":
        if n_layer == 3:
            starting_model = np.zeros([n_layer+3, 4, iteration])
        
            starting_model[:,:,0]= np.array([[   10.,  1820.,  1750.,   340.],
                                         [  700.,  2550.,  5000.,  2700.],
                                         [ 1540.,  2850.,  6500.,  3700.],
                                         [ 4750.,  3050.,  7100.,  4050.],
                                         [10000.,  3100.,  7530.,  4190.],
                                         [10000.,  3100.,  7530.,  4190.]])
       
        
        elif n_layer == 6:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            
            starting_model[:,:,0] = np.array([[   10.,  1820.,  1750.,   340.],
                                              [  350.,  2550.,  5000.,  2700.],
                                              [  350.,  2550.,  5000.,  2700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.]])

        
        elif n_layer == 9:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[   10.,  1820.,  1750.,   340.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [10000.,  3230.,  7530.,  4370.],
                                              [10000.,  3230.,  7530.,  4370.]])
        
        else :
            n_layer = 12
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[   10.,  1820.,  1750.,   340.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [10000.,  3230.,  7530.,  4370.],
                                              [10000.,  3230.,  7530.,  4370.]])
            
    elif sta == "RR36":
        if n_layer == 3:
            starting_model = np.zeros([n_layer+3, 4, iteration])
        
            starting_model[:,:,0]= np.array([[  700.,  2550.,  5000.,  2700.],
                                             [ 1540.,  2850.,  6500.,  3700.],
                                             [ 4750.,  3050.,  7100.,  4050.],
                                             [10000.,  3230.,  7830.,  4370.],
                                             [10000.,  3230.,  7830.,  4370.],
                                             [10000.,  3230.,  7830.,  4370.]])
       
        
        elif n_layer == 6:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            
            starting_model[:,:,0] = np.array([[  350.,  2550.,  5000.,  2700.],
                                              [  350.,  2550.,  5000.,  2700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [10000.,  3230.,  7830.,  4370.],
                                              [10000.,  3230.,  7830.,  4370.],
                                              [10000.,  3230.,  7830.,  4370.]])

        
        elif n_layer == 9:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [10000.,  3230.,  7830.,  4370.],
                                              [10000.,  3230.,  7830.,  4370.],
                                              [10000.,  3230.,  7830.,  4370.]])
        
        else :
            n_layer = 12
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [10000.,  3230.,  7830.,  4370.],
                                              [10000.,  3230.,  7830.,  4370.],
                                              [10000.,  3230.,  7830.,  4370.]])
                        
    elif sta == "RR38":
        if n_layer == 3:
            starting_model = np.zeros([n_layer+3, 4, iteration])
        
            starting_model[:,:,0]= np.array([[  700.,  2550.,  5000.,  2700.],
                                             [ 1540.,  2850.,  6500.,  3700.],
                                             [ 4750.,  3050.,  7100.,  4050.],
                                             [10000.,  3100.,  7530.,  4190.],
                                             [10000.,  3100.,  7530.,  4190.],
                                             [10000.,  3230.,  7830.,  4370.]])
       
        
        elif n_layer == 6:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            
            starting_model[:,:,0] = np.array([[  350.,  2550.,  5000.,  2700.],
                                              [  350.,  2550.,  5000.,  2700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3230.,  7830.,  4370.]])

        
        elif n_layer == 9:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3230.,  7830.,  4370.]])
        
        else :
            n_layer = 12
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3230.,  7830.,  4370.]])

    elif sta == "RR40":
        if n_layer == 3:
            starting_model = np.zeros([n_layer+3, 4, iteration])
        
            starting_model[:,:,0]= np.array([[   10.,  1820.,  1750.,   340.],
                                             [  700.,  2550.,  5000.,  2700.],
                                             [ 1540.,  2850.,  6500.,  3700.],
                                             [ 4750.,  3050.,  7100.,  4050.],
                                             [10000.,  3100.,  7530.,  4190.],
                                             [10000.,  3100.,  7530.,  4190.]])
       
        
        elif n_layer == 6:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            
            starting_model[:,:,0] = np.array([[   10.,  1820.,  1750.,   340.],
                                              [  350.,  2550.,  5000.,  2700.],
                                              [  350.,  2550.,  5000.,  2700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.]])

        
        elif n_layer == 9:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[   10.,  1820.,  1750.,   340.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.]])
        
        else :
            n_layer = 12
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[   10.,  1820.,  1750.,   340.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.]])
            
    elif sta == "RR50":
        if n_layer == 3:
            starting_model = np.zeros([n_layer+3, 4, iteration])
        
            starting_model[:,:,0]= np.array([[  700.,  2550.,  5000.,  2700.],
                                         [ 1540.,  2850.,  6500.,  3700.],
                                         [ 4750.,  3050.,  7100.,  4050.],
                                         [10000.,  3100.,  7530.,  4190.],
                                         [10000.,  3100.,  7530.,  4190.],
                                         [10000.,  3230.,  7830.,  4370.]])
       
        
        elif n_layer == 6:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            
            starting_model[:,:,0] = np.array([[  350.,  2550.,  5000.,  2700.],
                                              [  350.,  2550.,  5000.,  2700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [  770.,  2850.,  6500.,  3700.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [ 2375.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3230.,  7830.,  4370.]])

        
        elif n_layer == 9:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  233.,  2550.,  5000.,  2700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [  513.,  2850.,  6500.,  3700.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [ 1583.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3230.,  7830.,  4370.]])
        
        else :
            n_layer = 12
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  175.,  2550.,  5000.,  2700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [  385.,  2850.,  6500.,  3700.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [ 1187.,  3050.,  7100.,  4050.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3100.,  7530.,  4190.],
                                              [10000.,  3230.,  7830.,  4370.]])
    elif sta == "RR52":
        sediment_thickness = 0 
        if n_layer == 3:
            starting_model = np.zeros([n_layer+3, 4, iteration])
        
            starting_model[:,:,0]= np.array([[  1170.,  2400.,  4700.,  2540.],
                                         [ 1170.,  2760.,  6300.,  3590.],
                                         [ 4560.,  2960.,  6900.,  3940.],
                                         [10000.,  3190.,  7730.,  4310.],
                                         [10000.,  3190.,  7730.,  4310.],
                                         [10000.,  3230.,  7830.,  4370.]])
       
        
        elif n_layer == 6:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            
            starting_model[:,:,0] = np.array([[  585.,  2400.,  4700.,  2540.],
                                              [  585.,  2400.,  4700.,  2540.],
                                              [  585.,  2760.,  6300.,  3590.],
                                              [  585.,  2760.,  6300.,  3590.],
                                              [ 2280.,  2960.,  6900.,  3940.],
                                              [ 2280.,  2960.,  6900.,  3940.],
                                              [10000.,  3190.,  7730.,  4310.],
                                              [10000.,  3190.,  7730.,  4310.],
                                              [10000.,  3230.,  7830.,  4370.]])

        
        elif n_layer == 9:
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[  390.,  2400.,  4700.,  2540.],
                                              [  390.,  2400.,  4700.,  2540.],
                                              [  390.,  2400.,  4700.,  2540.],
                                              [  390.,  2760.,  6300.,  3590.],
                                              [  390.,  2760.,  6300.,  3590.],
                                              [  390.,  2760.,  6300.,  3590.],
                                              [ 1520.,  2960.,  6900.,  3940.],
                                              [ 1520.,  2960.,  6900.,  3940.],
                                              [ 1520.,  2960.,  6900.,  3940.],
                                              [10000.,  3190.,  7730.,  4310.],
                                              [10000.,  3190.,  7730.,  4310.],
                                              [10000.,  3230.,  7830.,  4370.]])
        
        else :
            n_layer = 12
            starting_model = np.zeros([n_layer+3, 4, iteration])
            starting_model[:,:,0] = np.array([[  292.,  2400.,  4700.,  2540.],
                                              [  292.,  2400.,  4700.,  2540.],
                                              [  292.,  2400.,  4700.,  2540.],
                                              [  292.,  2400.,  4700.,  2540.],
                                              [  292.,  2760.,  6300.,  3590.],
                                              [  292.,  2760.,  6300.,  3590.],
                                              [  292.,  2760.,  6300.,  3590.],
                                              [  292.,  2760.,  6300.,  3590.],
                                              [ 1140.,  2960.,  6900.,  3940.],
                                              [ 1140.,  2960.,  6900.,  3940.],
                                              [ 1140.,  2960.,  6900.,  3940.],
                                              [ 1140.,  2960.,  6900.,  3940.],
                                              [10000.,  3190.,  7730.,  4310.],
                                              [10000.,  3190.,  7730.,  4310.],
                                              [10000.,  3230.,  7830.,  4370.]])
    # elif sta == "A422A":
    else:
        # sediment_thickness = 0 
        if n_layer == 3:
            starting_model = np.zeros([n_layer+4, 4, iteration])
        
        
            starting_model[:,:,0]= np.array([[ 2000.,  2400.,  4700.,  2000.],
                                             [ 2000.,  2760.,  6300.,  2500.],
                                             [ 3000.,  2760.,  6300.,  3200.],
                                             [ 3000.,  2960.,  6900.,  3800.],
                                             [10000.,  3190.,  7730.,  4310.],
                                             [10000.,  3190.,  7730.,  4310.],
                                             [10000.,  3230.,  7830.,  4370.]])
       
       
        
        elif n_layer == 6:
            starting_model = np.zeros([n_layer+5, 4, iteration])
            
            starting_model[:,:,0] = np.array([[  1000.,  2400.,  4700.,  2000.],
                                              [  1000.,  2400.,  4700.,  2000.],
                                              [  1000.,  2760.,  6300.,  2500.],
                                              [  1000.,  2760.,  6300.,  2500.],
                                              [ 1500.,  2960.,  6900.,  3200.],
                                              [ 1500.,  2960.,  6900.,  3200.],
                                              [ 1500.,  2960.,  6900.,  3800.],
                                              [ 1500.,  2960.,  6900.,  3800.],
                                              [10000.,  3190.,  7730.,  4310.],
                                              [10000.,  3190.,  7730.,  4310.],
                                              [10000.,  3230.,  7830.,  4370.]])

            
            
                    
    vs0 = np.zeros([1, int(np.sum(starting_model[:, 0, 0])), 1])
    for i in range(0, len(starting_model)):
        vs0[0, int(np.sum(starting_model[:, :, 0][:, 0][0:i])):int(np.sum(
            starting_model[:, :, 0][:, 0][0:i+1])), 0] = starting_model[:, 3][i][0]

    vsi = np.zeros([1, int(np.sum(starting_model[:, 0, 0])), 1])
    
    # print("Depth above half-space is "+ str(dep) + " m")
    # print("Depth above half-space is "+ str(dep - starting_model[len(starting_model)-1][0][0]) + " m")
    return(starting_model,vs0,vsi)

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
def plot_inversion_v2(starting_model,vs,mis_fit,ncompl,Data,likelihood_data,freq,sta,iteration,s,sigma_v,sigma_h,n_layer,alpha,burnin = 50000,mis_fit_trsh = 1):

    depth = np.arange(0, -int(np.sum(starting_model[0:-1, 0, 0])), -1)
    # var_vs = np.zeros([1,vs.shape[0]])

    # for i in range(0, vs.shape[0]):
    #     var_vs[0,i] = np.var(vs[i])
    layer_number  = -2
    Vs_final = np.zeros(vs[0].shape)
    N = 0
    vs_filtered = []
    
    for i in range(burnin,iteration):
        if mis_fit[0][i] < mis_fit_trsh:
                # if starting_model[layer_number][3][i] < 50000:
                Vs_final = Vs_final + vs[i]
                N = N + 1
                vs_filtered.append(vs[i])
    vs_filtered = np.array(vs_filtered)
    
    Vs_final = Vs_final/N
    print(N)

    plt.figure(dpi=300, figsize=(15, 25))
    plt.plot(Vs_final, depth, color='black', label='Final Result',linewidth=3)

    plt.xlabel('Shear Velocity [m/s]')
    plt.ylabel('Depth [m Below Seafloor]')
    plt.ylim([-12000, 0])
    plt.grid(True)
    plt.legend(loc='lower left')
    plt.tight_layout()


    dd = np.zeros([len(vs[0]),1000])

    for i in range(0,len(vs[0])):
            dd[i] = np.histogram(vs_filtered[:,i],bins=5000,range=([0,5000]))[0]
            dd[i] = dd[i] / np.max(dd[i])
    plt.rcParams.update({'font.size': 40})
    nn = int((iteration - burnin)/1000)+1 # I want to see just 1000 points of the data
    # nn = 1
    plt.figure(dpi=300, figsize=(35, 25))
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
    plt.ylim([10e-13,10e-11])

    plt.grid(True)
    plt.legend(loc='upper left',fontsize=30)
    # plt.ylim([1e-11,6e-11])
    # plt.yscale('log')
    
    plt.subplot(122)
    plt.title(sta)
    plt.imshow(dd,aspect="auto")
    plt.colorbar()
    # plt.plot(vs[0], depth, color='black', label='Start Model',linewidth= 5 ,linestyle='dashed')

    # Add labels, legend, and colorbar as in your code
    plt.xlabel('Shear Velocity [m/s]')
    plt.ylabel('Depth [m Below Seafloor]')
    plt.ylim([12000,0])
    plt.grid(True)
    plt.legend(loc='lower left')
    plt.tight_layout()





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
def plot_inversion(starting_model,vs,mis_fit,ncompl,Data,likelihood_data,freq,sta,iteration,s,sigma_v,sigma_h,n_layer,alpha,burnin = 50000,mis_fit_trsh = 1):
    linewidth_models = 0.05
    linewidth_compliance = 0.5 
    depth = np.arange(0, -int(np.sum(starting_model[0:-1, 0, 0])), -1)
    # var_vs = np.zeros([1,vs.shape[0]])

    # for i in range(0, vs.shape[0]):
    #     var_vs[0,i] = np.var(vs[i])
    layer_number  = -2
    Vs_final = np.zeros(vs[0].shape)
    N = 0
    for i in range(burnin,iteration):
        if mis_fit[0][i] < mis_fit_trsh:
                # if starting_model[layer_number][3][i] < 50000:
                Vs_final = Vs_final + vs[i]
                N = N + 1
    
    Vs_final = Vs_final/N
    print(N)
    
    plt.figure(dpi=300, figsize=(15, 25))
    plt.title(sta)

    plt.plot(Vs_final, depth, color='black', label='Final Result',linewidth=3)

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
            
    
    plt.rcParams.update({'font.size': 40})
    nn = int((iteration - burnin)/1000) # I want to see just 1000 points of the data
    # nn = 1
    plt.figure(dpi=300, figsize=(35, 25))
    plt.suptitle("Sigma V = " + str(sigma_v) + ", Sigma h = " + str(sigma_h) + ", Alpha = " + str(alpha) + ",N of Layer =" +str(n_layer))
    plt.subplot(121)

    for i in range(burnin, ncompl.shape[0], nn):
        if mis_fit[0][i] < mis_fit_trsh:
            plt.plot(freq, ncompl[i], color='green', linewidth=linewidth_compliance)
        
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Normalized Compliance')
        
        # plt.plot(freq, Data, color='black', label='Measured Compliance')
    plt.errorbar(freq, Data, yerr=s, ecolor=('black'), color='black',fmt='none',
                    linewidth= 5, label='Measured Compliance',capsize=15)
        
    # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue', 
    #             label='Median of Brun-in')
    plt.plot(freq, ncompl[1],color='black', label='Start Compliance',linewidth= 5,linestyle='dashed')
    # plt.ylim([10e-14,10e-10])
    # plt.ylim([10e-13,10e-11])

    plt.grid(True)
    plt.legend(loc='upper left',fontsize=30)
    # plt.ylim([1e-11,6e-11])
    # plt.yscale('log')
    
    plt.subplot(122)
    plt.title("YV."+sta)
    # plt.figure(dpi=300, figsize=(35, 25))
    for i in range(burnin, vs.shape[0]):
        if mis_fit[0][i] < mis_fit_trsh:
            plt.plot(vs[i], depth, color = 'green', linewidth = linewidth_models)
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
    plt.suptitle("YV."+sta)

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
def plot_inversion_beta(starting_model,vs,vs0,mis_fit,ncompl,Data,likelihood_data,freq,sta,iteration,s,sigma_v,sigma_h,n_layer,alpha,burnin = 50000,mis_fit_trsh = 1):
    linewidth_models = 0.05
    linewidth_compliance = 0.5 
    depth = np.arange(0, -int(np.sum(starting_model[0:-1, 0, 0])), -1)
    # var_vs = np.zeros([1,vs.shape[0]])

    # for i in range(0, vs.shape[0]):
    #     var_vs[0,i] = np.var(vs[i])
    layer_number  = -2
    Vs_final = np.zeros(vs[0].shape)
    N = 0
    for i in range(burnin,iteration):
        if mis_fit[0][i] < mis_fit_trsh:
                # if starting_model[layer_number][3][i] < 50000:
                Vs_final = Vs_final + vs[i]
                N = N + 1
    
    Vs_final = Vs_final/N
    print(N)
    
    plt.figure(dpi=300, figsize=(15, 25))
    plt.title(sta)

    plt.plot(Vs_final, depth, color='black', label='Final Result',linewidth=3)

    plt.xlabel('Shear Velocity [m/s]')
    plt.ylabel('Depth [m Below Seafloor]')
    plt.ylim([-12000, 0])
    plt.xlim([0, 5000])
    plt.grid(True)
    plt.legend(loc='lower left')
    plt.tight_layout()

    # plt.show()  # Don't forget to show the plot

    # plt.rcParams.update({'font.size': 30})
    # plt.figure(dpi=300, figsize=(10, 10))
    # for i in range(burnin, vs.shape[0], nn):
    #     if mis_fit[0][i] < mis_fit_trsh:
    #         plt.plot(vs[i], depth, color='grey', linewidth= 1)
            
    
    plt.rcParams.update({'font.size': 40})
    
    nn = int((iteration - burnin)/1000) # I want to see just 1000 points of the data
    # nn = 1
    plt.figure(dpi=300, figsize=(35, 25))
    plt.suptitle("Sigma V = " + str(sigma_v) + ", Sigma h = " + str(sigma_h) + ", Alpha = " + str(alpha) + ",N of Layer =" +str(n_layer))
    plt.subplot(121)

    for i in range(burnin, ncompl.shape[0], nn):
        if mis_fit[0][i] < mis_fit_trsh:
            if np.mean(vs[i]) < np.mean(vs0[:,0:len(vs[0])]):
                plt.plot(freq, ncompl[i], color='red', linewidth=linewidth_compliance)
            else:
                plt.plot(freq, ncompl[i], color='green', linewidth=linewidth_compliance)

        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Normalized Compliance')
        
        # plt.plot(freq, Data, color='black', label='Measured Compliance')
    plt.errorbar(freq, Data, yerr=s, ecolor=('black'), color='black',fmt='none',
                    linewidth= 5, label='Measured Compliance',capsize=15)
        
    # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue', 
    #             label='Median of Brun-in')
    plt.plot(freq, ncompl[1],color='black', label='Start Compliance',linewidth= 5,linestyle='dashed')
    # plt.ylim([10e-14,10e-10])
    # plt.ylim([10e-13,10e-11])

    plt.grid(True)
    plt.legend(loc='upper left',fontsize=30)
    # plt.ylim([1e-11,6e-11])
    # plt.yscale('log')
    
    plt.subplot(122)
    plt.title("YV."+sta)
    # plt.figure(dpi=300, figsize=(35, 25))
    for i in range(burnin, vs.shape[0]):
        if mis_fit[0][i] < mis_fit_trsh:
            if np.mean(vs[i]) < np.mean(vs0[:,0:len(vs[0])]):
                plt.plot(vs[i], depth, color = 'red', linewidth = linewidth_models)
            else:
                plt.plot(vs[i], depth, color = 'green', linewidth = linewidth_models)

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
    plt.suptitle("YV."+sta)

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

def plot_inversion_density(vs,vs0,mis_fit,Data,s,freq,sta,burnin,ncompl,iteration,mis_fit_trsh = 8):
    start_model = start_model_plot(sta)/50
    # start_model = start_model_plot_mean(sta)/50
    bins = 100
    jj = 0
    vs_good = []

    for ii in range(0,len(mis_fit[0])):
        if mis_fit[0][ii] < mis_fit_trsh:
            # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[ii])
                jj = jj+1
        
    vs_good = np.array(vs_good)
      
    linewidth_compliance = 0.5 
    nn = int((iteration - burnin)/1000) 
    
    # a = []
    # b = plt.hist(vs[:,100,0],bins=100,range=([0,5000]))[1]
    # vs_good = np.array(vs_good)
    
    # for i in range(0,len(vs[0])):
    #     a.append(plt.hist(vs_good[:,i,0],bins=100,range=([0,5000]))[0])
    #     print(len(vs[0]) - i)
        
        
        
    # plt.figure(dpi=300,figsize=(20,20))
    # selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
    # selected_labels = b[selected_indices]
    # selected_labels = selected_labels.astype(int)
    
    # plt.imshow(a / 1000, aspect='auto', cmap='jet', norm=plt.Normalize(vmin=0, vmax=1.5))
    
    # # Setting custom y-axis ticks and labels
    # plt.xticks(ticks=selected_indices, labels=selected_labels)
    
    # plt.colorbar()
    # plt.show()



    vs_good = vs_good[:,0:10000,0]
    
    downsample_factor = vs_good.shape[1] // 100
    
    # Reshape the array to prepare for averaging
    # New shape will be (20000, 100, downsample_factor, 1)
    vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)
    
    # Take the mean along the downsample_factor dimension
    vs_downsampled = vs_reshaped.mean(axis=2)
    
    print(vs_downsampled.shape)  # This should print (20000, 100, 1)
    
    # a = np.zeros([100,100])
    c = np.zeros([bins,bins])
    b = plt.hist(vs_reshaped[:,10,0],bins=bins,range=([0,5000]))[1]
    
    for i in range(0,len(vs_downsampled[0])):
        
        c[i] = plt.hist(vs_downsampled[:,i],bins=bins,range=([0,5000]),density=True,log=False)[0]
        c[i] = c[i]/np.max(c[i])
        print(i)
        
    
    from matplotlib.colors import LinearSegmentedColormap
    
    # Define the colors for the colormap (from white to red to black)
    colors = ["white", "red", "black"]
    
    # Define the transition points for the colors
    n_bins = [0, 0.4, 1]  # You can adjust these thresholds based on your data
    
    # Create the custom colormap
    custom_colormap = LinearSegmentedColormap.from_list("custom_cmap", list(zip(n_bins, colors)))
    
        
    plt.figure(dpi=300,figsize=(30,20))
    plt.suptitle("YV."+str(sta))
    plt.subplot(121)
    for i in range(burnin, ncompl.shape[0], nn):
        if mis_fit[0][i] < mis_fit_trsh:
            plt.plot(freq, ncompl[i], color='green', linewidth=linewidth_compliance)

    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Normalized Compliance')
        
        # plt.plot(freq, Data, color='black', label='Measured Compliance')
    plt.errorbar(freq, Data, yerr=s, ecolor=('black'), color='black',fmt='none',
                    linewidth= 5, label='Measured Compliance',capsize=15)
    # plt.plot(freq, ncompl[1],color='blue', label='Starting Compliance',linewidth= 5,linestyle='dashed')
    plt.legend(loc='upper left',fontsize=30)
    
    plt.subplot(122)
    selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
    selected_labels = b[selected_indices]
    selected_labels = selected_labels.astype(int)
    
    depth = np.arange(0,100,1 )
    
    ff = plt.imshow(c , aspect='auto', cmap=custom_colormap, norm=plt.Normalize(vmin=0, vmax=1))
    # plt.plot(start_model,depth,color='blue',linewidth=8,linestyle='dashed',label="Starting Model ")
    
    
    # plt.imshow(a , aspect='auto', cmap='jet')

    # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue', 
    #             label='Median of Brun-in')
    # plt.ylim([10e-14,10e-10])
    # plt.ylim([10e-13,10e-11])

    plt.grid(True)
    # plt.ylim([1e-11,6e-11])
    # plt.yscale('log')
    
    # Setting custom y-axis ticks and labels
    plt.xticks(ticks=selected_indices, labels=selected_labels)
    cbar = plt.colorbar(ff)
    # depth = np.arange(0, -10000, -1)

    # plt.plot(vs[0][0:10000], depth/100, color='black', label='Start Model',linewidth= 5 ,linestyle='dashed')

    selected_indices_y = [0,12,24,36,48,60,72,84]
    selected_labels_y = [0,-2000,-4000,-6000,-8000,-10000,-12000,-14000]
    plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m]")
    cbar.set_label('Probability')
    plt.legend(loc='lower left',fontsize=35)
    plt.ylim(72,0)
    # plt.colorbar()
    plt.tight_layout()
    plt.show()
    
#%%

def plot_inversion_density_all(Inversion_container):
    
    plt.figure(dpi=300,figsize=(40,25))
    
    for ii in range(0,len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        # vs0 = Inversion_container[ii]["Shear Velocity Starting"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        # Data = Inversion_container[ii]["compliance Measured"]
        # s = Inversion_container[ii]["uncertainty"]
        # freq = Inversion_container[ii]["compliance Frequency"]
        sta = Inversion_container[ii]["Station"]
        # burnin = Inversion_container[ii]["burnin"]
        # ncompl = Inversion_container[ii]["compliance Forward"]
        # iteration = Inversion_container[ii]["iteration"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]
    
    
        start_model = start_model_plot(sta)/50
        
        start_model = start_model_plot_mean(sta)/50 #mean of all models,sedimental models differ from rocky models
        
        # start_model = start_model_plot_mean(sta)/50
        bins = 100
        jj = 0
        vs_good = []
        
        for i in range(0,len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:
    
                    vs_good.append(vs[i])
                    jj = jj+1
            
        vs_good = np.array(vs_good)
          
        vs_good = vs_good[:,0:10000,0]
        
        downsample_factor = vs_good.shape[1] // 100
        
        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)
        
        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)
        
        print(vs_downsampled.shape)  # This should print (20000, 100, 1)
        
        # a = np.zeros([100,100])
        c = np.zeros([bins,bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]
        
        for i in range(0,len(vs_downsampled[0])):
            
            # c[i] = plt.hist(vs_downsampled[:,i],bins=bins,range=([0,5000]),density=True,log=False)[0]
            c[i] = np.histogram(vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True)[0]
            c[i] = c[i]/np.max(c[i])
            print(i)
            
        
        from matplotlib.colors import LinearSegmentedColormap
        
        # Define the colors for the colormap (from white to red to black)
        colors = ["white", "red", "black"]
        
        # Define the transition points for the colors
        n_bins = [0, 0.4, 1]  # You can adjust these thresholds based on your data
        
        # Create the custom colormap
        custom_colormap = LinearSegmentedColormap.from_list("custom_cmap", list(zip(n_bins, colors)))
   
        plt.subplot(2,4,int(ii+1))
        plt.title(str("YV.")+sta)
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)
        
        depth = np.arange(0,100,1 )
        
        ff = plt.imshow(c , aspect='auto', cmap=custom_colormap, norm=plt.Normalize(vmin=0, vmax=1))
        plt.plot(start_model,depth,color='blue',linewidth=8,linestyle='dashed',label="Mean Starting Model ")
        # plt.imshow(a , aspect='auto', cmap='jet')
    
        # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue', 
        #             label='Median of Brun-in')
        # plt.ylim([10e-14,10e-10])
        # plt.ylim([10e-13,10e-11])
    
        plt.grid(True)
        # plt.ylim([1e-11,6e-11])
        # plt.yscale('log')
        
        # Setting custom y-axis ticks and labels
        plt.xticks(ticks=selected_indices, labels=selected_labels)
        cbar = plt.colorbar(ff)
        # depth = np.arange(0, -10000, -1)
    
        # plt.plot(vs[0][0:10000], depth/100, color='black', label='Start Model',linewidth= 5 ,linestyle='dashed')
    
        selected_indices_y = [0,12,24,36,48,60,72,84]
        selected_labels_y = [0,-2000,-4000,-6000,-8000,-10000,-12000,-14000]
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        if ii == 0 or ii == 4:
            plt.ylabel("Depth [m]")
        if ii == 4 or ii == 5 or ii == 6 or ii == 7:
            plt.xlabel("Shear Velocity [m/s]")
        if ii == 3 or ii == 7:
            cbar.set_label('Probability')
        if ii == 7:
            plt.legend(loc='lower left',fontsize=30)
        plt.ylim(48,0)
        # plt.ylim(72,0)
        # plt.colorbar()
    plt.tight_layout()
    # file_path_save = "/Users/mohammadamin/Desktop"
    # plt.savefig(file_path_save + "Inversion_Table.pdf")

#%%
def start_model_plot_mean(sta):
    start_model = np.zeros(100)
    if sta == "RR28" or sta == "RR29" or sta == "RR34":
        start_model[0:2] = 340
        start_model[2:7] = 2540
        start_model[7:16] = 3590
        start_model[16:44] = 3940
        start_model[44:100] = 4310
        
    # elif  sta == "RR29":
    #     start_model[0:1] = 340
    #     start_model[1:6] = 2540
    #     start_model[6:15] = 3590
    #     start_model[15:43] = 3940
    #     start_model[43:100] = 4310

    elif  sta == "RR52" or sta == "RR50" or sta == "RR40" or sta == "RR38" or sta == "RR36"  :
        start_model[0:5] = 2540
        start_model[5:14] = 3590
        start_model[14:42] = 3940
        start_model[42:100] = 4310
        
    return(start_model)
#%%
def start_model_plot(sta):
    start_model = np.zeros(100)
    if sta == "RR28":
        start_model[0:1] = 340
        start_model[1:6] = 2700
        start_model[6:16] = 3700
        start_model[16:45] = 4050
        start_model[45:99] = 4510
        
    elif  sta == "RR29":
        start_model[0:1] = 340
        start_model[1:5] = 2700
        start_model[5:14] = 3700
        start_model[14:44] = 4050
        start_model[44:99] = 4510
        
    elif  sta == "RR34":
        start_model[0:2] = 340
        start_model[2:7] = 2700
        start_model[7:15] = 3700
        start_model[15:42] = 4050
        start_model[42:99] = 4500  

    elif  sta == "RR36":
        start_model[0:4] = 2700
        start_model[4:15] = 3700
        start_model[15:42] = 4050
        start_model[42:99] = 4370  
        
    elif  sta == "RR38":
        start_model[0:4] = 2700
        start_model[4:15] = 3700
        start_model[15:42] = 4050
        start_model[42:99] = 4190  
                
    elif  sta == "RR40":
        start_model[0:4] = 2700
        start_model[4:13] = 3700
        start_model[13:42] = 4050
        start_model[42:100] = 4360      

    elif  sta == "RR50":
        start_model[0:7] = 2540
        start_model[7:14] = 3590
        start_model[14:42] = 3940
        start_model[42:100] = 4190 

    elif  sta == "RR52":
        start_model[0:7] = 2660
        start_model[7:14] = 3670
        start_model[14:42] = 4022
        start_model[42:100] = 4330
        
    else:
    # elif  sta == "A422A":
        start_model[0:13] = 2200
        start_model[13:25] = 2700
        start_model[25:42] = 3500
        start_model[42:60] = 3800
        start_model[60:100] = 4310
    return(start_model)

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
