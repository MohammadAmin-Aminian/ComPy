#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 14:04:23 2023

@author: mohammadamin
"""

import numpy as np
import matplotlib.pyplot as plt
def invert_compliace(Data,f,depth_s,uncertainty,starting_model = None,iteration = 20000,perturbation = 1e-2):
# With Liklihood thickness changing + Constrain for thickness
# Accept Probablity = delta s / s ^ 2 (L[i]/L[i+1])
    # Data = scipy.signal.savgol_filter(Data,12,2)
    
    #Data uncertainty 
    s = np.sqrt(np.var(Data))
    # s = np.std(Data)
    # s = np.mean(uncertainty) 
    starting_model,vs0,vsi = model_exp(iteration,first_layer=100,n_layer=7,power_factor=2)
    #Constrains
    const_vs_lower = 10 # +-20% of the vs for model constrains
    const_vs_higher = 10 # +-20% of the vs for model constrains

    const_th_lower = 1 # +-20% of the vs for model constrains
    const_th_higher = 10 # +-20% of the vs for model constrains

    likelihood_Model = np.zeros([1, iteration])

    likeli_hood = np.zeros([1, iteration])
    likelihood_data = np.zeros([1, iteration])
    likeli_hood[0] = 0
    mis_fit_Model = np.zeros([1, iteration])
    mis_fit = np.zeros([1, iteration])
    accept = 0


# s = np.sqrt(np.sum((Data - calc_norm_compliance(
    # depth, f, starting_model[:, :, 0]))**2) / len(Data))

    ncompl = np.zeros([iteration, np.shape(f)[0]])

    for i in range(1, iteration):
        print(iteration - i)
        r = np.random.uniform(-1, 1)*perturbation  # 1 % perturbation
        starting_model[:, :, i] = starting_model[:, :, i-1]
        layer = np.random.randint(0, starting_model.shape[0] -1)

        if np.random.randint(0,2) == 0: # Choosing between layer thickness or layer VS
    
            # Vs Changes based on perturbation
            starting_model[layer, 3, i] = starting_model[layer, 3, i-1]*(1+r)
    
            #Constrains VS
            if starting_model[layer, 3, i-1]*(1+r) < starting_model[layer, 3, 0]* (1-const_vs_lower):
                print("Lower limit of Velocity")
                r = -r/10
                starting_model[layer, 3, i] = starting_model[layer, 3, i-1]*(1+r)
     
            if starting_model[layer, 3, i-1]*(1+r) > starting_model[layer, 3, 0]* (1+const_vs_higher):
                print("Upper limit of Velocity")
                r = -r/10
                starting_model[layer, 3, i] = starting_model[layer, 3, i-1]*(1+r)        

            # Vp changes based on Vs perturbation (poisson ratio)
            starting_model[layer, 2, i] = velp(starting_model[layer, 3, i])
            
            # Density Changes based on Vs perturbation(Gardner's Eq)
            starting_model[layer, 1, i] = density(starting_model[layer, 3, i])
            
            #updating the S with every iteration
        
        else:
            # r = np.random.randint(0, 10)/100
            
            layer = np.random.randint(0, starting_model.shape[0])
            
            layer = 0 # Just sediment layer can change in thickness
            starting_model[layer, 0, i] = starting_model[layer, 0, i-1]*(1+(1*r)) # added to rectified the bug cause by model issue

            # layer1 = starting_model.shape[0]-1

            starting_model[starting_model.shape[0] -2, 0, i] = starting_model[starting_model.shape[0] -2, 0, i] + starting_model[layer, 0, i-1]*((1*-r))
            


             # starting_model[layer1, 0, i] = starting_model[layer1, 0, i-1]*(1+(1*r))
         
            if starting_model[layer, 0, i-1]*(1+1*r) < starting_model[layer, 0, 0]* (1-const_th_lower):
             
                    r = -r/10  # 1 % perturbation
             
                    starting_model[layer, 0, i] = starting_model[layer, 0, i-1]*(1+1*r)
             

                    # starting_model[layer1, 0, i] = starting_model[layer1, 0, i-1]*(1+(1*r))
                    print("Lower limit of thickness")

            if starting_model[layer, 0, i-1]*(1+1*r) > starting_model[layer, 0, 0]* (1+const_th_higher):
                 
                 r = -r/10  # 1 % perturbation
                 
                 starting_model[layer, 0, i] = starting_model[layer, 0, i-1]*(1+1*r)
                    
                    # starting_model[layer1, 0, i] = starting_model[layer1, 0, i-1]*(1+(1*r))
                 print("Upper limit of thickness")
        ncompl[i, :] = calc_norm_compliance(depth_s, f, starting_model[:, :, i])

        for ii in range(0, len(starting_model)):
            vsi[0, int(np.sum(starting_model[:, :, i][:, 0][0:ii])):int(np.sum(
                starting_model[:, :, i][:, 0][0:ii+1])), 0] = starting_model[:, 3][ii][i]

        likelihood_Model[0, i] = liklihood(vs0, vsi,s = np.std(vs0[0]))  # Model Likelihood 
    
        likelihood_data[0, i] = liklihood(Data, ncompl[i, :], s=s)
    
        # likeli_hood[0, i] = likelihood_data[0, i]/2 + (likelihood_Model[0, i])  # Likihood Data + Model
        # Likelihood Model+Data Based on Adreas Fischner

        likeli_hood[0, i] = liklihood(Data, ncompl[i, :], s=s)  # Liklihood Data
        
        likeli_hood[0, i] = liklihood_roughness(Data, ncompl[i, :],vs0, s=s,alpha=2)  # Liklihood Data + Roughness

        # Likelihood Based on Mariano 


        mis_fit[0, i] = misfit(Data, ncompl[i, :],s=s)  # Misfit Data
        # mis_fit_l1[0, i] = misfit(Data, ncompl[i, :],s=s, l=1)  # Misfit Data L^1
    
        # s = np.sqrt(np.sum((Data - ncompl[i, :])**2)/ len(Data))
        # Updating S based on Thomas Bodin 
        # print(s)
        print(mis_fit[0, i])
    
        if likeli_hood[0, i] >= likeli_hood[0, i-1]:
            # starting_model[:, :, i] = starting_model[:, :, i]
            accept+=1

        else:
            # Probablity of Acceptance
            p_candidate = np.random.rand(1)[0]
            # print((likeli_hood[0,i]/likeli_hood[0,i-1]))

            if p_candidate < (likeli_hood[0,i]/likeli_hood[0,i-1])  :
                accept+=1
                # starting_model[:, :, i] = starting_model[:, :, i]
            else:                                   #New Line
                  starting_model[:, :, i] = starting_model[:, :, i-1]
                  likeli_hood[0, i] = likeli_hood[0, i-1]
                  mis_fit[0, i] = mis_fit[0, i-1]
                  likelihood_Model[0, i] = likelihood_Model[0, i-1]
                  likelihood_data[0, i] = likelihood_data[0, i-1]

    vs = np.zeros([iteration, int(np.sum(starting_model[0:-1, 0, 0])), 1])

    for j in range(0, iteration):
        print(iteration - j)
        for i in range(0, len(starting_model)):
            vs[j, int(np.sum(starting_model[:, :, j][:, 0][0:i])):int(np.sum(
                starting_model[:, :, j][:, 0][0:i+1])), 0] = starting_model[:, 3,j][i]
            
            accept/iteration
    plot_inversion(starting_model,vs,mis_fit,ncompl,Data,likelihood_data=likeli_hood,freq=f,
                   iteration=iteration,s=s,burnin = 50000)

#%%

def model_exp(iteration, first_layer = 200, n_layer = 15,power_factor = 1.15):
    starting_model = np.zeros([n_layer, 4, iteration])
    dep = 0
    for i in range(0,len(starting_model)):
        starting_model[i][0][0] = np.float16(first_layer*(power_factor)**i) # Thickness
    
        starting_model[0][3][0] = 350 # VS
        # starting_model[1][3][0] = 350 # VS
    
        print(starting_model[i][0][0])
    
        dep = starting_model[i][0][0] + dep
    
        starting_model[i][3][0] = 1200 + ((i/len(starting_model))*3000) # VS
    
        # starting_model[i][3][0] = 2000

        starting_model[i][2][0] = velp(starting_model[i][3][0]) # Vp
    
        starting_model[i][1][0] = density(starting_model[i][3][0]) # Density
        #RR38
    
    #Half space
    starting_model[len(starting_model)-1][0][0] = 6000000 # Thickness



    
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
    print(dep)
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
     # misfit = np.sum(((d-m)/s)**l)
     misfit = np.sum(((d-m)**l)/(s**2))
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
     L =   k * np.exp(-0.5*np.linalg.norm(d-m)**2/(s**2))

     return(L)
 #%%
def liklihood_roughness(d,m,vs,k=1,s=1,alpha=1):
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
      R_m = alpha*Roughness(vs.flatten())
      L =   k * np.exp(-0.5*(np.linalg.norm(d-m)**2/(s**2)+R_m))

      return(L)
  
#%%
def Roughness(vs):
    R_m = np.sqrt(np.sum(np.gradient(vs,2)))
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
def raydep(P,om,d,ro,vp2,vs2):
    
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
def plot_inversion(starting_model,vs,mis_fit,ncompl,Data,likelihood_data,freq,iteration,s,burnin = 50000):

    depth = np.arange(0, -int(np.sum(starting_model[0:-1, 0, 0])), -1)

    # var_vs = np.zeros([1,vs.shape[0]])

    # for i in range(0, vs.shape[0]):
    #     var_vs[0,i] = np.var(vs[i])
    
    Vs_final = np.zeros(vs[0].shape)
    N = 0
    for i in range(burnin,iteration):
        Vs_final = Vs_final + vs[i]
        N = N + 1
    
    Vs_final = Vs_final/N
    print(N)

    plt.rcParams.update({'font.size': 30})

    nn = int((iteration - burnin)/500) # I want to see just 1000 points of the data

    plt.figure(dpi=300, figsize=(25, 20))
    # plt.subplot(221)

    # # plt.plot(likeli_hood[0, 1:iteration-1],color='Blue',label="Total Likelihood")
    # # plt.plot(likelihood_Model[0, 1:iteration-1],color='red',label="Model Likelihood")
    # plt.plot(likelihood_data[0, 1:iteration-1],color='green',label="Data Likelihood")
    # plt.xscale('log')
    # plt.vlines(x=burnin, ymin=0, ymax=1, color='r',
    #            label='Burn-in Region', linestyles='dashed')
    # plt.legend(loc='upper left')
    
    # plt.title('Likelihood')
    # plt.xlabel('Iteration')
    # plt.ylabel('Liklihood')
    # plt.grid(True)
    
    plt.subplot(221)
    for i in range(burnin, ncompl.shape[0], nn):
        plt.plot(freq, ncompl[i], color='red', linewidth=0.25)
        
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Normalized Compliance')
        
        # plt.plot(freq, Data, color='black', label='Measured Compliance')
    plt.errorbar(freq, Data, yerr=s, ecolor=('r'), color='black',
                    linewidth=3, label='Measured Compliance')
        
    plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue', 
                label='Median of Brun-in')
    plt.plot(freq, ncompl[1],color='green', label='Start Compliance',linewidth= 3)
        
    plt.grid(True)
    plt.legend(loc='upper left')
    # plt.ylim([1e-11,6e-11])
        
    plt.subplot(223)
    plt.plot(mis_fit[0,1:iteration-1])
    plt.xscale('log')
    # plt.yscale('log')
    
    plt.vlines(x=burnin, ymin=0, ymax=np.max(mis_fit[0,1:iteration-1]), color='r',
               label='Burn-in Region', linestyles='dashed')
    plt.legend(loc='upper left')
    
    plt.title('Data Misfit (L^2)')
    plt.xlabel('Iteration')
    plt.ylabel('Error')
    plt.grid(True)
    

        
    plt.subplot(122)     
    for i in range(burnin, vs.shape[0], nn):
        plt.plot(vs[i], depth, color='grey', linewidth=0.5)
            
        # plt.plot(np.mean(vs[burnin:iteration],axis=0), depth, color='blue', 
        #          label='Median of Burn-in ')
    plt.plot(Vs_final, depth, color='black', label='Final Result',linewidth=3)
    plt.plot(vs[0], depth, color='green', label='Start Model',linewidth=3,linestyle='dashed')
    
    plt.grid(True)
    plt.xlabel('Shear Velocity [m/s]')
    plt.ylabel('Depth [m]')
    # plt.ylim([-int(np.sum(starting_model[:, 0, 0]))+2000, 0])
    
    plt.ylim([-7000, 0])
    # plt.xlim([0,4500])
    
    vs_burnin = np.zeros([iteration, int(np.sum(starting_model[:, 0, 0])), 1])
    plt.legend(loc='lower left')
    plt.tight_layout()
    
    plt.rcParams.update({'font.size': 30})
    plt.figure(dpi=300, figsize=(10, 10))
    plt.plot(Vs_final, depth, color='black', label='Final Result',linewidth=3)
    plt.plot(vs[0], depth, color='green', label='Start Model',linewidth=3,linestyle='dashed')
    
    plt.grid(True)
    plt.xlabel('Shear Velocity [m/s]')
    plt.ylabel('Depth [m]')
    # plt.ylim([-int(np.sum(starting_model[:, 0, 0]))+2000, 0])
    
    plt.ylim([-10000, 0])
    # plt.xlim([0,4500])
    
    vs_burnin = np.zeros([iteration, int(np.sum(starting_model[:, 0, 0])), 1])
    plt.legend(loc='lower left')
    plt.tight_layout()