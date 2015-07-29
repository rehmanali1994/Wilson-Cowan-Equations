# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 17:25:06 2015

@author: rehmanali
"""

import time
from WilsonCowanSystemsFFT import *

def WCSimulationDiagram(TE):
    " CONSTANTS: DO NOT CHANGE! "
    BETA = 50; TI = 0.4; SE = 15; #TE=0.125
    AEE = 1;  AIE = 1; AEI = 1.5; AII = 0.25;
    
    " Defining ODE Parameters "
    L = 400; # Length of mesh
    span = 180; # Point to Left and Right of Central Maximum in Kernel
    nx = 2*span + 1; # Number of finite elements in mesh
    dx = L/(nx-1); # Spacing of mesh
    dx_kern = dx; # Spacing of Points in the Kernel
    dt = 0.1; tmore = 200; # Time Step and Initlal Integration Time
    tshow = 10; # Amount of Time Data to Display
    kernType='Exponential' # 'Gaussian' or 'Exponential'
    mode = 'wrap'; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
    
    " Initial Conditions "
    t0 = 0; u0 = 0.03*np.random.randn((nx))+0.41; 
    v0 = 0.03*np.random.randn((nx))+0.21;
    
    def WC1DFFTSim(TAU,SI_SE_ratio):
        " Creating the WilsonCowan1D object "
        TAU = TAU; # Time Constant: Remaining ODE Parameters SI=8; TAU=0.6
        SI = SI_SE_ratio*SE; # Ensuring Spatial Coupling Proportionality
        pardict = {'BETA':BETA, 'TE':TE, 'TI':TI, 'SE':SE, 'AEE': AEE, 'AIE':AIE,
                   'AEI':AEI, 'AII':AII, 'L':L, 'span':span, 'dx_kern':dx_kern,
                   'SI':SI, 'TAU': TAU, 'dt':dt, 'kernType':kernType, 'mode':mode};
        WC1D = WilsonCowan1D_FFT(pardict=pardict);
                  
        " Initial Conditions for ODE "
        WC1D.setInitConds(t0,u0,v0,tmore=tmore,tshow=tshow);
        
        #WC1D.interactiveIO_img();
        tstart = time.time();
        WC1D.WilsonCowanIntegrator();
        locs = np.arange(WC1D.tvals.size); 
        locs = locs[WC1D.tvals>(WC1D.tvals[-1]-WC1D.tshow)];
        ydisp = WC1D.yvals[:,locs];
        udisp = ydisp[:WC1D.nx,:]; vdisp = ydisp[WC1D.nx:(2*WC1D.nx),:]; 
        patVar = np.mean(np.var(udisp,axis=0)+np.var(vdisp,axis=0));
        telapsed = time.time()-tstart; 
        print('Elapsed time is: '+str(telapsed));
        print('Pattern Variance is: '+str(patVar));
        #WC1D.imgWilsonCowan1D(tdisp, WC1D.xmesh, ydisp, 0); 
        return patVar;
    
    def patFormArr(SI_range,tau_range):
        patVarArr = np.zeros((tau_range.size,SI_range.size));
        for SI_index in np.arange(SI_range.size):
            for tau_index in np.arange(tau_range.size):
                SI = SI_range[SI_index]; tau = tau_range[tau_index];
                print('SI = '+str(SI)+', TAU = '+str(tau));
                patVarArr[tau_index,SI_index] = WC1DFFTSim(tau,SI)
        return patVarArr
                   
    tau_range = np.linspace(0.25,1.05,81); SI_range = np.linspace(0.02,1.6,80);  
    #tau_range = np.linspace(0.2,0.8,81); # ONLY FOR TE = 0.125
    tstart = time.time(); patVarArr = patFormArr(SI_range,tau_range); 
    telapsed = time.time()-tstart; np.savetxt('TE_'+str(TE)+'.txt',patVarArr); 
    print('TOTAL Elapsed time is: '+str(telapsed));

#tstart = time.time();    
#WC1D = WC1DFFTSim(0.6,2/3);
#telapsed = time.time()-tstart;
#disp(telapsed)