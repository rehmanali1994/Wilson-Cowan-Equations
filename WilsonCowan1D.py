# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 22:14:46 2015

@author: rehmanali
"""

from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage import convolve as conv
from scipy import signal
from scipy.integrate import ode
import matplotlib.pyplot as plt
from matplotlib import cm
import visvis as vv
import numpy as np
import pandas as pd

class WilsonCowan1D:
    'Organized Structure for 1D Wilson Cowan Equations Parameters and Results'
    tshow_default = 10; tmore_default = 10;    
    
    " Object Constructor "
    def __init__(self,pardict=None,filename=None):
        if (filename is None) and (pardict is None):
            raise Exception('Must specify filename or pardict but not both')
        if filename is None:
            " Constant Params for ODE "
            self.BETA = pardict['BETA']; self.TE = pardict['TE']; 
            self.TI = pardict['TI']; self.SE = pardict['SE'];
            self.AEE = pardict['AEE']; self.AIE = pardict['AIE']; 
            self.AEI = pardict['AEI']; self.AII = pardict['AII'];
    
            " Defining Other ODE Parameters "
            self.L = pardict['L']; # Length of mesh
            self.nx = pardict['nx']; # Spaces in mesh
            self.dx = self.L/(self.nx-1); # Spacing of mesh
            self.xmesh = np.linspace(0,self.L,self.nx); # the mesh itself
            self.span = pardict['span']; # Point to Left and Right of Central Maximum in Kernel
            self.dx_kern = pardict['dx_kern']; # Spacing of Points in the Kernel
            self.SI = pardict['SI']; self.TAU = pardict['TAU']; # Remaining ODE Parameters
            self.dt = pardict['dt']; self.WCSystem = None;
            self.kernType= pardict['kernType']; # 'Gaussian' or 'Exponential'
            self.mode = pardict['mode']; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
        elif pardict is None:
            " Retrieving saved state from file "
            xl = pd.ExcelFile(filename);
            dfs = {sheet: xl.parse(sheet) for sheet in xl.sheet_names}
            paramFields = {'BETA','TE','TI','SE','SI','TAU','L','nx','dt',
                'AEE','AIE','AEI','AII','span','dx_kern','mode','kernType'}
            pardict = {'BETA':None,'TE':None,'TI':None,'SE':None,'AEE': None, 
                       'AIE':None,'AEI':None,'AII':None, 'L':None, 'nx':None, 
                       'span':None,'dx_kern':None,'SI':None,'TAU':None,
                       'dt':None, 'kernType':None, 'mode':None};
            for names in paramFields:
                pardict[names] = dfs['params'][names]['Parameters'];
            self.BETA = pardict['BETA']; self.TE = pardict['TE']; 
            self.TI = pardict['TI']; self.SE = pardict['SE'];
            self.AEE = pardict['AEE']; self.AIE = pardict['AIE']; 
            self.AEI = pardict['AEI']; self.AII = pardict['AII'];
            self.L = pardict['L']; # Length of mesh
            self.nx = pardict['nx']; # Spaces in mesh
            self.dx = self.L/(self.nx-1); # Spacing of mesh
            self.xmesh = np.linspace(0,self.L,self.nx); # the mesh itself
            self.span = pardict['span']; # Point to Left and Right of Central Maximum in Kernel
            self.dx_kern = pardict['dx_kern']; # Spacing of Points in the Kernel
            self.SI = pardict['SI']; self.TAU = pardict['TAU']; # Remaining ODE Parameters
            self.dt = pardict['dt']; self.WCSystem = None;
            self.kernType= pardict['kernType']; # 'Gaussian' or 'Exponential'
            self.mode = pardict['mode']; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
            self.xmesh = dfs['xmesh'][0].values;
            self.tvals = dfs['tvals'][0].values;
            uvals = np.transpose(dfs['uvals'].values); 
            vvals = np.transpose(dfs['vvals'].values); 
            self.yvals = np.concatenate((uvals,vvals));
            self.WCSystem = ode(self.WilsonCowanODE).set_integrator('dopri');
            self.WCSystem.set_initial_value(self.yvals[:,-1],self.tvals[-1]);
        else: raise Exception('Must specify filename or pardict but not both');
            
    " Must Set Initial Conditions before Integration (WilsonCowanIntegrator) "    
    def setInitConds(self,t0,u0,v0,tmore=tmore_default,tshow=tshow_default):
        self.y0 = np.concatenate((u0,v0)); self.t0 = t0; 
        self.tmore = tmore; self.tshow = tshow;
        
    " Defining the Wilson Cowan Equations "
    def WilsonCowanODE(self,t,y):
        u = y[0:self.nx]; v = y[self.nx:(2*self.nx)];
        kerE = self.kern(self.SE,self.span,self.dx_kern); 
        kerI = self.kern(self.SI,self.span,self.dx_kern);
        Lu = self.AEE*conv(u,kerE,mode=self.mode)-  \
            self.AEI*conv(v,kerI,mode=self.mode)-self.TE;
        Lv = self.AIE*conv(u,kerE,mode=self.mode)-  \
            self.AII*conv(v,kerI,mode=self.mode)-self.TI;
        du = -u + self.f(Lu); dv = (-v + self.f(Lv))/(self.TAU);
        return np.concatenate((du,dv));    
    def WilsonCowanIntegrator(self):
        self.WCSystem = ode(self.WilsonCowanODE).set_integrator('dopri');
        self.WCSystem.set_initial_value(self.y0,self.t0);
        self.tvals, self.yvals = [self.t0], self.y0;
        while self.WCSystem.successful() and self.WCSystem.t < self.tmore:
            self.WCSystem.integrate(self.WCSystem.t+self.dt);
            self.tvals = np.concatenate((self.tvals,[self.WCSystem.t]));
            self.yvals = np.column_stack((self.yvals,self.WCSystem.y));
        return self.tvals, self.yvals, self.WCSystem    
    def WCIntegrateLast(self):
        tlast = self.tvals[-1]; tend = tlast+self.tmore;
        while self.WCSystem.successful() and self.WCSystem.t < tend:
            self.WCSystem.integrate(self.WCSystem.t+self.dt);
            self.tvals = np.concatenate((self.tvals,[self.WCSystem.t]));
            self.yvals = np.column_stack((self.yvals,self.WCSystem.y));
        return self.tvals, self.yvals, self.WCSystem

    " Plotting the Results of the Wilson Cowan Equations "
    def imgWilsonCowan1D(self,tvals, xmesh, yvals, pltType):
        dt = np.mean(np.diff(tvals)); dx = np.mean(np.diff(xmesh));
        uvals = yvals[0:self.nx,:]; vvals = yvals[self.nx:(2*self.nx),:];
        if pltType == 1:
            zeniths = tvals #+ np.max(tvals) - 2*np.min(tvals);
            azimuths = np.radians(np.linspace(0,360,np.size(xmesh)));
            r, theta = np.meshgrid(zeniths, azimuths);
            fig, axs = plt.subplots(1, 2, 
                subplot_kw=dict(projection='polar'),figsize=(16,6))
            c1 = axs[0].contourf(theta,r,uvals,50); 
            axs[0].set_xlabel('Radial Variable: Time [t]');
            axs[0].set_title('Excitatory Firing [u]\n'); plt.colorbar(c1,ax=axs[0]);
            c2 = axs[1].contourf(theta,r,vvals,50); 
            axs[1].set_xlabel('Radial Variable: Time [t]');
            axs[1].set_title('Inhibitory Firing [v]\n'); plt.colorbar(c2,ax=axs[1]);
            plt.show();
        elif pltType == 2:
            tt = np.append(tvals, tvals[-1]+dt)-dt/2; 
            rr = np.append(xmesh, xmesh[-1]+dx)-dx/2;
            TT, RR = np.meshgrid(tt, rr);
            XX = np.cos(2*np.pi*RR/self.L); YY = np.sin(2*np.pi*RR/self.L);
            norm_uvals = (uvals-np.min(uvals))/(np.max(uvals)-np.min(uvals));
            norm_vvals = (vvals-np.min(vvals))/(np.max(uvals)-np.min(uvals));
            cmap_uvals = cm.jet(norm_uvals)[:,:,0:3]; 
            cmap_vvals = cm.jet(norm_vvals)[:,:,0:3];
            ax1 = vv.subplot(211); ax1.daspect=(10/(np.max(tvals)-np.min(tvals)),1,1); 
            ax2 = vv.subplot(212); ax2.daspect=(10/(np.max(tvals)-np.min(tvals)),1,1); 
            vv.surf(TT, XX, YY, cmap_uvals,axes=ax1); vv.axis('tight',axes=ax1); 
            ax1.axis.xLabel = 'Time [t]'; vv.title('Excitatory Firing [u]',axes=ax1);
            vv.surf(TT, XX, YY, cmap_vvals,axes=ax2); vv.axis('tight',axes=ax2);
            ax2.axis.xLabel = 'Time [t]'; vv.title('Inhibitory Firing [v]',axes=ax2);
        elif pltType == 0:
            plt.subplot(2,1,1); # For Excitatory Firing Plot
            plt.imshow(uvals,interpolation='nearest',origin='lower',
                   extent=[np.min(tvals)-dt/2,np.max(tvals)+dt/2,
                           np.min(xmesh)-dx/2,np.max(xmesh)+dx/2],
                           cmap = cm.jet, aspect = 'auto');
            plt.title('Excitatory Firing [u]'); plt.colorbar();
            plt.xlabel('Time [t]'); plt.ylabel('Spatial Variable [x]'); 
            plt.subplot(2,1,2); # For Inhibitory FIring Plot
            plt.imshow(vvals,interpolation='nearest',origin='lower',
                   extent=[np.min(tvals)-dt/2,np.max(tvals)+dt/2,
                           np.min(xmesh)-dx/2,np.max(xmesh)+dx/2],
                           cmap = cm.jet, aspect = 'auto');
            plt.title('Inhibitory Firing [v]'); plt.colorbar();
            plt.xlabel('Time [t]'); plt.ylabel('Spatial Variable [x]'); 
            plt.subplots_adjust(hspace = 0.4); plt.show();
    def pltWilsonCowan1D(self,xmesh,yvals):
        uvals = yvals[0:self.nx]; vvals = yvals[self.nx:(2*self.nx)];
        plt.subplot(1,2,1); # For Firing Profiles Over Space
        plt.plot(xmesh,uvals[:,-1],'r',label='Excitatory Firing [u]');
        plt.plot(xmesh,vvals[:,-1],'b',label='Inhibitory Firing [v]');
        plt.xlabel('Spatial Variable [x]'); plt.ylabel('Firing Activity'); 
        plt.title('Spatial Firing Profiles'); plt.legend();
        plt.subplot(1,2,2); plt.plot(uvals[:,-1],vvals[:,-1],'k');
        plt.xlabel('Excitatory Firing [u]');
        plt.ylabel('Inhibitory Firing [v]');
        plt.title('Spatial Firing Phase Diagram'); 
        plt.subplots_adjust(wspace = 0.4); plt.show();
    def pltWilsonCowanVsT(self,tvals,yvals,index):
        uvals = yvals[0:self.nx]; vvals = yvals[self.nx:(2*self.nx)];
        plt.subplot(1,2,1); # For Firing Profiles Over Space
        plt.plot(tvals,uvals[index,:],'r',label='Excitatory Firing [u]');
        plt.plot(tvals,vvals[index,:],'b',label='Inhibitory Firing [v]');
        plt.xlabel('Time [t]'); plt.ylabel('Firing Activity'); 
        plt.title('Spatial Firing Profiles'); plt.legend();
        plt.subplot(1,2,2); plt.plot(uvals[index,:],vvals[index,:],'k');
        plt.xlabel('Excitatory Firing [u]');
        plt.ylabel('Inhibitory Firing [v]');
        plt.title('Temporal Firing Phase Diagram'); 
        plt.subplots_adjust(wspace = 0.4); plt.show(); 
        
    " Saving the results to file "
    def saveCurrentState(self,filename):
        writer = pd.ExcelWriter(filename)
        uvals = self.yvals[0:self.nx,:]; vvals = self.yvals[self.nx:(2*self.nx),:];
        uvalsDF = pd.DataFrame(np.transpose(uvals),columns=self.xmesh,index=self.tvals)
        vvalsDF = pd.DataFrame(np.transpose(vvals),columns=self.xmesh,index=self.tvals)
        tvalsDF = pd.DataFrame(self.tvals); xmeshDF = pd.DataFrame(self.xmesh);
        uvalsDF.to_excel(writer,'uvals'); vvalsDF.to_excel(writer,'vvals');
        tvalsDF.to_excel(writer,'tvals'); xmeshDF.to_excel(writer,'xmesh');
        paramsDF = pd.DataFrame({'BETA':self.BETA,'TE':self.TE,'TI':self.TI,
            'SE':self.SE,'SI':self.SI,'TAU':self.TAU,'L':self.L,'nx':self.nx,
            'AEE':self.AEE,'AIE':self.AIE,'AEI':self.AEI,'AII':self.AII,
            'span':self.span,'dx_kern':self.dx_kern,'mode':self.mode,
            'kernType':self.kernType,'dt':self.dt},index = ['Parameters']);
        paramsDF.to_excel(writer,'params');
    
    " Interactive real-time plotting by I/O dialog  "
    def interactiveIO_img(self,pltType=0):
        if self.WCSystem == None: 
            self.WilsonCowanIntegrator();
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs];
            self.imgWilsonCowan1D(tdisp, self.xmesh, ydisp, pltType); 
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        self.tmore = float(input("How much more integration time: ")); 
        self.tshow = float(input("How much time to display?: "));
        pltType = int(input("Plot Type (0-Normal, 1-Polar, 2-Cylinder): "));
        while(self.tshow>0):
            self.WCIntegrateLast(); locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs];
            self.imgWilsonCowan1D(tdisp, self.xmesh, ydisp, pltType); 
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename);
            self.tmore = float(input("How much more integration time: "));
            self.tshow = float(input("How much time to display?: "));
            pltType = int(input("Plot Type (0-Normal, 1-Polar, 2-Cylinder): "));
    def interactiveIO_pltX(self):
        if self.WCSystem == None: 
            self.WilsonCowanIntegrator();
            self.pltWilsonCowan1D(self.xmesh, self.yvals); 
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        self.tmore = float(input("How much more integration time: "));
        while(self.tmore>0):
            self.WCIntegrateLast();
            self.pltWilsonCowan1D(self.xmesh, self.yvals); 
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename);
            self.tmore = float(input("How much more integration time: "));  
    def interactiveIO_pltT(self,index=0):
        if self.WCSystem == None: 
            self.WilsonCowanIntegrator(); locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs];
            self.pltWilsonCowanVsT(tdisp, ydisp, index); 
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        self.tmore = float(input("How much more integration time: "));
        self.tshow = float(input("How much time to display?: "));
        index = int(input("Which neural element (0 to "+str(self.nx-1)+"): ")); 
        while(self.tshow>0):
            self.WCIntegrateLast(); locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs];
            self.pltWilsonCowanVsT(tdisp, ydisp, index);
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename);
            self.tmore = float(input("How much more integration time: "));
            self.tshow = float(input("How much time to display?: "));
            index = int(input("Which neural element (0 to "+str(self.nx-1)+"): "));    

    " Defining sigmoidal activation function "
    def f(self,x): return 1/(1+np.exp(-self.BETA*x))
    def df(self,x): return self.BETA*self.f(x)*(1-self.f(x))
    
    " Defining Convolution Kernel "
    def kern(self,sigma,span,dx): 
        x = dx*np.arange(-span,span+1);
        if self.kernType == 'Gaussian':
            return (1/sigma)*np.exp(-np.pi*(x**2)/(sigma**2))*dx;
        if self.kernType == 'Exponential':
            return (1/(2*sigma))*np.exp(-np.abs(x)/sigma)*dx;
        
    " Various Fourier Transforms of the above kernel "
    def kernFT(self,sig,L,nx,span,cycles):
        dx = L/(nx-1); fs = 1/dx; # Spacial Sampling
        numFTpts = 1000; numFTcycles = cycles; # For Fourier Transform of Spatial Mesh
        w = np.linspace(-numFTcycles*fs*np.pi,numFTcycles*fs*np.pi,numFTpts); 
        wnorm = w/fs; ker = self.kern(sig,span,dx); ww = np.array(wnorm)    
        WW, kerFT = signal.freqz(ker,worN=ww);
        return w, np.abs(kerFT), ker, dx
    def kernAnalyticFT(self,sig,w): 
        if self.kernType == 'Gaussian': return np.exp(-(sig**2)*(w**2)/(4*np.pi));
        if self.kernType == 'Exponential': return 1/(1+(sig*w)**2);
    def compareFTs(self,SI,L,nx,span,cycle):
        w, kerFTnumerical, ker, dx = self.kernFT(SI,L,nx,span,cycle);
        kerFTanalytical = self.kernAnalyticFT(SI,w);
        plt.plot(w,kerFTanalytical,'r', label = 'Continuous-Space Fourier Transform')
        plt.plot(w,kerFTnumerical,'b', label = 'Discrete-Space Fourier Transform');
        plt.xlabel('Spatial Frequency [Inverse Length]');
        plt.ylabel('Fourier Amplitude'); plt.legend(); plt.show();
    # Example Code for this (WC1D is a WilsonCowan1D object):
    # WC1D.compareFTs(0.1,2*np.pi,501,25,3)
    # WC1D.compareFTs(0.1,2*np.pi,400,60,2)
    
""" THE REST OF THIS FILE IS JUST A SCRIPT TO TEST THE CLASS-LEVEL
IMPLEMENTATIONS OF MY ALGORITHMS. ANYTHING BELOW THIS POINT ONLY 
SERVES AS EXAMPLE CODE AND NOTHING ELSE. EVERYTHING ABOVE THIS POINT 
IS SELF-SUFFICIENT CODE WITH MAY BE IMPORTED INTO ANOTHER PYTHON FILE"""
        
" FIX CLKWISE/CTRCLKWISE PROBLEM IN POLAR PLOTS "

" Zeroth: ANIMATE THESE PLOTS "      
        
" Second: Look at Harmonic Excitation based on Fourier Transforms "

" Third: Do what Bard did on wcsmooth200.ode for pattern tracking "

" Fourth: Make a GUI for this "


#" CONSTANTS: DO NOT CHANGE! "
#BETA = 50; TE = 0.125; TI = 0.4; SE = 12;
#AEE = 1;  AIE = 1; AEI = 1.5; AII = 0.25;
#
#" Defining ODE Parameters "
#L = 400; # Length of mesh
#nx = 400; # Spaces in mesh
#dx = L/(nx-1); # Spacing of mesh
#xx = np.linspace(0,L,nx); # the mesh itself
#span = 60; # Point to Left and Right of Central Maximum in Kernel
#dx_kern = 1; # Spacing of Points in the Kernel
#SI = 8; TAU = 0.6; # Remaining ODE Parameters
#dt = 0.025; tmore = 10; # Time Step and Initlal Integration Time
#tshow = 10; # Amount of Time Data to Display
#kernType='Gaussian' # 'Gaussian' or 'Exponential'
#mode = 'wrap'; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
#
#pardict = {'BETA':BETA, 'TE':TE, 'TI':TI, 'SE':SE, 'AEE': AEE, 'AIE':AIE,
#           'AEI':AEI, 'AII':AII, 'L':L, 'nx':nx, 'span':span, 'dx_kern':dx_kern,
#           'SI':SI, 'TAU': TAU, 'dt':dt, 'kernType':kernType, 'mode':mode};
#           
#WC1D = WilsonCowan1D(pardict=pardict); xx = WC1D.xmesh
#          
#" Initial Conditions for ODE "
## For Synchronous Solution:
##u0 = 0.4*np.ones((nx)); v0 = 0.2*np.ones((nx));
## For Travelling Wave in One Direction:
## SI = 4 and SE = 6 with travelling wave in one direction
## produces both a honeycomb and travelling wave
##u0 = 0.01*np.random.randn((nx))+0.1*np.sin(2*np.pi*xx/L)+0.4; 
##v0 = 0.01*np.random.randn((nx))+0.1*np.cos(2*np.pi*xx/L)+0.2;
## For Travelling Wave in Opposite Direction:
#u0 = 0.01*np.random.randn((nx))+0.1*np.cos(2*np.pi*xx/L)+0.4; 
#v0 = 0.01*np.random.randn((nx))+0.1*np.sin(2*np.pi*xx/L)+0.2;
## For Standing Wave (AKA: Honeycomb)
##u0 = 0.01*np.random.randn((nx))+0.41; 
##v0 = 0.01*np.random.randn((nx))+0.21;
#t0 = 0; WC1D.setInitConds(t0,u0,v0);
#WC1D.interactiveIO_img();
#WC1D2 = WilsonCowan1D(filename='test.xlsx');
#print("NOW REOPENING FILE INTO NEW OBJECT")
#WC1D2.interactiveIO_img();

    
    
    

