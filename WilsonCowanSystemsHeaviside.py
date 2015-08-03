# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 18:35:58 2015

@author: rehmanali
"""

from matplotlib.colors import colorConverter, Normalize
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage import convolve as conv
from scipy.signal import fftconvolve as fftconv
from scipy import signal
from scipy.integrate import ode
import matplotlib.pyplot as plt
from matplotlib import cm
import visvis as vv
import numpy as np
import pandas as pd


class WilsonCowan1D:
    'Organized Structure for 1D Wilson Cowan Equations Parameters and Results'
    tshow_default = 10; tmore_default = 10; twin_default = 5;
    
    " Object Constructor "
    def __init__(self,pardict=None,filename=None):
        if (filename is None) and (pardict is None):
            raise Exception('Must specify filename or pardict but not both')
        if filename is None:
            " Constant Params for ODE "
            self.TE = pardict['TE']; 
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
            paramFields = {'TE','TI','SE','SI','TAU','L','nx','dt',
                'AEE','AIE','AEI','AII','span','dx_kern','mode','kernType'}
            pardict = {'TE':None,'TI':None,'SE':None,'AEE': None, 
                       'AIE':None,'AEI':None,'AII':None, 'L':None, 'nx':None, 
                       'span':None,'dx_kern':None,'SI':None,'TAU':None,
                       'dt':None, 'kernType':None, 'mode':None};
            for names in paramFields:
                pardict[names] = dfs['params'][names]['Parameters'];
            self.TE = pardict['TE']; 
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
            self.WCSystem = ode(self.WilsonCowanODE).set_integrator('dopri5',nsteps=10000);
            self.WCSystem.set_initial_value(self.yvals[:,-1],self.tvals[-1]);
        else: raise Exception('Must specify filename or pardict but not both');
            
    " Retriving uvals and vvals after integration "
    def getVals(self):
        uvals = self.yvals[0:self.nx,:]; 
        vvals = self.yvals[self.nx:(2*self.nx),:];
        return uvals, vvals
    
    " Must Set Initial Conditions before Integration (WilsonCowanIntegrator) "    
    def setInitConds(self,t0,u0,v0,tmore=tmore_default,tshow=tshow_default,twin=twin_default):
        self.y0 = np.concatenate((u0,v0)); self.t0 = t0; 
        self.tmore = tmore; self.tshow = tshow; self.twin = twin
        
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
        self.WCSystem = ode(self.WilsonCowanODE).set_integrator('dopri5',nsteps=10000);
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

    " Plotting/Animating the Results of the Wilson Cowan Equations "
    def imgWilsonCowan1D(self, tvals, xmesh, yvals, pltType):
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
        plt.title('Spatial Firing Profiles'); plt.legend(framealpha=0);
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
        plt.title('Firing Profiles Over Time'); plt.legend(framealpha=0);
        plt.subplot(1,2,2); plt.plot(uvals[index,:],vvals[index,:],'k');
        plt.xlabel('Excitatory Firing [u]');
        plt.ylabel('Inhibitory Firing [v]');
        plt.title('Temporal Firing Phase Diagram'); 
        plt.subplots_adjust(wspace = 0.4); plt.show();
    def animWilsonCowan1D(self, tdisp, twin, ydisp, polar, tspeed=(4,1e-9)):
        dt = np.mean(np.diff(tdisp)); dx = np.mean(np.diff(self.xmesh));
        udisp = ydisp[0:self.nx,:]; vdisp = ydisp[self.nx:(2*self.nx),:]; 
        tshift, tpause = tspeed; tindex = 0; tsamples = np.floor(twin/dt)
        if polar:
            plt.ion(); plt.figure(figsize=(16,6))
            while (tindex + tsamples) < tdisp.size:
                plt.clf(); tvals = tdisp[tindex:(tindex+tsamples)];
                uvals = udisp[:,tindex:(tindex+tsamples)];
                vvals = vdisp[:,tindex:(tindex+tsamples)];
                zeniths = tvals #+ np.max(tvals) - 2*np.min(tvals);
                azimuths = np.radians(np.linspace(0,360,np.size(self.xmesh)));
                r, theta = np.meshgrid(zeniths, azimuths);
                ax0 = plt.subplot(121,polar=True);
                c1 = ax0.contourf(theta,r,uvals,100); 
                ax0.set_xlabel('Radial Variable: Time [t]');
                ax0.set_title('Excitatory Firing [u]\n'); plt.colorbar(c1,ax=ax0);
                ax1 = plt.subplot(122,polar=True);
                c2 = ax1.contourf(theta,r,vvals,100); 
                ax1.set_xlabel('Radial Variable: Time [t]');
                ax1.set_title('Inhibitory Firing [v]\n'); plt.colorbar(c2,ax=ax1);
                tindex += tshift; plt.draw(); plt.pause(tpause); 
        else:
            while (tindex + tsamples) < tdisp.size:
                plt.clf(); tvals = tdisp[tindex:(tindex+tsamples)];
                uvals = udisp[:,tindex:(tindex+tsamples)];
                vvals = vdisp[:,tindex:(tindex+tsamples)];
                plt.subplot(2,1,1); # For Excitatory Firing Plot
                plt.imshow(uvals,interpolation='nearest',origin='lower',
                       extent=[np.min(tvals)-dt/2,np.max(tvals)+dt/2,
                               np.min(self.xmesh)-dx/2,np.max(self.xmesh)+dx/2],
                               cmap = cm.jet, aspect = 'auto');
                plt.title('Excitatory Firing [u]'); plt.colorbar();
                plt.xlabel('Time [t]'); plt.ylabel('Spatial Variable [x]'); 
                plt.subplot(2,1,2); # For Inhibitory FIring Plot
                plt.imshow(vvals,interpolation='nearest',origin='lower',
                       extent=[np.min(tvals)-dt/2,np.max(tvals)+dt/2,
                               np.min(self.xmesh)-dx/2,np.max(self.xmesh)+dx/2],
                               cmap = cm.jet, aspect = 'auto');
                plt.title('Inhibitory Firing [v]'); plt.colorbar();
                plt.xlabel('Time [t]'); plt.ylabel('Spatial Variable [x]'); 
                plt.subplots_adjust(hspace = 0.4); 
                tindex += tshift; plt.draw(); plt.pause(tpause); 
    def animWilsonCowanVsX(self,tdisp,ydisp,tspeed=(4,1e-9)):
        tshift, tpause = tspeed; 
        shiftIndices = np.arange(np.floor(tdisp.size/tshift))
        for kk in shiftIndices:
            plt.clf(); tval = tdisp[kk*tshift]; yvals = ydisp[:,kk*tshift];
            uvals = yvals[0:self.nx]; vvals = yvals[self.nx:(2*self.nx)];
            plt.subplot(1,2,1); # For Firing Profiles Over Space
            plt.plot(self.xmesh,uvals[:],'r',label='Excitatory Firing [u]');
            plt.plot(self.xmesh,vvals[:],'b',label='Inhibitory Firing [v]');
            plt.xlabel('Spatial Variable [x]'); plt.ylabel('Firing Activity'); 
            plt.title('Spatial Firing Profiles at Time [t = '+str(tval)+']'); 
            plt.legend(framealpha=0); plt.subplot(1,2,2); plt.plot(uvals[:],vvals[:],'k');
            plt.xlabel('Excitatory Firing [u]'); plt.ylabel('Inhibitory Firing [v]');
            plt.title('Spatial Firing Phase Diagram at Time [t = '+str(tval)+']'); 
            plt.subplots_adjust(wspace = 0.4); 
            plt.draw(); plt.pause(tpause); 
    def animWilsonCowanVsT(self,tdisp,twin,ydisp,index,tspeed=(4,1e-9)):
        tshift, tpause = tspeed; 
        tindex = 0; tsamples = np.floor(twin/self.dt)
        udisp = ydisp[0:self.nx,:]; vdisp = ydisp[self.nx:(2*self.nx),:];
        while (tindex + tsamples) < tdisp.size:
            plt.clf(); tvals = tdisp[tindex:(tindex+tsamples)];
            uvals = udisp[:,tindex:(tindex+tsamples)];
            vvals = vdisp[:,tindex:(tindex+tsamples)];
            plt.subplot(1,2,1); # For Firing Profiles Over Space
            plt.plot(tvals,uvals[index,:],'r',label='Excitatory Firing [u]');
            plt.plot(tvals,vvals[index,:],'b',label='Inhibitory Firing [v]');
            plt.xlabel('Time [t]'); plt.ylabel('Firing Activity'); 
            plt.title('Firing Profiles Over Time'); plt.legend(framealpha=0);
            plt.subplot(1,2,2); plt.plot(uvals[index,:],vvals[index,:],'k');
            plt.xlabel('Excitatory Firing [u]');
            plt.ylabel('Inhibitory Firing [v]');
            plt.title('Temporal Firing Phase Diagram'); 
            plt.subplots_adjust(wspace = 0.4);
            tindex += tshift; plt.draw(); plt.pause(tpause); 
        
    " Saving the results to file "
    def saveCurrentState(self,filename):
        writer = pd.ExcelWriter(filename)
        uvals = self.yvals[0:self.nx,:]; vvals = self.yvals[self.nx:(2*self.nx),:];
        uvalsDF = pd.DataFrame(np.transpose(uvals),columns=self.xmesh,index=self.tvals)
        vvalsDF = pd.DataFrame(np.transpose(vvals),columns=self.xmesh,index=self.tvals)
        tvalsDF = pd.DataFrame(self.tvals); xmeshDF = pd.DataFrame(self.xmesh);
        uvalsDF.to_excel(writer,'uvals'); vvalsDF.to_excel(writer,'vvals');
        tvalsDF.to_excel(writer,'tvals'); xmeshDF.to_excel(writer,'xmesh');
        paramsDF = pd.DataFrame({'TE':self.TE,'TI':self.TI,
            'SE':self.SE,'SI':self.SI,'TAU':self.TAU,'L':self.L,'nx':self.nx,
            'AEE':self.AEE,'AIE':self.AIE,'AEI':self.AEI,'AII':self.AII,
            'span':self.span,'dx_kern':self.dx_kern,'mode':self.mode,
            'kernType':self.kernType,'dt':self.dt},index = ['Parameters']);
        paramsDF.to_excel(writer,'params');
        
    " Clearing all integration history and reinitialize to final values in last integration "
    def refreshFromLast(self,dt=None,tmore=tmore_default,tshow=tshow_default):
        if dt == None: dt = self.dt
        " Creating the WilsonCowan2D object "
        pardict = {'TE':self.TE,'TI':self.TI,'SE':self.SE, 
                   'AEE':self.AEE,'AIE':self.AIE,'AEI':self.AEI,'AII':self.AII, 
                   'L':self.L,'nx':self.nx,'span':self.span,
                   'dx_kern':self.dx_kern,'SI':self.SI,'TAU':self.TAU,
                   'dt':dt,'kernType':self.kernType,'mode':self.mode};
        WC1D = WilsonCowan1D(pardict=pardict); 
        uvals = self.yvals[0:(self.nx),:]; vvals = self.yvals[(self.nx):(2*self.nx),:]
        WC1D.setInitConds(0,uvals[:,-1],vvals[:,-1],tmore=tmore,tshow=tshow);
        print("Final Conditions from Saved Data made into Initial Conditions of New Simulation");
        return WC1D
    
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
    def interactiveIO_anim(self,polar=False,tshift=4):
        tpause = 1e-9; tspeed = (tshift,tpause);
        if self.WCSystem == None: 
            self.WilsonCowanIntegrator();
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs]; 
            self.animWilsonCowan1D(tdisp,self.twin,ydisp,polar,tspeed=tspeed)
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        self.tmore = float(input("How much more integration time: ")); 
        self.tshow = float(input("How much total time in animation?: "));
        self.twin = float(input("How much display time?: "));
        tshift = int(input("Time Shift Step for Animation (Integer >= 1): "));
        polar = bool(int(input("Plot Type (0-Normal, 1-Polar): ")));
        while(self.tshow>0):
            self.WCIntegrateLast(); locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs]; tspeed = (tshift,tpause);
            self.animWilsonCowan1D(tdisp,self.twin,ydisp,polar,tspeed=tspeed);
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename);
            self.tmore = float(input("How much more integration time: "));
            self.tshow = float(input("How much total time in animation?: "));
            self.twin = float(input("How much display time?: "));
            tshift = int(input("Time Shift Step for Animation (Integer >= 1): "));
            polar = bool(int(input("Plot Type (0-Normal, 1-Polar): ")));
    def interactiveIO_animX(self,tshift=4):
        tpause = 1e-9; tspeed = (tshift,tpause);        
        if self.WCSystem == None: 
            self.WilsonCowanIntegrator();
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs]; 
            self.animWilsonCowanVsX(tdisp,ydisp,tspeed=tspeed); 
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        self.tmore = float(input("How much more integration time: "));
        self.tshow = float(input("How much total time in animation?: "));
        tshift = int(input("Time Shift Step for Animation (Integer >= 1): "));
        while(self.tmore>0):
            self.WCIntegrateLast(); tspeed = (tshift,tpause);
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs]; 
            self.animWilsonCowanVsX(tdisp,ydisp,tspeed=tspeed); 
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename);
            self.tmore = float(input("How much more integration time: "));
            self.tshow = float(input("How much total time in animation?: "));
            tshift = int(input("Time Shift Step for Animation (Integer >= 1): "));
    def interactiveIO_animT(self,index=0,tshift=4):
        tpause = 1e-9; tspeed = (tshift,tpause);
        if self.WCSystem == None: 
            self.WilsonCowanIntegrator(); locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs];
            self.animWilsonCowanVsT(tdisp,self.twin,ydisp,index,tspeed); 
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        self.tmore = float(input("How much more integration time: "));
        self.tshow = float(input("How much total time in animation?: "));
        self.twin = float(input("How much display time?: "));
        tshift = int(input("Time Shift Step for Animation (Integer >= 1): "));
        index = int(input("Which neural element (0 to "+str(self.nx-1)+"): ")); 
        while(self.tshow>0):
            self.WCIntegrateLast(); locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs];
            self.animWilsonCowanVsT(tdisp,self.twin,ydisp,index,tspeed); 
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename);
            self.tmore = float(input("How much more integration time: "));
            self.tshow = float(input("How much total time in animation?: "));
            self.twin = float(input("How much display time?: "));
            tshift = int(input("Time Shift Step for Animation (Integer >= 1): "));
            index = int(input("Which neural element (0 to "+str(self.nx-1)+"): "));
            
    " Defining Heaviside or Unit Step Function "
    def heaviside(self,x): return 0.5 * (np.sign(x) + 1)
        
    " Defining sigmoidal activation function "
    def f(self,x): return self.heaviside(x) #1/(1+np.exp(-self.BETA*x))
    
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
        plt.ylabel('Fourier Amplitude'); plt.legend(framealpha=0); plt.show();
    # Example Code for this (WC1D is a WilsonCowan1D object):
    # WC1D.compareFTs(0.1,2*np.pi,501,25,3)
    # WC1D.compareFTs(0.1,2*np.pi,400,60,2)
    
class WilsonCowan2D:
    'Organized Structure for 1D Wilson Cowan Equations Parameters and Results'
    tshow_default = 10; tmore_default = 10; twin_default = 5;
    
    " Object Constructor "
    def __init__(self,pardict=None,filename=None):
        self.nthPowOfTwo = 2; # Radial Scaling of Visual Hallucinogram
        if (filename is None) and (pardict is None):
            raise Exception('Must specify filename or pardict but not both')
        if filename is None:
            " Constant Params for ODE "
            self.TE = pardict['TE']; 
            self.TI = pardict['TI']; self.SE = pardict['SE'];
            self.AEE = pardict['AEE']; self.AIE = pardict['AIE']; 
            self.AEI = pardict['AEI']; self.AII = pardict['AII'];
    
            " Defining Other ODE Parameters "
            self.Lx = pardict['Lx']; self.Ly = pardict['Ly'] # Length of mesh
            self.nx = pardict['nx']; self.ny = pardict['ny'] # Spaces in mesh
            self.dx = self.Lx/(self.nx-1); self.dy = self.Ly/(self.ny-1) # Spacing of mesh
            self.xmesh = np.linspace(0,self.Lx,self.nx); 
            self.ymesh = np.linspace(0,self.Ly,self.ny);
            self.XMesh, self.YMesh = np.meshgrid(self.xmesh,self.ymesh);
            self.span_x = pardict['span_x']; self.span_y = pardict['span_y']; 
            self.dx_kern = pardict['dx_kern']; self.dy_kern = pardict['dy_kern'];
            self.SI = pardict['SI']; self.TAU = pardict['TAU']; # Remaining ODE Parameters
            self.dt = pardict['dt']; self.WCSystem = None;
            self.kernType= pardict['kernType']; # 'Gaussian' or 'Exponential'
            self.mode = pardict['mode']; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
        elif pardict is None:
            " Retrieving saved state from file "
            xl = pd.ExcelFile(filename);
            dfs = {sheet: xl.parse(sheet) for sheet in xl.sheet_names}
            paramFields = {'TE','TI','SE','SI','TAU','Lx','Ly','nx','ny',
                'AEE','AIE','AEI','AII','span_x','span_y','dx_kern','dy_kern',
                'dt','mode','kernType'}
            pardict = {'TE':None,'TI':None,'SE':None,'AEE': None, 
                       'AIE':None,'AEI':None,'AII':None, 'Lx':None, 'Ly':None,
                       'nx':None,'ny':None,'span_x':None,'span_y':None,
                       'dx_kern':None,'dy_kern':None,'SI':None,'TAU':None,
                       'dt':None, 'kernType':None, 'mode':None};
            for names in paramFields:
                pardict[names] = dfs['params'][names]['Parameters'];
            self.TE = pardict['TE']; 
            self.TI = pardict['TI']; self.SE = pardict['SE'];
            self.AEE = pardict['AEE']; self.AIE = pardict['AIE']; 
            self.AEI = pardict['AEI']; self.AII = pardict['AII'];
            self.Lx = pardict['Lx']; self.Ly = pardict['Ly'] # Length of mesh
            self.nx = pardict['nx']; self.ny = pardict['ny']# Spaces in mesh
            self.dx = self.Lx/(self.nx-1); self.dy = self.Ly/(self.ny-1); # Spacing of mesh
            self.xmesh = np.linspace(0,self.Lx,self.nx); 
            self.ymesh = np.linspace(0,self.Ly,self.ny);
            self.XMesh, self.YMesh = np.meshgrid(self.xmesh,self.ymesh);
            self.span_x = pardict['span_x']; self.span_y = pardict['span_y'];
            self.dx_kern = pardict['dx_kern']; self.dy_kern = pardict['dy_kern'];
            self.SI = pardict['SI']; self.TAU = pardict['TAU']; # Remaining ODE Parameters
            self.dt = pardict['dt']; self.WCSystem = None;
            self.kernType= pardict['kernType']; # 'Gaussian' or 'Exponential'
            self.mode = pardict['mode']; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
            self.xmesh = dfs['xmesh'][0].values; self.ymesh = dfs['ymesh'][0].values;
            self.XMesh, self.YMesh = np.meshgrid(self.xmesh,self.ymesh);
            self.tvals = dfs['tvals'][0].values;
            uvals = dfs['uvals'].values; vvals = dfs['vvals'].values;
            uvals = uvals.reshape((self.tvals.size,self.nx*self.ny));
            vvals = vvals.reshape((self.tvals.size,self.nx*self.ny));
            self.yvals = np.concatenate((np.transpose(uvals),np.transpose(vvals)));
            self.WCSystem = ode(self.WilsonCowanODE).set_integrator('dopri5',nsteps=10000);
            self.WCSystem.set_initial_value(self.yvals[:,-1],self.tvals[-1]);
        else: raise Exception('Must specify filename or pardict but not both');
            
    " Retriving uvals and vvals after integration "
    def getVals(self):
        uvals = self.yvals[0:(self.nx*self.ny),:]; 
        vvals = self.yvals[(self.nx*self.ny):(2*self.nx*self.ny),:]
        uvals = np.reshape(uvals,(self.ny,self.nx,self.tvals.size)); 
        vvals = np.reshape(vvals,(self.ny,self.nx,self.tvals.size));
        return uvals, vvals    
    
    " Must Set Initial Conditions before Integration (WilsonCowanIntegrator) "    
    def setInitConds(self,t0,u0,v0,tmore=tmore_default,tshow=tshow_default,twin=twin_default):
        uu0 = np.reshape(u0,self.nx*self.ny); vv0 = np.reshape(v0,self.nx*self.ny)
        self.yy0 = np.concatenate((uu0,vv0)); self.t0 = t0; self.twin = twin
        self.tmore = tmore; self.tshow = tshow;
    
    " Defining the Wilson Cowan Equations "
    def WilsonCowanODE(self,t,y):
        u = y[0:(self.nx*self.ny)]; v = y[(self.nx*self.ny):(2*self.nx*self.ny)];
        U = np.reshape(u,(self.ny,self.nx)); V = np.reshape(v,(self.ny,self.nx));
        kerE = self.kern(self.SE,(self.span_x,self.span_y),(self.dx_kern,self.dy_kern)); 
        kerI = self.kern(self.SI,(self.span_x,self.span_y),(self.dx_kern,self.dy_kern)); 
        Lu = self.AEE*conv(U,kerE,mode=self.mode)-self.AEI*conv(V,kerI,mode=self.mode)-self.TE;
        Lv = self.AIE*conv(U,kerE,mode=self.mode)-self.AII*conv(V,kerI,mode=self.mode)-self.TI;
        dU = -U + self.f(Lu); dV = (-V + self.f(Lv))/self.TAU;
        du = np.reshape(dU,self.ny*self.nx); dv = np.reshape(dV,self.ny*self.nx);
        return np.concatenate((du,dv));
    def WilsonCowanIntegrator(self):
        self.WCSystem = ode(self.WilsonCowanODE).set_integrator('dopri5',nsteps=10000);
        self.WCSystem.set_initial_value(self.yy0,self.t0);
        self.tvals, self.yvals = [self.t0], self.yy0;
        while self.WCSystem.successful() and self.WCSystem.t < self.tmore:
            self.WCSystem.integrate(self.WCSystem.t+self.dt);
            self.tvals = np.concatenate((self.tvals,[self.WCSystem.t]));
            self.yvals = np.column_stack((self.yvals,self.WCSystem.y));
        return self.tvals, self.yvals, self.WCSystem
    def WCIntegrateLast(self):
        tlast = self.tvals[-1]; 
        while self.WCSystem.successful() and self.WCSystem.t < (tlast+self.tmore):
            self.WCSystem.integrate(self.WCSystem.t+self.dt);
            self.tvals = np.concatenate((self.tvals,[self.WCSystem.t]));
            self.yvals = np.column_stack((self.yvals,self.WCSystem.y));
        return self.tvals, self.yvals, self.WCSystem
    
    " Plotting/Animating the Results of the Wilson Cowan Equations "
    def imgWilsonCowan3D(self,tvals,yvals,video):
        uvals = yvals[0:(self.nx*self.ny),:]; 
        vvals = yvals[(self.nx*self.ny):(2*self.nx*self.ny),:]
        uvals = np.reshape(uvals,(self.ny,self.nx,tvals.size)); 
        vvals = np.reshape(vvals,(self.ny,self.nx,tvals.size)); vv.clf();
        xrange = (np.min(self.xmesh),np.max(self.xmesh));
        yrange = (np.min(self.ymesh),np.max(self.ymesh));
        trange = (np.min(tvals),np.max(tvals));    
        if video:
            uval_list = [kk for kk in np.swapaxes(np.swapaxes(uvals,0,2),1,2)];
            vval_list = [kk for kk in np.swapaxes(np.swapaxes(vvals,0,2),1,2)];
            speedFactor = float(input('Speed-Up Factor: '));
            vv.subplot(121); ax1 = vv.gca(); 
            vv.movieShow(uval_list,duration=self.dt/speedFactor,axesAdjust=True,axes=ax1);
            vv.colorbar(axes=ax1); ax1.daspect = (1,-1,1); ax1.axis.tickFontSize = 0;
            vv.xlabel('Space [x='+str(xrange[0])+'..'+str(xrange[1])+']',axes=ax1); 
            vv.ylabel('Space [y='+str(yrange[0])+'..'+str(yrange[1])+']',axes=ax1);
            vv.title('Excitatory Firing [u]',axes=ax1); vv.axis('tight',axes=ax1);
            vv.subplot(122); ax2 = vv.gca(); 
            vv.movieShow(vval_list,duration=self.dt/speedFactor,axesAdjust=True,axes=ax2);
            vv.colorbar(axes=ax2); ax2.daspect = (1,-1,1); ax2.axis.tickFontSize = 0;
            vv.xlabel('Space [x='+str(xrange[0])+'..'+str(xrange[1])+']',axes=ax2); 
            vv.ylabel('Space [y='+str(yrange[0])+'..'+str(yrange[1])+']',axes=ax2);
            vv.title('Inhibitory Firing [v]',axes=ax2); vv.axis('tight',axes=ax2); 
        else:
            ax1 = vv.subplot(121); ax2 = vv.subplot(122);
            ax1.axis.tickFontSize = 0; ax2.axis.tickFontSize = 0;
            ax1.daspect=(2,1,1); ax2.daspect=(2,1,1); 
            vv.ylabel('Space [x='+str(xrange[0])+'..'+str(xrange[1])+']',axes=ax1); 
            vv.zlabel('Space [y='+str(yrange[0])+'..'+str(yrange[1])+']',axes=ax1);
            vv.xlabel('Time [t='+str(trange[0])+'..'+str(trange[1])+']',axes=ax1); 
            vv.title('Excitatory Firing [u]',axes=ax1); vv.axis('tight',axes=ax1);
            vv.ylabel('Space [x='+str(xrange[0])+'..'+str(xrange[1])+']',axes=ax2); 
            vv.zlabel('Space [y='+str(yrange[0])+'..'+str(yrange[1])+']',axes=ax2);
            vv.xlabel('Time [t='+str(trange[0])+'..'+str(trange[1])+']',axes=ax2);
            vv.title('Inhibitory Firing [v]',axes=ax2); vv.axis('tight',axes=ax2); 
            vv.volshow(uvals,axes=ax1); vv.volshow(vvals,axes=ax2);
            vv.ColormapEditor(ax1); vv.ColormapEditor(ax2);
    def imgWilsonCowan2D(self, tval, yvals, pltType):
        uvals = yvals[0:(self.nx*self.ny)]; 
        vvals = yvals[(self.nx*self.ny):(2*self.nx*self.ny)];
        uvals = np.reshape(uvals,(self.ny,self.nx)); 
        vvals = np.reshape(vvals,(self.ny,self.nx));
        if pltType == 2: 
            R = 5; r = 2; # Major and Minor Radii For Torus Geometry
            ttx = np.append(self.xmesh,self.xmesh[-1]+self.dx)-self.dx/2; 
            tty = np.append(self.ymesh, self.ymesh[-1]+self.dy)-self.dy/2;
            TTx, TTy = np.meshgrid(ttx, tty);
            rangeX = self.Lx + self.dx; rangeY = self.Ly+self.dy
            XX1 = (R+r*np.cos(2*np.pi*TTx/rangeX))*np.cos(2*np.pi*TTy/rangeY); 
            YY1 = (R+r*np.cos(2*np.pi*TTx/rangeX))*np.sin(2*np.pi*TTy/rangeY);
            ZZ1 = r*np.sin(2*np.pi*TTx/rangeX);
            XX2 = (R+r*np.cos(2*np.pi*TTy/rangeY))*np.cos(2*np.pi*TTx/rangeX); 
            YY2 = (R+r*np.cos(2*np.pi*TTy/rangeY))*np.sin(2*np.pi*TTx/rangeX);
            ZZ2 = r*np.sin(2*np.pi*TTy/rangeY);
            norm_uvals = (uvals-np.min(uvals))/(np.max(uvals)-np.min(uvals));
            norm_vvals = (vvals-np.min(vvals))/(np.max(uvals)-np.min(uvals));
            cmap_uvals = cm.jet(norm_uvals)[:,:,0:3]; 
            cmap_vvals = cm.jet(norm_vvals)[:,:,0:3];
            ax1 = vv.subplot(221); #ax1.daspect=(10/(np.max(tvals)-np.min(tvals)),1,1); 
            ax2 = vv.subplot(222); #ax2.daspect=(10/(np.max(tvals)-np.min(tvals)),1,1);
            ax3 = vv.subplot(223); #ax1.daspect=(10/(np.max(tvals)-np.min(tvals)),1,1); 
            ax4 = vv.subplot(224); #ax2.daspect=(10/(np.max(tvals)-np.min(tvals)),1,1);
            vv.surf(XX1, YY1, ZZ1, cmap_uvals,axes=ax1); vv.axis('tight',axes=ax1);
            vv.title('Excitatory Firing [u] Over Torus at Time [t = '+str(tval)+']',axes=ax1);
            vv.surf(XX1, YY1, ZZ1, cmap_vvals,axes=ax2); vv.axis('tight',axes=ax2); 
            vv.title('Inhibitory Firing [v] Over Torus Time [t = '+str(tval)+']',axes=ax2);
            vv.surf(XX2, YY2, ZZ2, cmap_uvals,axes=ax3); vv.axis('tight',axes=ax3);
            vv.title('Excitatory Firing [u] Over Ortho-Torus Time [t = '+str(tval)+']',axes=ax3);
            vv.surf(XX2, YY2, ZZ2, cmap_vvals,axes=ax4); vv.axis('tight',axes=ax4);
            vv.title('Inhibitory Firing [v] Over Ortho-Torus Time [t = '+str(tval)+']',axes=ax4);
        elif pltType == 0:
            plt.subplot(1,2,1); # For Excitatory Firing Plot
            plt.imshow(uvals,interpolation='nearest',origin='lower',
                   extent=[np.min(self.xmesh)-self.dx/2,np.max(self.xmesh)+self.dx/2,
                           np.min(self.ymesh)-self.dy/2,np.max(self.ymesh)+self.dy/2],
                           cmap = cm.jet, aspect = 'auto'); plt.colorbar();
            plt.title('Excitatory Firing [u] at Time [t = '+str(tval)+']'); 
            plt.xlabel('Spatial Variable [x]'); plt.ylabel('Spatial Variable [y]'); 
            plt.subplot(1,2,2); # For Inhibitory Firing Plot
            plt.imshow(vvals,interpolation='nearest',origin='lower',
                   extent=[np.min(self.xmesh)-self.dx/2,np.max(self.xmesh)+self.dx/2,
                           np.min(self.ymesh)-self.dy/2,np.max(self.ymesh)+self.dy/2],
                           cmap = cm.jet, aspect = 'auto'); plt.colorbar();
            plt.title('Inhibitory Firing [v] at Time [t = '+str(tval)+']'); 
            plt.xlabel('Spatial Variable [x]'); plt.ylabel('Spatial Variable [y]'); 
            plt.subplots_adjust(wspace = 0.4); plt.show();
        elif pltType == 1:
            eps_0 = 1; alpha = (self.Lx+self.dx)/(self.nthPowOfTwo*np.log(2));
            XMesh_disp = np.concatenate((self.XMesh,self.XMesh[-1,:][np.newaxis,:]),axis=0);
            eps = eps_0*np.exp(XMesh_disp/alpha); 
            a = 2*np.pi*self.YMesh/(self.Ly+self.dy);
            a = np.vstack((a,2*np.pi*np.ones(a.shape[1])));
            uvals_disp = np.concatenate((uvals,uvals[0,:][np.newaxis,:]),axis=0);
            vvals_disp = np.concatenate((vvals,vvals[0,:][np.newaxis,:]),axis=0);
            fig, axs = plt.subplots(1, 2, 
                subplot_kw=dict(projection='polar'),figsize=(16,7))
            c1 = axs[0].contourf(a,eps,uvals_disp,100); 
            axs[0].set_title('Visual Hallucination Due to \nExcitatory Firing [u]\nat Time [t = '+str(tval)+']\n'); 
            plt.colorbar(c1,ax=axs[0]);
            c2 = axs[1].contourf(a,eps,vvals_disp,100); 
            axs[1].set_title('Visual Hallucination Due to \nInhibitory Firing [v]\nat Time [t = '+str(tval)+']\n'); 
            plt.colorbar(c2,ax=axs[1]); plt.show();
    def imgWilsonCowan1DvsT(self,tvals,yvals,info):
        axis, indx, pltType = info; 
        if axis == 1: 
            mesh = self.ymesh; diff = self.dy;
            uvals = yvals[0:self.ny,:]; 
            vvals = yvals[self.ny:(2*self.ny),:];
        elif axis == 0: 
            mesh = self.xmesh; diff = self.dx;
            uvals = yvals[0:self.nx,:]; 
            vvals = yvals[self.nx:(2*self.nx),:];
        if pltType == 0: # Normal Plot
            plt.subplot(2,1,1); # For Excitatory Firing Profiles Over Space
            plt.imshow(uvals,interpolation='nearest',origin='lower',
                   extent=[np.min(tvals)-self.dt/2,np.max(tvals)+self.dt/2,
                           np.min(mesh)-diff/2,np.max(mesh)+diff/2],
                           cmap = cm.jet, aspect = 'auto'); 
            plt.colorbar(); plt.xlabel('Time [t]'); 
            if axis == 1: 
                plt.ylabel('Spatial Variable [y]'); 
                plt.title('Excitatory Firing [u] at [x = '+str(indx)+']'); 
            elif axis == 0: 
                plt.ylabel('Spatial Variable [x]');
                plt.title('Excitatory Firing [u] at [y = '+str(indx)+']');
            plt.subplot(2,1,2); # For Inhibitory Firing Profiles Over Space
            plt.imshow(vvals,interpolation='nearest',origin='lower',
                   extent=[np.min(tvals)-self.dt/2,np.max(tvals)+self.dt/2,
                           np.min(mesh)-diff/2,np.max(mesh)+diff/2],
                           cmap = cm.jet, aspect = 'auto'); 
            plt.colorbar(); plt.xlabel('Time [t]'); 
            if axis == 1: 
                plt.ylabel('Spatial Variable [y]'); 
                plt.title('Inhibitory Firing [v] at [x = '+str(indx)+']'); 
            elif axis == 0: 
                plt.ylabel('Spatial Variable [x]');
                plt.title('Inhibitory Firing [v] at [y = '+str(indx)+']'); 
            plt.subplots_adjust(hspace = 0.4); plt.show();
        if pltType == 1: # Polar Plot
            zeniths = tvals #+ np.max(tvals) - 2*np.min(tvals);
            azimuths = np.radians(np.linspace(0,360,np.size(mesh)));
            r, theta = np.meshgrid(zeniths, azimuths);
            fig, axs = plt.subplots(1, 2, 
                subplot_kw=dict(projection='polar'),figsize=(16,6))
            c1 = axs[0].contourf(theta,r,uvals,50); 
            if axis == 1:
                axs[0].set_xlabel('Radial Variable: Time [t]\nAngular Variable [y]');
                axs[0].set_title('Excitatory Firing [u] at [x = '+str(indx)+']\n');
            elif axis == 0:
                axs[0].set_xlabel('Radial Variable: Time [t]\nAngular Variable [x]');
                axs[0].set_title('Excitatory Firing [u] at [y = '+str(indx)+']\n');
            plt.colorbar(c1,ax=axs[0]);
            c2 = axs[1].contourf(theta,r,vvals,50); 
            if axis == 1:
                axs[0].set_xlabel('Radial Variable: Time [t]\nAngular Variable [y]');
                axs[0].set_title('Inhibitory Firing [v] at [x = '+str(indx)+']\n');
            elif axis == 0:
                axs[1].set_xlabel('Radial Variable: Time [t]\nAngular Variable [x]');
                axs[1].set_title('Inhibitory Firing [v] at [y = '+str(indx)+']\n');
            plt.colorbar(c2,ax=axs[1]); plt.show();
        if pltType == 2: # Cylinder Plot
            tt = np.append(tvals, tvals[-1]+self.dt)-self.dt/2; 
            rr = np.append(mesh, mesh[-1]+diff)-diff/2;
            TT, RR = np.meshgrid(tt, rr);
            XX = np.cos(2*np.pi*RR/self.Lx); YY = np.sin(2*np.pi*RR/self.Ly);
            norm_uvals = (uvals-np.min(uvals))/(np.max(uvals)-np.min(uvals));
            norm_vvals = (vvals-np.min(vvals))/(np.max(uvals)-np.min(uvals));
            cmap_uvals = cm.jet(norm_uvals)[:,:,0:3]; 
            cmap_vvals = cm.jet(norm_vvals)[:,:,0:3];
            ax1 = vv.subplot(211); ax2 = vv.subplot(212); 
            ax1.daspect=(10/(np.max(tvals)-np.min(tvals)),1,1); 
            ax2.daspect=(10/(np.max(tvals)-np.min(tvals)),1,1); 
            vv.surf(TT, XX, YY, cmap_uvals,axes=ax1); 
            vv.axis('tight',axes=ax1); ax1.axis.xLabel = 'Time [t]'; 
            vv.surf(TT, XX, YY, cmap_vvals,axes=ax2); 
            vv.axis('tight',axes=ax2); ax2.axis.xLabel = 'Time [t]';
            if axis == 1:
                vv.title('Excitatory Firing [u] Over [y] at [x = '+str(indx)+']',axes=ax1);
                vv.title('Inhibitory Firing [v] Over [y] at [x = '+str(indx)+']',axes=ax2);
            elif axis == 0:
                vv.title('Excitatory Firing [u] Over [x] at [y = '+str(indx)+']',axes=ax1);
                vv.title('Inhibitory Firing [v] Over [x] at [y = '+str(indx)+']',axes=ax2);
    def pltWilsonCowan1DvsT(self,tvals,uvals,vvals):
        plt.subplot(1,2,1); # For Firing Profiles Over Space
        plt.plot(tvals,uvals,'r',label='Excitatory Firing [u]');
        plt.plot(tvals,vvals,'b',label='Inhibitory Firing [v]');
        plt.xlabel('Time [t]'); plt.ylabel('Firing Activity'); 
        plt.title('Firing Profiles Over Time'); plt.legend(framealpha=0);
        plt.subplot(1,2,2); plt.plot(uvals,vvals,'k');
        plt.xlabel('Excitatory Firing [u]');
        plt.ylabel('Inhibitory Firing [v]');
        plt.title('Temporal Firing Phase Diagram'); 
        plt.subplots_adjust(wspace = 0.4); plt.show();
    def animWilsonCowan2D(self, tdisp, ydisp, polar, tspeed=(4,1e-9)):
        udisp = ydisp[0:(self.nx*self.ny),:]; 
        vdisp = ydisp[(self.nx*self.ny):(2*self.nx*self.ny),:];
        udisp = np.reshape(udisp,(self.ny,self.nx,tdisp.size)); 
        vdisp = np.reshape(vdisp,(self.ny,self.nx,tdisp.size));
        tshift, tpause = tspeed; tindex = 0;
        if polar:
            eps_0 = 1; alpha = (self.Lx+self.dx)/(self.nthPowOfTwo*np.log(2));
            XMesh_disp = np.concatenate((self.XMesh,self.XMesh[-1,:][np.newaxis,:]),axis=0);
            eps = eps_0*np.exp(XMesh_disp/alpha); 
            a = 2*np.pi*self.YMesh/(self.Ly+self.dy);
            a = np.vstack((a,2*np.pi*np.ones(a.shape[1])));
            plt.ion(); plt.figure(figsize=(16,7))
            while tindex < tdisp.size:
                plt.clf(); uvals = udisp[:,:,tindex]; vvals = vdisp[:,:,tindex]; 
                uvals_disp = np.concatenate((uvals,uvals[0,:][np.newaxis,:]),axis=0);
                vvals_disp = np.concatenate((vvals,vvals[0,:][np.newaxis,:]),axis=0);
                ax0 = plt.subplot(121,polar=True); ax1 = plt.subplot(122,polar=True);
                c1 = ax0.contourf(a,eps,uvals_disp,100); tval = tdisp[tindex];
                ax0.set_title('Visual Hallucination Due to \nExcitatory Firing [u]\nat Time [t = '+str(tval)+']\n'); 
                plt.colorbar(c1,ax=ax0); c2 = ax1.contourf(a,eps,vvals_disp,100); 
                ax1.set_title('Visual Hallucination Due to \nInhibitory Firing [v]\nat Time [t = '+str(tval)+']\n'); 
                plt.colorbar(c2,ax=ax1); plt.draw(); plt.pause(tpause); tindex += tshift;
        else:
            while tindex < tdisp.size:
                plt.clf(); tval = tdisp[tindex]; plt.subplot(1,2,1); # For Excitatory Firing Plot
                plt.imshow(udisp[:,:,tindex],interpolation='nearest',origin='lower',
                       extent=[np.min(self.xmesh)-self.dx/2,np.max(self.xmesh)+self.dx/2,
                               np.min(self.ymesh)-self.dy/2,np.max(self.ymesh)+self.dy/2],
                               cmap = cm.jet, aspect = 'auto'); plt.colorbar();
                plt.title('Excitatory Firing [u] at Time [t = '+str(tval)+']'); 
                plt.xlabel('Spatial Variable [x]'); plt.ylabel('Spatial Variable [y]'); 
                plt.subplot(1,2,2); # For Inhibitory Firing Plot
                plt.imshow(vdisp[:,:,tindex],interpolation='nearest',origin='lower',
                       extent=[np.min(self.xmesh)-self.dx/2,np.max(self.xmesh)+self.dx/2,
                               np.min(self.ymesh)-self.dy/2,np.max(self.ymesh)+self.dy/2],
                               cmap = cm.jet, aspect = 'auto'); plt.colorbar();
                plt.title('Inhibitory Firing [v] at Time [t = '+str(tval)+']'); 
                plt.xlabel('Spatial Variable [x]'); plt.ylabel('Spatial Variable [y]'); 
                plt.subplots_adjust(wspace = 0.4); plt.draw(); plt.pause(tpause); tindex += tshift;
    def animWilsonCowan1DvsT(self,tdisp,twin,ydisp,axis,indx,polar,tspeed=(4,1e-9)):
        tshift, tpause = tspeed; tindex = 0; tsamples = np.floor(twin/self.dt)
        if axis == 1: 
            mesh = self.ymesh; diff = self.dy;
            udisp = ydisp[0:self.ny,:]; 
            vdisp = ydisp[self.ny:(2*self.ny),:];
        elif axis == 0: 
            mesh = self.xmesh; diff = self.dx;
            udisp = ydisp[0:self.nx,:]; 
            vdisp = ydisp[self.nx:(2*self.nx),:];
        if polar: # Polar Plot
            plt.ion(); plt.figure(figsize=(16,7))
            while (tindex + tsamples) < tdisp.size:
                plt.clf(); zeniths = tdisp[tindex:(tindex+tsamples)];
                uvals = udisp[:,tindex:(tindex+tsamples)];
                vvals = vdisp[:,tindex:(tindex+tsamples)];
                azimuths = np.radians(np.linspace(0,360,np.size(mesh)));
                r, theta = np.meshgrid(zeniths, azimuths); 
                ax0 = plt.subplot(121,polar=True); ax1 = plt.subplot(122,polar=True);
                c1 = ax0.contourf(theta,r,uvals,50); 
                if axis == 1:
                    ax0.set_xlabel('Radial Variable: Time [t]\nAngular Variable [y]');
                    ax0.set_title('Excitatory Firing [u] at [x = '+str(indx)+']\n');
                elif axis == 0:
                    ax0.set_xlabel('Radial Variable: Time [t]\nAngular Variable [x]');
                    ax0.set_title('Excitatory Firing [u] at [y = '+str(indx)+']\n');
                plt.colorbar(c1,ax=ax0);
                c2 = ax1.contourf(theta,r,vvals,50); 
                if axis == 1:
                    ax0.set_xlabel('Radial Variable: Time [t]\nAngular Variable [y]');
                    ax0.set_title('Inhibitory Firing [v] at [x = '+str(indx)+']\n');
                elif axis == 0:
                    ax1.set_xlabel('Radial Variable: Time [t]\nAngular Variable [x]');
                    ax1.set_title('Inhibitory Firing [v] at [y = '+str(indx)+']\n');
                plt.colorbar(c2,ax=ax1); plt.draw(); plt.pause(tpause); tindex += tshift;
        else: # Normal Plot
            while (tindex + tsamples) < tdisp.size:
                plt.clf(); tvals = tdisp[tindex:(tindex+tsamples)];
                uvals = udisp[:,tindex:(tindex+tsamples)];
                vvals = vdisp[:,tindex:(tindex+tsamples)];
                plt.subplot(2,1,1); # For Excitatory Firing Profiles Over Space
                plt.imshow(uvals,interpolation='nearest',origin='lower',
                       extent=[np.min(tvals)-self.dt/2,np.max(tvals)+self.dt/2,
                               np.min(mesh)-diff/2,np.max(mesh)+diff/2],
                               cmap = cm.jet, aspect = 'auto'); 
                plt.colorbar(); plt.xlabel('Time [t]'); 
                if axis == 1: 
                    plt.ylabel('Spatial Variable [y]'); 
                    plt.title('Excitatory Firing [u] at [x = '+str(indx)+']'); 
                elif axis == 0: 
                    plt.ylabel('Spatial Variable [x]');
                    plt.title('Excitatory Firing [u] at [y = '+str(indx)+']');
                plt.subplot(2,1,2); # For Inhibitory Firing Profiles Over Space
                plt.imshow(vvals,interpolation='nearest',origin='lower',
                       extent=[np.min(tvals)-self.dt/2,np.max(tvals)+self.dt/2,
                               np.min(mesh)-diff/2,np.max(mesh)+diff/2],
                               cmap = cm.jet, aspect = 'auto'); 
                plt.colorbar(); plt.xlabel('Time [t]'); 
                if axis == 1: 
                    plt.ylabel('Spatial Variable [y]'); 
                    plt.title('Inhibitory Firing [v] at [x = '+str(indx)+']'); 
                elif axis == 0: 
                    plt.ylabel('Spatial Variable [x]');
                    plt.title('Inhibitory Firing [v] at [y = '+str(indx)+']'); 
                plt.subplots_adjust(hspace = 0.4); plt.draw(); 
                plt.pause(tpause); tindex += tshift;
    def animWilsonCowanPltT(self,tdisp,twin,udisp,vdisp,tspeed=(4,1e-9)):
        tshift, tpause = tspeed; tindex = 0; tsamples = np.floor(twin/self.dt);
        while (tindex + tsamples) < tdisp.size:
            plt.clf(); tvals = tdisp[tindex:(tindex+tsamples)];
            uvals = udisp[tindex:(tindex+tsamples)];
            vvals = vdisp[tindex:(tindex+tsamples)];
            plt.subplot(1,2,1); # For Firing Profiles Over Space
            plt.plot(tvals,uvals,'r',label='Excitatory Firing [u]');
            plt.plot(tvals,vvals,'b',label='Inhibitory Firing [v]');
            plt.xlabel('Time [t]'); plt.ylabel('Firing Activity'); 
            plt.title('Firing Profiles Over Time'); plt.legend(framealpha=0);
            plt.subplot(1,2,2); plt.plot(uvals,vvals,'k');
            plt.xlabel('Excitatory Firing [u]');
            plt.ylabel('Inhibitory Firing [v]');
            plt.title('Temporal Firing Phase Diagram'); 
            plt.subplots_adjust(wspace = 0.4); plt.draw(); 
            plt.pause(tpause); tindex += tshift;
    
    " Saving the results to file "
    def saveCurrentState(self,filename):
        writer = pd.ExcelWriter(filename)
        uvals = self.yvals[0:(self.nx*self.ny),:]; 
        vvals = self.yvals[(self.nx*self.ny):(2*self.nx*self.ny),:];
        uvals = np.transpose(uvals); vvals = np.transpose(vvals);
        YMESH, TVALS = np.meshgrid(self.ymesh,self.tvals);
        TVALS = TVALS.reshape((self.tvals.size*self.ymesh.size));
        YMESH = YMESH.reshape((self.tvals.size*self.ymesh.size));
        ykeys, tkeys = [ii for ii in YMESH], [jj for jj in TVALS];
        tupleKeys = list(zip(ykeys,tkeys));
        strKeys = ["%f, %f" % x for x in tupleKeys]
        uvals = uvals.reshape((self.ny*self.tvals.size,self.nx)); 
        vvals = vvals.reshape((self.ny*self.tvals.size,self.nx));
        uvalsDF = pd.DataFrame(uvals,columns=self.xmesh,index=strKeys)
        vvalsDF = pd.DataFrame(vvals,columns=self.xmesh,index=strKeys)
        # NEED TO FIGURE OUT HOW TO GIVE UVALS AND VVALS THE CORRECT KEYS
        tvalsDF = pd.DataFrame(self.tvals); xmeshDF = pd.DataFrame(self.xmesh); 
        ymeshDF = pd.DataFrame(self.ymesh); uvalsDF.to_excel(writer,'uvals'); 
        vvalsDF.to_excel(writer,'vvals'); tvalsDF.to_excel(writer,'tvals'); 
        xmeshDF.to_excel(writer,'xmesh'); ymeshDF.to_excel(writer,'ymesh');
        paramsDF = pd.DataFrame({'TE':self.TE,'TI':self.TI,
            'SE':self.SE,'SI':self.SI,'TAU':self.TAU,'mode':self.mode,
            'kernType':self.kernType,'span_x':self.span_x,'span_y':self.span_y,
            'AEE':self.AEE,'AIE':self.AIE,'AEI':self.AEI,'AII':self.AII,
            'dx_kern':self.dx_kern,'dy_kern':self.dy_kern,'dt':self.dt,
            'Lx':self.Lx,'Ly':self.Ly,'nx':self.nx,'ny':self.ny},index=['Parameters']);
        paramsDF.to_excel(writer,'params'); 
        
        " Clearing all integration history and reinitialize to final values in last integration "
    def refreshFromLast(self,dt=None,tmore=tmore_default,tshow=tshow_default):
        if dt == None: dt = self.dt
        " Creating the WilsonCowan2D object "
        pardict = {'TE':self.TE,'TI':self.TI,'SE':self.SE, 
                   'AEE': self.AEE,'AIE':self.AIE,'AEI':self.AEI,'AII':self.AII, 
                   'Lx':self.Lx,'Ly':self.Ly,'nx':self.nx,'ny':self.ny, 
                   'span_x':self.span_x,'span_y':self.span_y, 
                   'dx_kern':self.dx_kern,'dy_kern':self.dy_kern,
                   'SI':self.SI,'TAU':self.TAU,'dt':dt, 
                   'kernType':self.kernType, 'mode':self.mode};
        WC2D = WilsonCowan2D(pardict=pardict); 
        uvals = self.yvals[0:(self.nx*self.ny),:]; 
        vvals = self.yvals[(self.nx*self.ny):(2*self.nx*self.ny),:]
        uvals = np.reshape(uvals,(self.ny,self.nx,self.tvals.size)); 
        vvals = np.reshape(vvals,(self.ny,self.nx,self.tvals.size));
        WC2D.setInitConds(0,uvals[:,:,-1],vvals[:,:,-1],tmore=tmore,tshow=tshow);
        print("Final Conditions from Saved Data made into Initial Conditions of New Simulation");
        return WC2D
    
    " Interactive real-time plotting by I/O dialog  "
    def interactiveIO_img_3D(self): 
        if self.WCSystem == None: 
            video = bool(int(input('Video? (1:YES; 0:NO): ')));
            self.WilsonCowanIntegrator();
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs];
            self.imgWilsonCowan3D(tdisp, ydisp, video);
        self.tmore = float(input("How much more integration time: ")); 
        self.tshow = float(input("How much time to display?: "));
        video = bool(int(input('Video? (1:YES; 0:NO): ')));
        cont= bool(int(input('Continue? (1:YES; 0:NO): ')));
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        while(cont):
            self.WCIntegrateLast(); locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs];
            self.imgWilsonCowan3D(tdisp, ydisp, video);
            self.tmore = float(input("How much more integration time: ")); 
            self.tshow = float(input("How much time to display?: "));
            video = bool(int(input('Video? (1:YES; 0:NO): ')));
            cont= bool(int(input('Continue? (1:YES; 0:NO): ')));
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename);
    def interactiveIO_img_2D_mesh(self,pltType=0):
        if self.WCSystem == None:
            self.WilsonCowanIntegrator();
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs[0]]; ydisp = self.yvals[:,locs[0]];
            self.imgWilsonCowan2D(tdisp, ydisp, pltType); 
        self.tmore = float(input("How much more integration time: ")); 
        self.tshow = float(input("How much time to display?: "));
        pltType = int(input('Plot Type? (0:Normal, 1:Polar, 2:Torus): '));
        cont= bool(int(input('Continue? (1:YES; 0:NO): ')));
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        while(cont):
            self.WCIntegrateLast();
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs[0]]; ydisp = self.yvals[:,locs[0]];
            self.imgWilsonCowan2D(tdisp, ydisp, pltType); 
            self.tmore = float(input("How much more integration time: ")); 
            self.tshow = float(input("How much time to display?: "));
            pltType = int(input('Plot Type? (0:Normal, 1:Polar, 2:Torus): '));
            cont= bool(int(input('Continue? (1:YES; 0:NO): ')));
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename);
    def interactiveIO_img_vs_T(self,axis=0,indx=0,pltType=0): # HANDLE THIS FUNCTION WITH CARE
        if self.WCSystem == None:             
            self.WilsonCowanIntegrator();
            yvals_mesh = np.reshape(self.yvals,(2*self.ny,self.nx,np.size(self.tvals)));
            uvals_mesh = yvals_mesh[0:self.ny,:,:]; 
            vvals_mesh = yvals_mesh[self.ny:(2*self.ny),:,:];
            if axis == 0:
                uvals_vs_t = uvals_mesh[indx,:,:]; 
                vvals_vs_t = vvals_mesh[indx,:,:]; 
            elif axis == 1: 
                uvals_vs_t = uvals_mesh[:,indx,:]; 
                vvals_vs_t = vvals_mesh[:,indx,:]; 
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            udisp = uvals_vs_t[:,locs]; vdisp = vvals_vs_t[:,locs];
            tdisp = self.tvals[locs]; ydisp = np.concatenate((udisp,vdisp))
            self.imgWilsonCowan1DvsT(tdisp,ydisp,(axis, indx, pltType)); 
        self.tmore = float(input("How much more integration time: "));
        self.tshow = float(input("How much time to display?: "));
        axis = int(input('Y-Slice(0) or X-Slice(1)?: '));
        if axis == 1:
            indx = int(input("Which X-Slice (0 to "+str(self.nx-1)+"): ")); 
        elif axis == 0:
            indx = int(input("Which Y-Slice (0 to "+str(self.ny-1)+"): "));
        pltType = int(input("Plot Type (0-Normal, 1-Polar, 2-Cylinder): "));
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        while(self.tshow>0):
            self.WCIntegrateLast();
            yvals_mesh = np.reshape(self.yvals,(2*self.ny,self.nx,np.size(self.tvals)));
            uvals_mesh = yvals_mesh[0:self.ny,:,:]; 
            vvals_mesh = yvals_mesh[self.ny:(2*self.ny),:,:];
            if axis == 0:
                uvals_vs_t = uvals_mesh[indx,:,:]; 
                vvals_vs_t = vvals_mesh[indx,:,:];
            elif axis == 1: 
                uvals_vs_t = uvals_mesh[:,indx,:]; 
                vvals_vs_t = vvals_mesh[:,indx,:]; 
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            udisp = uvals_vs_t[:,locs]; vdisp = vvals_vs_t[:,locs];
            tdisp = self.tvals[locs]; ydisp = np.concatenate((udisp,vdisp))
            self.imgWilsonCowan1DvsT(tdisp, ydisp, (axis, indx, pltType)); 
            self.tmore = float(input("How much more integration time: "));
            self.tshow = float(input("How much time to display?: "));
            axis = int(input('Y-Slice(0) or X-Slice(1)?: '));
            if axis == 1:
                indx = int(input("Which X-Slice (0 to "+str(self.nx-1)+"): ")); 
            elif axis == 0:
                indx = int(input("Which Y-Slice (0 to "+str(self.ny-1)+"): ")); 
            pltType = int(input("Plot Type (0-Normal, 1-Polar, 2-Cylinder): "));
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename);
    def interactiveIO_pltT(self,locIndx=(0,0)):
        if self.WCSystem == None: 
            x_indx, y_indx = locIndx; self.WilsonCowanIntegrator();
            yvals_mesh = np.reshape(self.yvals,(2*self.ny,self.nx,np.size(self.tvals)));
            uvals_mesh = yvals_mesh[0:self.ny,:,:]; 
            vvals_mesh = yvals_mesh[self.ny:(2*self.ny),:,:];
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            uvals_disp = uvals_mesh[y_indx,x_indx,locs]; tdisp = self.tvals[locs];
            vvals_disp = vvals_mesh[y_indx,x_indx,locs];
            self.pltWilsonCowan1DvsT(tdisp, uvals_disp, vvals_disp); 
        self.tmore = float(input("How much more integration time: "));
        self.tshow = float(input("How much time to display?: "));
        x_indx = int(input("Which X-Slice (0 to "+str(self.nx-1)+"): ")); 
        y_indx = int(input("Which Y-Slice (0 to "+str(self.ny-1)+"): "));
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        while(self.tshow>0):
            self.WCIntegrateLast();
            yvals_mesh = np.reshape(self.yvals,(2*self.ny,self.nx,np.size(self.tvals)));
            uvals_mesh = yvals_mesh[0:self.ny,:,:]; 
            vvals_mesh = yvals_mesh[self.ny:(2*self.ny),:,:];
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            uvals_disp = uvals_mesh[y_indx,x_indx,locs]; tdisp = self.tvals[locs];
            vvals_disp = vvals_mesh[y_indx,x_indx,locs];
            self.pltWilsonCowan1DvsT(tdisp, uvals_disp, vvals_disp); 
            self.tmore = float(input("How much more integration time: "));
            self.tshow = float(input("How much time to display?: "));
            x_indx = int(input("Which X-Slice (0 to "+str(self.nx-1)+"): ")); 
            y_indx = int(input("Which Y-Slice (0 to "+str(self.ny-1)+"): ")); 
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename);   
    def interactiveIO_anim_2D_mesh(self, polar=False, tshift=4):
        tpause = 1e-9; tspeed = (tshift,tpause);
        if self.WCSystem == None:
            self.WilsonCowanIntegrator();
            locs = np.arange(self.tvals.size);
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs];
            self.animWilsonCowan2D(tdisp, ydisp, polar, tspeed); 
        self.tmore = float(input("How much more integration time: ")); 
        self.tshow = float(input("How much time to display?: "));
        tshift = int(input("Time Shift Step for Animation (Integer >= 1): "));
        polar = bool(int(input('Polar Plot? (1:YES, 0:NO): ')));
        cont= bool(int(input('Continue? (1:YES; 0:NO): ')));
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        while(cont):
            self.WCIntegrateLast();
            locs = np.arange(self.tvals.size); tspeed = (tshift,tpause);
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            tdisp = self.tvals[locs]; ydisp = self.yvals[:,locs];
            self.animWilsonCowan2D(tdisp, ydisp, polar, tspeed); 
            self.tmore = float(input("How much more integration time: ")); 
            self.tshow = float(input("How much time to display?: "));
            tshift = int(input("Time Shift Step for Animation (Integer >= 1): "));
            polar = bool(int(input('Polar Plot? (1:YES, 0:NO): ')));
            cont= bool(int(input('Continue? (1:YES; 0:NO): ')));
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename);
    def interactiveIO_anim_vs_T(self,axis=0,indx=0,polar=False,tshift=4): 
        if self.WCSystem == None:             
            self.WilsonCowanIntegrator(); tspeed = (tshift,1e-9);
            yvals_mesh = np.reshape(self.yvals,(2*self.ny,self.nx,np.size(self.tvals)));
            uvals_mesh = yvals_mesh[0:self.ny,:,:]; 
            vvals_mesh = yvals_mesh[self.ny:(2*self.ny),:,:];
            if axis == 0:
                uvals_vs_t = uvals_mesh[indx,:,:]; 
                vvals_vs_t = vvals_mesh[indx,:,:]; 
            elif axis == 1: 
                uvals_vs_t = uvals_mesh[:,indx,:]; 
                vvals_vs_t = vvals_mesh[:,indx,:]; 
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            udisp = uvals_vs_t[:,locs]; vdisp = vvals_vs_t[:,locs];
            tdisp = self.tvals[locs]; ydisp = np.concatenate((udisp,vdisp))
            self.animWilsonCowan1DvsT(tdisp,self.twin,ydisp,axis,indx,polar,tspeed=tspeed); 
        self.tmore = float(input("How much more integration time: "));
        self.tshow = float(input("How much total time in animation?: "));
        self.twin = float(input("How much display time?: "));
        tshift = int(input("Time Shift Step for Animation (Integer >= 1): "));
        polar = bool(int(input("Plot Type (0-Normal, 1-Polar): ")));
        axis = int(input('Y-Slice(0) or X-Slice(1)?: ')); tspeed = (tshift,1e-9);
        if axis == 1:
            indx = int(input("Which X-Slice (0 to "+str(self.nx-1)+"): ")); 
        elif axis == 0:
            indx = int(input("Which Y-Slice (0 to "+str(self.ny-1)+"): "));
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        while(self.tshow>0):
            self.WCIntegrateLast();
            yvals_mesh = np.reshape(self.yvals,(2*self.ny,self.nx,np.size(self.tvals)));
            uvals_mesh = yvals_mesh[0:self.ny,:,:]; 
            vvals_mesh = yvals_mesh[self.ny:(2*self.ny),:,:];
            if axis == 0:
                uvals_vs_t = uvals_mesh[indx,:,:]; 
                vvals_vs_t = vvals_mesh[indx,:,:];
            elif axis == 1: 
                uvals_vs_t = uvals_mesh[:,indx,:]; 
                vvals_vs_t = vvals_mesh[:,indx,:]; 
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            udisp = uvals_vs_t[:,locs]; vdisp = vvals_vs_t[:,locs];
            tdisp = self.tvals[locs]; ydisp = np.concatenate((udisp,vdisp))
            self.animWilsonCowan1DvsT(tdisp,self.twin,ydisp,axis,indx,polar,tspeed=tspeed); 
            self.tmore = float(input("How much more integration time: "));
            self.tshow = float(input("How much total time in animation?: "));
            self.twin = float(input("How much display time?: "));
            tshift = int(input("Time Shift Step for Animation (Integer >= 1): "));
            polar = bool(int(input("Plot Type (0-Normal, 1-Polar): ")));
            axis = int(input('Y-Slice(0) or X-Slice(1)?: ')); tspeed = (tshift,1e-9);
            if axis == 1:
                indx = int(input("Which X-Slice (0 to "+str(self.nx-1)+"): ")); 
            elif axis == 0:
                indx = int(input("Which Y-Slice (0 to "+str(self.ny-1)+"): "));
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename);
    # animWilsonCowanPltT(self,tdisp,twin,udisp,vdisp,tspeed=(4,1e-9))
    def interactiveIO_anim_pltT(self,locIndx=(0,0),tshift=4):
        if self.WCSystem == None: 
            x_indx, y_indx = locIndx; self.WilsonCowanIntegrator();
            yvals_mesh = np.reshape(self.yvals,(2*self.ny,self.nx,np.size(self.tvals)));
            uvals_mesh = yvals_mesh[0:self.ny,:,:]; 
            vvals_mesh = yvals_mesh[self.ny:(2*self.ny),:,:];
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            uvals_disp = uvals_mesh[y_indx,x_indx,locs]; tdisp = self.tvals[locs];
            vvals_disp = vvals_mesh[y_indx,x_indx,locs]; tspeed = (tshift,1e-9);
            self.animWilsonCowanPltT(tdisp,self.twin,uvals_disp,vvals_disp,tspeed=tspeed); 
        self.tmore = float(input("How much more integration time: "));
        self.tshow = float(input("How much total time in animation?: "));
        self.twin = float(input("How much display time?: "));
        tshift = int(input("Time Shift Step for Animation (Integer >= 1): "));
        x_indx = int(input("Which X-Slice (0 to "+str(self.nx-1)+"): ")); 
        y_indx = int(input("Which Y-Slice (0 to "+str(self.ny-1)+"): "));
        shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
        if shouldSave: 
            filename = str(input("Excel file to save data to: "));
            self.saveCurrentState(filename);
        while(self.tshow>0):
            self.WCIntegrateLast();
            yvals_mesh = np.reshape(self.yvals,(2*self.ny,self.nx,np.size(self.tvals)));
            uvals_mesh = yvals_mesh[0:self.ny,:,:]; 
            vvals_mesh = yvals_mesh[self.ny:(2*self.ny),:,:];
            locs = np.arange(self.tvals.size); 
            locs = locs[self.tvals>(self.tvals[-1]-self.tshow)];
            uvals_disp = uvals_mesh[y_indx,x_indx,locs]; tdisp = self.tvals[locs];
            vvals_disp = vvals_mesh[y_indx,x_indx,locs]; tspeed = (tshift,1e-9);
            self.animWilsonCowanPltT(tdisp,self.twin,uvals_disp,vvals_disp,tspeed=tspeed);  
            self.tmore = float(input("How much more integration time: "));
            self.tshow = float(input("How much total time in animation?: "));
            self.twin = float(input("How much display time?: "));
            tshift = int(input("Time Shift Step for Animation (Integer >= 1): "));
            x_indx = int(input("Which X-Slice (0 to "+str(self.nx-1)+"): ")); 
            y_indx = int(input("Which Y-Slice (0 to "+str(self.ny-1)+"): ")); 
            shouldSave = bool(int(input("Save Data? (1:YES or 0:NO): ")));
            if shouldSave: 
                filename = str(input("Excel file to save data to: "));
                self.saveCurrentState(filename); 
    
    " Defining Heaviside or Unit Step Function "
    def heaviside(self,x): return 0.5 * (np.sign(x) + 1)
        
    " Defining sigmoidal activation function "
    def f(self,x): return self.heaviside(x) #1/(1+np.exp(-self.BETA*x))
    
    " Defining Convolution Kernel "
    def kern(self,sigma,span,diff): 
        span_x, span_y = span; dx, dy = diff;
        x = dx*np.arange(-span_x,span_x+1);
        y = dy*np.arange(-span_y,span_y+1);
        X, Y = np.meshgrid(x,y);
        if self.kernType == 'Exponential':
            return (1/(2*np.pi*sigma**2))*np.exp(-np.sqrt(X**2+Y**2)/sigma)*dx*dy;
        if self.kernType == 'Gaussian':
            return (1/(sigma**2))*np.exp(-np.pi*(X**2+Y**2)/(sigma**2))*dx*dy;
        
    " Various Fourier Transforms of the above kernel "
    def freqz2(self,hh,ww_cols,ww_rows):
        kk_rows = np.arange(np.shape(hh)[0]); kk_cols = np.arange(np.shape(hh)[1]);
        KK_rows = np.matrix(kk_rows); KK_cols = np.transpose(np.matrix(kk_cols));
        WW_rows = np.transpose(np.matrix(ww_rows)); WW_cols = np.matrix(ww_cols)
        frqz_indx_rows = WW_rows*KK_rows; frqz_indx_cols = KK_cols*WW_cols;
        frqz_mtx_rows = np.matrix(np.exp(-1j*np.array(frqz_indx_rows)));
        frqz_mtx_cols = np.matrix(np.exp(-1j*np.array(frqz_indx_cols)));
        HH = frqz_mtx_rows*np.matrix(hh)*frqz_mtx_cols; 
        W_cols, W_rows = np.meshgrid(ww_cols, ww_rows); 
        return np.array(np.abs(HH)), W_cols, W_rows;
    def kernFT(self,sig,LL,nn,numFTpts,span,cycles):
        nx, ny = nn; Lx, Ly = LL; dx = Lx/(nx-1); dy = Ly/(ny-1); fs_x = 1/dx; 
        fs_y = 1/dy; numFTpts_x, numFTpts_y = numFTpts; span_x, span_y = span;
        numFTcycles_x, numFTcycles_y = cycles; # For Fourier Transform of Spatial Mesh
        w_x = np.linspace(-numFTcycles_x*fs_x*np.pi,numFTcycles_x*fs_x*np.pi,numFTpts_x); 
        w_y = np.linspace(-numFTcycles_y*fs_y*np.pi,numFTcycles_y*fs_y*np.pi,numFTpts_y); 
        wnorm_x = w_x/fs_x; wnorm_y = w_y/fs_y; 
        diff = (dx,dy); ker = self.kern(sig,span,diff);  
        kerFT, ww_x, ww_y = self.freqz2(ker,wnorm_x,wnorm_y); ww = (ww_x*fs_x,ww_y*fs_y);
        return ww, np.abs(kerFT), ker, diff
    def kernAnalyticFT(self,sig,ww): 
        ww_x, ww_y = ww; 
        if self.kernType == 'Exponential':
            return 1/((1+(sig**2)*(ww_x**2+ww_y**2))**(3/2));
        if self.kernType == 'Gaussian':
            return np.exp(-(sig**2)*(ww_x**2 + ww_y**2)/(4*np.pi));  
    def compareFTs(self,SI,LL,nn,numFTpts,span,cycle):
        ww, kerFT, ker, diff = self.kernFT(SI,LL,nn,numFTpts,span,cycle);
        kerFTanalytical = self.kernAnalyticFT(SI,ww); 
        fig = plt.figure(); ax = fig.gca(projection='3d');
        ax.plot_wireframe(ww[0],ww[1],kerFT,rstride=3,cstride=3,
                          color='r',label = 'Discrete-Space Fourier Transform');
        ax.plot_wireframe(ww[0],ww[1],kerFTanalytical,rstride=2,cstride=2,
                          color='b',label = 'Continuous-Space Fourier Transform');
        ax.set_xlabel('Spatial Frequency [Inverse Length] /nin X Direction');
        ax.set_ylabel('Spatial Frequency [Inverse Length] /nin Y Direction');
        ax.set_zlabel('Fourier Amplitude')
        plt.legend(framealpha=0); plt.show(); return ww, kerFT, ker, diff;    
    # Example Code for this:
    # ww, kerFT, ker, diff = compareFTs(0.05,(4,4),(90,90),(90,90),(90,90),(2,2))
    # ww, kerFT, ker, diff = compareFTs(0.1,(4,4),(100,100),(200,200),(50,50),(3,3))

""" THE REST OF THIS FILE IS JUST A SCRIPT TO TEST THE CLASS-LEVEL
IMPLEMENTATIONS OF MY ALGORITHMS. ANYTHING BELOW THIS POINT ONLY 
SERVES AS EXAMPLE CODE AND NOTHING ELSE. EVERYTHING ABOVE THIS POINT 
IS SELF-SUFFICIENT CODE WHiCH MAY BE IMPORTED INTO ANOTHER PYTHON FILE"""
        
" FIX CLKWISE/CTRCLKWISE PROBLEM IN 1D POLAR PLOTS "

" First: ANIMATE THESE PLOTS "       
        
" Second: Look at Harmonic Excitation based on Fourier Transforms in 1D "

" Third: Look at Lattice shapes in 2D pattern formation  "

" Fourth: Make a GUI for this "


" RECYCLING CODE "
#dftmtx = lambda N: np.fft.fft(np.eye(N))
#TO COMPARE MY CUSTOM BUILT FREQZ TO SCIPY'S:
#def myFrqz(b,wpts):
#    ww = np.linspace(-np.pi,np.pi,wpts); kk = np.arange(np.size(b)); 
#    w, h = signal.freqz(b,worN=ww)
#    WW = np.matrix(ww); KK = np.matrix(kk); KK = np.transpose(KK);
#    frqz_indx = KK*WW; frqz_mtx = np.matrix(np.exp(-1j*np.array(frqz_indx)));
#    HH = np.matrix(b)*frqz_mtx; return np.array(HH)[0,:], h, ww;    
    

