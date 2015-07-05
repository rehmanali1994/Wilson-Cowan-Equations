# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 18:36:00 2015

@author: rehman
"""

from WilsonCowanSystems import * 
# I'm importing numpy, scipy, matplotlib, pandas, visvis, etc. in addition to 
# WilsonCowan1D (from WilsonCowan1D import WilsonCowan1D). The Wilson Cowan 1D 
# above is really just the python filename (WilsonCowan1D.py) but not the class

def WC1Dmain():
    " CONSTANTS: DO NOT CHANGE! "
    BETA = 50; TE = 0.05; TI = 0.4; SE = 12; #TE=0.125
    AEE = 1;  AIE = 1; AEI = 1.5; AII = 0.25;
    
    " Defining ODE Parameters "
    L = 400; # Length of mesh
    nx = 400; # Spaces in mesh
    dx = L/(nx-1); # Spacing of mesh
    xx = np.linspace(0,L,nx); # the mesh itself
    span = 60; # Point to Left and Right of Central Maximum in Kernel
    dx_kern = 1; # Spacing of Points in the Kernel
    SI = 10.8; TAU = 0.7; # Remaining ODE Parameters SI=8; TAU=0.6
    dt = 0.025; tmore = 10; # Time Step and Initlal Integration Time
    tshow = 10; # Amount of Time Data to Display
    kernType='Exponential' # 'Gaussian' or 'Exponential'
    mode = 'wrap'; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
    
    " Creating the WilsonCowan1D object "
    pardict = {'BETA':BETA, 'TE':TE, 'TI':TI, 'SE':SE, 'AEE': AEE, 'AIE':AIE,
               'AEI':AEI, 'AII':AII, 'L':L, 'nx':nx, 'span':span, 'dx_kern':dx_kern,
               'SI':SI, 'TAU': TAU, 'dt':dt, 'kernType':kernType, 'mode':mode};
    WC1D = WilsonCowan1D(pardict=pardict); xx = WC1D.xmesh
              
    " Initial Conditions for ODE "
    # For Synchronous Solution:
#    u0 = 0.4*np.ones((nx)); v0 = 0.2*np.ones((nx));
    # For Travelling Wave in One Direction:
    # SI = 4 and SE = 6 with travelling wave in one direction
    # produces both a honeycomb and travelling wave
#    u0 = 0.01*np.random.randn((nx))+0.1*np.sin(2*np.pi*xx/L)+0.4; 
#    v0 = 0.01*np.random.randn((nx))+0.1*np.cos(2*np.pi*xx/L)+0.2;
    # For Travelling Wave in Opposite Direction:
    u0 = 0.01*np.random.randn((nx))+0.1*np.cos(2*np.pi*xx/L)+0.4; 
    v0 = 0.01*np.random.randn((nx))+0.1*np.sin(2*np.pi*xx/L)+0.2;
    # For Standing Wave (AKA: Honeycomb)
#    u0 = 0.1*np.random.randn((nx))+0.41; 
#    v0 = 0.1*np.random.randn((nx))+0.21;
    t0 = 0; WC1D.setInitConds(t0,u0,v0);
    WC1D.interactiveIO_img();
#    WC1D2 = WilsonCowan1D(filename='test.xlsx');
#    print("JUST REOPENED FILE INTO NEW OBJECT")
#    WC1D2.interactiveIO_img();
    
def WC2Dmain():
    " Defining ODE Parameters "
    L = 100; Lx = L; Ly = L; # Dimensions of mesh
    n = 50; nx = n; ny = n; # Spaces in mesh
    #dx = Lx/(nx-1); dy = Ly/(ny-1); # Spacing of mesh
    #xx = np.linspace(0,L,nx); # the x-component of the mesh 
    #yy = np.linspace(0,L,ny); # the y-component of the mesh
    #XX, YY = np.meshgrid(xx, yy); # the complete 2D mesh
    span = 15; # Point to Left and Right of Central Maximum in 1D Kernel
    span_x = span; span_y = span; # Same thing but for the 2D Kernel
    dx_kern = 1; dy_kern = 1; # Spacing of Points in the Kernel
    SI = 2; TAU = 0.6; # Remaining ODE Parameters
    dt = 0.1; tmore = 10; # Time Step and Initlal Integration Time
    tshow = 1; # Amount of Time Data to Display
    kernType='Gaussian' # 'Gaussian' or 'Exponential'
    mode = 'wrap'; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
    
    " CONSTANTS: DO NOT CHANGE! "
    BETA = 50; TE = 0.125; TI = 0.4; SE = 3;
    AEE = 1;  AIE = 1; AEI = 1.5; AII = 0.25;
    
    " Creating the WilsonCowan2D object "
    pardict = {'BETA':BETA, 'TE':TE, 'TI':TI, 'SE':SE, 'AEE': AEE, 'AIE':AIE,
               'AEI':AEI, 'AII':AII, 'Lx':Lx, 'Ly':Ly, 'nx':nx, 'ny':ny, 
               'span_x':span_x, 'span_y':span_y, 'dx_kern':dx_kern, 'dy_kern':dy_kern,
               'SI':SI, 'TAU': TAU, 'dt':dt, 'kernType':kernType, 'mode':mode};
    WC2D = WilsonCowan2D(pardict=pardict); XX = WC2D.XMesh; YY = WC2D.YMesh; 
    
    " Initial Conditions for ODE "
    # For Synchronous Solution:
    #u0 = 0.4*np.ones((nx,ny)); v0 = 0.2*np.ones((nx,ny));
    # For Standing Wave (AKA: Bathroom Tiles in Random Directions):
    u0 = 0.05*np.random.randn(nx,ny)+0.41; 
    v0 = 0.05*np.random.randn(nx,ny)+0.21;
    # For Directional Standing Wave Patterns (Theta Degrees Above X-Axis)
    #theta=np.pi/6; zFreq = 4; ZZ = XX*np.cos(theta)-YY*np.sin(theta);
    #u0 = 0.1*np.sin(2*np.pi*zFreq*ZZ/L)+0.4;
    #v0 = 0.1*np.sin(2*np.pi*zFreq*ZZ/L)+0.2; #(can change a sin into a cos)
    # For Line impulse along diagonol (requires nx = ny) THIS LEAVE A DIAGONOL SINUSOID:
    #u0 = 0.4+0.1*np.identity(nx); v0 = 0.2+0.1*np.identity(nx);
    # For X-Shaped Line Impulses along diagonals (requires nx = ny) ANGLE ZERO PATTERN:
    #u0 = 0.4+0.05*np.identity(nx)+0.05*np.flipud(np.identity(nx)); 
    #v0 = 0.2+0.05*np.identity(nx)+0.05*np.flipud(np.identity(nx));
    # For Travelling Wave in 1D:
    #u0 = 0.05*np.sin(2*np.pi*XX/Lx)+0.4;
    #v0 = 0.05*np.cos(2*np.pi*XX/Lx)+0.2;
    # Ripples in a pond:
    #u0 = 0.4*np.ones((nx,ny)); v0 = 0.2*np.ones((nx,ny));
    #u0[np.round(ny/2),np.round(nx/2)] = 0.6;
    #v0[np.round(ny/2),np.round(nx/2)] = 0.1;
    t0 = 0; WC2D.setInitConds(t0,u0,v0,tmore=tmore,tshow=tshow);
    WC2D.interactiveIO_img_2D_mesh()
    WC2D.interactiveIO_img_3D()

WC1Dmain();
#WC2Dmain();