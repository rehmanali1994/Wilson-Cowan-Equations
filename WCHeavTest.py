# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 18:36:00 2015

@author: rehman
"""

from WilsonCowanSystemsHeaviside import *
from WilsonCowanSystemsFFTHeaviside import *
import time
# I'm importing numpy, scipy, matplotlib, pandas, visvis, etc. in addition to 
# WilsonCowan1D (from WilsonCowan1D import WilsonCowan1D). The Wilson Cowan 1D 
# above is really just the python filename (WilsonCowan1D.py) but not the class

def WC1Dmain():
    " CONSTANTS: DO NOT CHANGE! "
    TE = 0.08; TI = 0.4; SE = 10; #TE=0.125
    AEE = 1;  AIE = 1; AEI = 1.5; AII = 0.25;
    
    " Defining ODE Parameters "
    L = 400; # Length of mesh
    nx = 400; # Spaces in mesh
    dx = L/(nx-1); # Spacing of mesh
    xx = np.linspace(0,L,nx); # the mesh itself
    span = 60; # Point to Left and Right of Central Maximum in Kernel
    dx_kern = dx; # Spacing of Points in the Kernel
    SI = 9; TAU = 0.6; # Remaining ODE Parameters SI=8; TAU=0.6
    dt = 0.001; tmore = 60; # Time Step and Initlal Integration Time
    tshow = 60; # Amount of Time Data to Display
    kernType='Exponential' # 'Gaussian' or 'Exponential'
    mode = 'wrap'; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
    
    " Creating the WilsonCowan1D object "
    pardict = {'TE':TE, 'TI':TI, 'SE':SE, 'AEE': AEE, 'AIE':AIE,
               'AEI':AEI, 'AII':AII, 'L':L, 'nx':nx, 'span':span, 'dx_kern':dx_kern,
               'SI':SI, 'TAU': TAU, 'dt':dt, 'kernType':kernType, 'mode':mode};
    WC1D = WilsonCowan1D(pardict=pardict); xx = WC1D.xmesh
              
    " Initial Conditions for ODE (COMMENT AND UNCOMMENT THE FOLLOWING 4 ICs AS NEEDED) "
#    # For Synchronous Solution:
#    u0 = 0.44*np.ones((nx)); v0 = 0.23*np.ones((nx));
    
#    # For Standing Wave (AKA: Honeycomb)
#    u0 = 0.0001*np.random.randn((nx))+0.44; 
#    v0 = 0.0001*np.random.randn((nx))+0.23;
    
#    # For Travelling Wave in One Direction:
#    u0 = 0.001*np.random.randn((nx))+0.2*np.cos(2*np.pi*xx/L)+0.44; 
#    v0 = 0.001*np.random.randn((nx))+0.1*np.sin(2*np.pi*xx/L)+0.23;
    
    # For Travelling Wave in Opposite Direction:
    u0 = 0.001*np.random.randn((nx))+0.2*np.sin(2*np.pi*xx/L)+0.44; 
    v0 = 0.001*np.random.randn((nx))+0.1*np.cos(2*np.pi*xx/L)+0.23;

    " DON'T UNCOMMENT THIS! "
    t0 = 0; WC1D.setInitConds(t0,u0,v0,tmore=tmore,tshow=tshow);
    
    " These perform the actual integration "
    WC1D.interactiveIO_img();
#    WC1D.interactiveIO_anim();
#    WC1D.interactiveIO_animX();
#    WC1D.interactiveIO_animT();

    return WC1D
    
def WC1DmainFFT():
#    " CONSTANTS: DO NOT CHANGE! "
#    TE = 0.25; TI = 0.55; SE = 45; #TE=0.125
#    AEE = 1;  AIE = 1; AEI = 1.5; AII = 0.25;
#    
#    " Defining ODE Parameters "
#    L = 400; # Length of mesh
#    span = 200; # Point to Left and Right of Central Maximum in Kernel
#    nx = 2*span + 1; # Number of finite elements in mesh
#    dx = L/(nx-1); # Spacing of mesh
#    xx = np.linspace(0,L,nx); # the mesh itself
#    dx_kern = dx; # Spacing of Points in the Kernel
#    SI = 30; TAU = 0.80; # Remaining ODE Parameters SI=8; TAU=0.6
#    dt = 0.025; tmore = 10; # Time Step and Initlal Integration Time
#    tshow = 10; # Amount of Time Data to Display
#    kernType='Gaussian' # 'Gaussian' or 'Exponential'
#    mode = 'wrap'; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
    
    " CONSTANTS: DO NOT CHANGE! "
    TE = 0.05; TI = 0.4; SE = 12; #TE=0.125
    AEE = 1;  AIE = 1; AEI = 1.5; AII = 0.25;
    
    " Defining ODE Parameters "
    L = 400; # Length of mesh
    span = 200; # Point to Left and Right of Central Maximum in Kernel
    nx = 2*span + 1; # Number of finite elements in mesh
    dx = L/(nx-1); # Spacing of mesh
    xx = np.linspace(0,L,nx); # the mesh itself
    dx_kern = dx; # Spacing of Points in the Kernel
    SI = 8; TAU = 0.5; # Remaining ODE Parameters SI=8; TAU=0.6
    dt = 0.025; tmore = 10; # Time Step and Initlal Integration Time
    tshow = 10; # Amount of Time Data to Display
    kernType='Exponential' # 'Gaussian' or 'Exponential'
    mode = 'wrap'; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
    
    " Creating the WilsonCowan1D object "
    pardict = {'TE':TE, 'TI':TI, 'SE':SE, 'AEE': AEE, 'AIE':AIE,
               'AEI':AEI, 'AII':AII, 'L':L, 'span':span, 'dx_kern':dx_kern,
               'SI':SI, 'TAU': TAU, 'dt':dt, 'kernType':kernType, 'mode':mode};
    WC1D = WilsonCowan1D_FFT(pardict=pardict); xx = WC1D.xmesh
              
    " Initial Conditions for ODE (COMMENT AND UNCOMMENT THE FOLLOWING 4 ICs AS NEEDED) "
#    # For Synchronous Solution:
#    u0 = 0.44*np.ones((nx)); v0 = 0.23*np.ones((nx));
    
#    # For Standing Wave (AKA: Honeycomb)
#    u0 = 0.0001*np.random.randn((nx))+0.44; 
#    v0 = 0.0001*np.random.randn((nx))+0.23;
    
#    # For Travelling Wave in One Direction:
#    u0 = 0.001*np.random.randn((nx))+0.2*np.cos(2*np.pi*xx/L)+0.44; 
#    v0 = 0.001*np.random.randn((nx))+0.1*np.sin(2*np.pi*xx/L)+0.23;
    
    # For Travelling Wave in Opposite Direction:
    u0 = 0.001*np.random.randn((nx))+0.2*np.sin(2*np.pi*xx/L)+0.44; 
    v0 = 0.001*np.random.randn((nx))+0.1*np.cos(2*np.pi*xx/L)+0.23;

    " DON'T UNCOMMENT THIS! "
    t0 = 0; WC1D.setInitConds(t0,u0,v0,tmore=tmore,tshow=tshow);
    
    " These perform the actual integration "
    WC1D.interactiveIO_img();
#    WC1D.interactiveIO_anim();
#    WC1D.interactiveIO_animX();
#    WC1D.interactiveIO_animT();
    
    return WC1D;
    
def WC2Dmain():
    " Defining ODE Parameters "
    L = 50; Lx = L; Ly = L; # Dimensions of mesh
    n = 100; nx = n; ny = n; # Spaces in mesh
    #dx = Lx/(nx-1); dy = Ly/(ny-1); # Spacing of mesh
    #xx = np.linspace(0,L,nx); # the x-component of the mesh 
    #yy = np.linspace(0,L,ny); # the y-component of the mesh
    #XX, YY = np.meshgrid(xx, yy); # the complete 2D mesh
    span = 15; # Point to Left and Right of Central Maximum in 1D Kernel
    span_x = span; span_y = span; # Same thing but for the 2D Kernel
    dx_kern = 1; dy_kern = 1; # Spacing of Points in the Kernel
    SI = 4; TAU = 0.5; # Remaining ODE Parameters
    dt = 0.1; tmore = 10; # Time Step and Initlal Integration Time
    tshow = 10; # Amount of Time Data to Display
    kernType='Gaussian' # 'Gaussian' or 'Exponential'
    mode = 'wrap'; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
    
    " CONSTANTS: DO NOT CHANGE! "
    TE = 0.125; TI = 0.4; SE = 6;
    AEE = 1;  AIE = 1; AEI = 1.5; AII = 0.25;
    
    " Creating the WilsonCowan2D object "
    pardict = {'TE':TE, 'TI':TI, 'SE':SE, 'AEE': AEE, 'AIE':AIE,
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
    #print('2D Plotting Coming Up ...')
    #WC2D.interactiveIO_img_2D_mesh()
    #print('3D Plotting/Animation Coming Up ...')
    WC2D.interactiveIO_img_3D()
    #WC2D.interactiveIO_anim_2D_mesh();
    #WC2D.interactiveIO_anim_pltT();
    #WC2D.interactiveIO_anim_vs_T();
    return WC2D
    
def WC2DmainFFT():
    " Defining ODE Parameters "
    L = 50; Lx = L; Ly = L; # Dimensions of mesh
    #dx = Lx/(nx-1); dy = Ly/(ny-1); # Spacing of mesh
    #xx = np.linspace(0,L,nx); # the x-component of the mesh 
    #yy = np.linspace(0,L,ny); # the y-component of the mesh
    #XX, YY = np.meshgrid(xx, yy); # the complete 2D mesh
    span = 170; # Point to Left and Right of Central Maximum in 1D Kernel
    span_x = span; span_y = span; # Same thing but for the 2D Kernel
    n = 2*span+1; nx = n; ny = n; # Spaces in mesh
    dx_kern = L/span; dy_kern = L/span; # Spacing of Points in the Kernel
    SI = 2; TAU = 0.5; # Remaining ODE Parameters
    dt = 0.1; tmore = 10; # Time Step and Initlal Integration Time
    tshow = 10; # Amount of Time Data to Display
    kernType='Exponential' # 'Gaussian' or 'Exponential'
    mode = 'reflect'; #'wrap' for periodic boundary, 'reflect' for reflecting boundary
    
    " CONSTANTS: DO NOT CHANGE! "
    TE = 0.125; TI = 0.4; SE = 3;
    AEE = 1;  AIE = 1; AEI = 1.5; AII = 0.25;
    
    " Creating the WilsonCowan2D object "
    pardict = {'TE':TE, 'TI':TI, 'SE':SE, 'AEE': AEE, 'AIE':AIE,
               'AEI':AEI, 'AII':AII, 'Lx':Lx, 'Ly':Ly, 'span_x':span_x, 
               'span_y':span_y, 'dx_kern':dx_kern, 'dy_kern':dy_kern,
               'SI':SI, 'TAU': TAU, 'dt':dt, 'kernType':kernType, 'mode':mode};
    WC2D = WilsonCowan2D_FFT(pardict=pardict); 
    XX = WC2D.XMesh; YY = WC2D.YMesh; WC2D.aspectRatio = 20;
    
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
    #print('2D Plotting Coming Up ...')
    #WC2D.interactiveIO_img_2D_mesh()
    #print('3D Plotting/Animation Coming Up ...')
    WC2D.interactiveIO_img_3D()
    WC2DFromLast = WilsonCowan2D_FFT(pardict=pardict);
    uvals = WC2D.yvals[0:(WC2D.nx*WC2D.ny),:]; 
    vvals = WC2D.yvals[(WC2D.nx*WC2D.ny):(2*WC2D.nx*WC2D.ny),:]
    uvals = np.reshape(uvals,(WC2D.ny,WC2D.nx,WC2D.tvals.size)); 
    vvals = np.reshape(vvals,(WC2D.ny,WC2D.nx,WC2D.tvals.size));
    WC2DFromLast.setInitConds(t0,uvals[:,:,-1],vvals[:,:,-1],tmore=tmore,tshow=tshow);
    WC2DFromLast.interactiveIO_img_3D();
    #WC2D.interactiveIO_anim_2D_mesh();
    #WC2D.interactiveIO_anim_pltT();
    #WC2D.interactiveIO_anim_vs_T();
    return WC2D
    
    
#def WC2DFFTresume(filename,dt,tmore,tshow):
#    WClast = WilsonCowan2D_FFT(filename=filename);
#    print("Saved Data Loaded ...")
#    " Creating the WilsonCowan2D object "
#    pardict = {'BETA':WClast.BETA,'TE':WClast.TE,'TI':WClast.TI,'SE':WClast.SE, 
#               'AEE': WClast.AEE, 'AIE':WClast.AIE,'AEI':WClast.AEI,'AII':WClast.AII, 
#               'Lx':WClast.Lx,'Ly':WClast.Ly,'span_x':WClast.span_x,'span_y':WClast.span_y, 
#               'dx_kern':WClast.dx_kern, 'dy_kern':WClast.dy_kern,
#               'SI':WClast.SI,'TAU':WClast.TAU,'dt':dt,
#               'kernType':WClast.kernType, 'mode':WClast.mode};
#    WC2D = WilsonCowan2D_FFT(pardict=pardict); 
#    XX = WC2D.XMesh; YY = WC2D.YMesh; WC2D.aspectRatio = 10;
#    uvals = WClast.yvals[0:(WClast.nx*WClast.ny),:]; 
#    vvals = WClast.yvals[(WClast.nx*WClast.ny):(2*WClast.nx*WClast.ny),:]
#    uvals = np.reshape(uvals,(WClast.ny,WClast.nx,WClast.tvals.size)); 
#    vvals = np.reshape(vvals,(WClast.ny,WClast.nx,WClast.tvals.size));
#    WC2D.setInitConds(0,uvals[:,:,-1],vvals[:,:,-1],tmore=tmore,tshow=tshow);
#    print("Final Conditions from Saved Data made into Initial Conditions of New Simulation");
#    WC2D.interactiveIO_img_3D();
#    return WC2D
    

#def fft1Dconv(x,kern,mode):
#    if mode == 'wrap':
#        raw_conv = np.fft.ifftn(np.fft.fftn(x)*np.fft.fftn(kern));
#        return np.real(np.roll(np.fft.fftshift(raw_conv),1));
#    if mode == 'reflect':
#        x_withreflects = np.concatenate((x[::-1],x,x[::-1]));
#        kern_zeropadding = np.zeros(kern.size);
#        kern_zeropads = np.concatenate((kern_zeropadding,kern,kern_zeropadding));
#        raw_conv = self.fft1Dconv(x_withreflects,kern_zeropads,'wrap');
#        return raw_conv[kern.size:(2*kern.size)];
#def fft2Dconv(x,kern,mode):
#    if mode == 'wrap':
#        raw_conv = np.fft.ifftn(np.fft.fftn(x)*np.fft.fftn(kern));
#        conv_shifted = np.real(np.fft.fftshift(raw_conv));
#        return np.roll(np.roll(conv_shifted,1,axis=0),1,axis=1);
#    if mode == 'reflect':
#        x_flipedLR = np.fliplr(x);
#        x_concatLR = np.concatenate((x_flipedLR,x,x_flipedLR),axis=1);
#        x_flipedUD = np.flipud(x_concatLR);
#        x_withreflects = np.concatenate((x_flipedUD,x_concatLR,x_flipedUD));
#        zeropadsLR = np.zeros(kern.shape);
#        kern_zeropadsLR = np.concatenate((zeropadsLR,kern,zeropadsLR),axis=1);
#        zeropadsUD = np.zeros(kern_zeropadsLR.shape);
#        kern_zeropads = np.concatenate((zeropadsUD,kern_zeropadsLR,zeropadsUD));
#        raw_conv = fft2Dconv(x_withreflects,kern_zeropads,'wrap');
#        return raw_conv[kern.shape[0]:(2*kern.shape[0]),kern.shape[1]:(2*kern.shape[1])];
#c = np.random.randn(1001); d = np.random.randn(1001);
#init_t = time.time(); h2 = fft1Dconv(c,d,mode='wrap'); time.time()-init_t
#init_t = time.time(); h1 = conv(c,d,mode='wrap'); time.time()-init_t
#init_t = time.time(); h2 = fft1Dconv(c,d,mode='reflect'); time.time()-init_t
#init_t = time.time(); h1 = conv(c,d,mode='reflect'); time.time()-init_t
#a = np.random.randn(101,101); b = np.random.randn(101,101);
#init_t = time.time(); h2 = fft2Dconv(a,b,mode='wrap'); time.time()-init_t
#init_t = time.time(); h1 = conv(a,b,mode='wrap'); time.time()-init_t
#init_t = time.time(); h2 = fft2Dconv(a,b,mode='reflect'); time.time()-init_t
#init_t = time.time(); h1 = conv(a,b,mode='reflect'); time.time()-init_t


#" How to open simulation from file "
#WC1D2 = WilsonCowan1D(filename='test.xlsx');
#print("JUST REOPENED FILE INTO NEW OBJECT")
#WC1D2.interactiveIO_img();


#" WAVE SPEED CODE "
#uvals = WC1D.yvals[0:WC1D.nx,:];
#vvals = WC1D.yvals[WC1D.nx:(2*WC1D.nx),:];
#tindex_start = 2280; time_steps = 30
#x1 = WC1D.dx*np.argmax(uvals[:,tindex_start]);
#x2 = WC1D.dx*np.argmax(uvals[:,(tindex_start+time_steps)]);
#speed = (x2-x1)/(time_steps*WC1D.dt);


WC1D = WC1Dmain();
#WC2D = WC2Dmain();
#WC1D = WC1DmainFFT();
#WC2D = WC2DmainFFT();
#WC2D = WC2DFFTresume('NoPatFormFFT1.xlsx',0.025,5,5)


