#!/usr/bin/env python

"""
Seismic: 2D finite difference simulation of membrane wave propagation.

Difraction example in cylindrical wedge model. Based on:
R. M. Alford, K. R. Kelly and D. M. Boore -
Accuracy of finite-difference modeling of the acoustic wave equation.
Geophysics  1974
"""
import numpy as np
# from matplotlib import animation
import wavefd
# from fatiando.vis import mpl
import wavesolver2d as wave2d
## Check: spacing, V model, padding
shape = (481, 481)
# shape = (1001, 1001)
ds = 1000.  # spacing
duration = 300.0
padding=200
dt=0.1
fc=0.1; 
area = [0, (shape[0]-1) * ds, 0, (shape[1]-1) * ds]
# Velocity Model
myVecolity=wave2d.VelocityModel(shape, ds, V=3500)

# Xc=200;
# Zc=200;
# R=50;
# Va=2000;
# myVecolity.CircleGradualAnomaly( Xc=Xc, Zc=Zc, R=R, Va=Va)

# Source
mysources = wave2d.Sources(shape=shape, ds=ds)
# mysources.SingleSource(Nx=100, Nz=100, fc=fc)
mysources.SingleSource(amp=100., Nx=100, Nz=240, fc=fc)
# mysources.SingleSource(Nx=105, Nz=75, fc=fc)

#Receivers
myRec=wave2d.Receivers(shape=shape, ds=ds)
# myRec.SingleReceiver(Nx=75, Nz=125)
myRec.HomoStaLst(5,5)
# 
# # INData=wave2d.InputDataBase(shape=shape, ds=ds, dt=dt, duration=duration, padding=padding, sources=mysources, fmin=fc,
# #         receivers=myRec, velocityM=myVecolity, snapshots=10)
INData=wave2d.InputDataBase(shape=shape, ds=ds, dt=dt, duration=duration, padding=padding, sources=mysources, fmin=fc,
        receivers=myRec, velocityM=myVecolity)
membraneWS=wave2d.WaveSolverMembrane(INData)
stations=membraneWS.Run_sh(outmp4='/lustre/janus_scratch/life9360/001.mp4')
membraneWS.SaveSeismograms(outdir='/lustre/janus_scratch/life9360/MEM2D_homo/SAC_003')
