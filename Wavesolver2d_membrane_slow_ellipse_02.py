"""
Seismic: 2D finite difference simulation of membrane wave propagation.

Difraction example in cylindrical wedge model. Based on:
R. M. Alford, K. R. Kelly and D. M. Boore -
Accuracy of finite-difference modeling of the acoustic wave equation.
Geophysics  1974
"""
import numpy as np
from matplotlib import animation
import wavefd
from fatiando.vis import mpl
import wavesolver2d as ws2d
import matplotlib.pyplot as plt
from scipy.signal import hilbert, chirp

# Set the parameters of the finite difference grid
shape = (800, 800)
ds = 6000.  # spacing
area = [0, shape[0] * ds, 0, shape[1] * ds]
myVecolity=ws2d.VelocityModel(shape, ds, V=4000)
Xc=400*ds
Zc=400*ds
dv=-0.1
a=100*ds
b=20*ds
# myVecolity.CosineEllipse( Xc=25000, Zc=20000, a=2000, b=8000, va=3000)
# myVecolity.CosineEllipse( Xc=2500000, Zc=2000000, a=200000, b=800000, va=3000)
# myVecolity.CosineEllipse( Xc=Xc, Zc=Zc, a=b, b=a, va=3000)
myVecolity.CircleCosineAnomaly( Xc=Xc, Zc=Zc, R=a, va=3500)
velocity=myVecolity.velocity


fc = 0.1
sources = [wavefd.GaussSource(10 * ds, 400 * ds, area, shape,  1., fc)]
dt = wavefd.scalar_maxdt(area, shape, np.max(velocity))
duration = 1500.0
maxit = int(duration / dt)
# stations = [[400 * ds, 400 * ds], [285 * ds, 675 * ds], [700 * ds, 400 * ds], [697 * ds, 460 * ds], [658 * ds, 635 * ds]]  # x, z coordinate of the seismometer
# stations = [[400 * ds, 400 * ds], [285 * ds, 675 * ds], [700 * ds, 400 * ds], [699 * ds, 424 * ds],[658 * ds, 635 * ds]]  # x, z coordinate of the seismometer
# stations = [[400 * ds, 400 * ds], [450 * ds, 400 * ds], [500 * ds, 400 * ds], [550 * ds, 400 * ds], [600 * ds, 400 * ds], [650 * ds, 400 * ds]]
stations=[[400 * ds, 400 * ds], [425 * ds, 400 * ds], [450 * ds, 400 * ds], [475 * ds, 400 * ds], [500 * ds, 400 * ds], [700 * ds, 400 * ds]]
# xArr = np.arange(300, 601, 25)
# for x in xArr:
#     stations.append([x * ds, 400 * ds])

Nr = len(stations)
snapshots = 5  # every 3 iterations plots one
simulation = wavefd.membrane(
    velocity, area, dt, maxit, sources, stations, snapshots, padding=50)

# This part makes an animation using matplotlibs animation API
background = (velocity - 4000) * 5*10 ** 2
fig = mpl.figure(figsize=(16, 12))
mpl.subplot(Nr+1, 1, 1)
mpl.subplots_adjust(right=0.98, left=0.11, hspace=0.5, top=0.93)
# mpl.subplot2grid((4, 3), (0, 0), colspan=3, rowspan=3)
wavefield = mpl.imshow(np.zeros_like(velocity), extent=area,
                       cmap=mpl.cm.seismic_r, vmin=-1500000, vmax=1500000)
mpl.points(stations, '^b', size=8)
# mpl.points([[125 * ds, 75 * ds]], '*g', size=8)
mpl.ylim(area[2:][::-1])
mpl.xlabel('x (km)')
mpl.ylabel('z (km)')
mpl.m2km()
# mpl.subplot2grid((4, 3), (3, 0), colspan=3)

Nr = len(stations)
seisLst = []
for i in xrange(Nr):
    mpl.subplot(Nr+1, 1, i+2)
    # mpl.subplot(Nr, 1, i+1)
    # mpl.subplot2grid((4, 3), (0, 0), colspan=3, rowspan=1)
    seismogram1, = mpl.plot([], [], '-k')
    seisLst.append(seismogram1)
    mpl.xlim(0, duration)
    mpl.ylim(-2000000, 2000000)
    mpl.ylabel('Amplitude')
times = np.linspace(0, dt * maxit, maxit)
# This function updates the plot every few timesteps


def animate(i):
    global seismogram
    t, u, seismogram = simulation.next()
    for i in xrange(Nr):
        seisLst[i].set_data(times[:t + 1], seismogram[i][:t + 1])
    wavefield.set_array(background[::-1] + u[::-1])
    return wavefield


anim = animation.FuncAnimation(
    fig, animate, frames=maxit / snapshots, interval=200)
anim.save('ellipse_wavefield_circle_002.mp4', dpi=400)
mpl.show()

fig = plt.figure(figsize=(16, 12))
plt.subplot(3, 1, 1)
plt.plot(seismogram[0], 'b-', lw=5)
plt.subplot(3, 1, 2)
plt.plot(-hilbert(seismogram[0]), 'r-', lw=5)
plt.subplot(3, 1, 3)
plt.plot(seismogram[-1], 'k-', lw=5)


