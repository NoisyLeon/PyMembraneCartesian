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

# Set the parameters of the finite difference grid
shape = (500, 500)
ds = 100.  # spacing
area = [0, shape[0] * ds, 0, shape[1] * ds]
# Set the parameters of the finite difference grid
# velocity = np.zeros(shape) + 3000.
# velocity[75:125, 75:125] = 2500.


shape = (500, 500)
ds = 100.  # spacing
myVecolity=ws2d.VelocityModel(shape, ds, V=4000)
Xc=100
Zc=100
dv=-0.1
a=20*ds
b=10*ds
myVecolity.CosineEllipse( Xc=25000, Zc=20000, a=2000, b=8000, va=2500)
velocity=myVecolity.velocity


fc = 8.
sources = [wavefd.GaussSource(10 * ds, 250 * ds, area, shape,  1., fc)]
dt = wavefd.scalar_maxdt(area, shape, np.max(velocity))
duration = 15.0
maxit = int(duration / dt)
# stations = [[10 * ds, 100 * ds]]  # x, z coordinate of the seismometer
stations = [[200 * ds, 250 * ds], [146 * ds, 385 * ds], [400 * ds, 250 * ds], [398 * ds, 283 * ds], [347 * ds, 445 * ds]]  # x, z coordinate of the seismometer
Nr = len(stations)
snapshots = 5  # every 3 iterations plots one
simulation = wavefd.membrane(
    velocity, area, dt, maxit, sources, stations, snapshots, padding=50)

# This part makes an animation using matplotlibs animation API
background = (velocity - 4000) * 10 ** -1
fig = mpl.figure(figsize=(16, 12))
mpl.subplot(Nr+1, 1, 1)
mpl.subplots_adjust(right=0.98, left=0.11, hspace=0.5, top=0.93)
# mpl.subplot2grid((4, 3), (0, 0), colspan=3, rowspan=3)
wavefield = mpl.imshow(np.zeros_like(velocity), extent=area,
                       cmap=mpl.cm.seismic_r, vmin=-500, vmax=500)
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
    seismogram1, = mpl.plot([], [], '-k')
    seisLst.append(seismogram1)
    mpl.xlim(0, duration)
    mpl.ylim(-400, 400)
    mpl.ylabel('Amplitude')
times = np.linspace(0, dt * maxit, maxit)
# This function updates the plot every few timesteps


def animate(i):
    t, u, seismogram = simulation.next()
    for i in xrange(Nr):
        seisLst[i].set_data(times[:t + 1], seismogram[i][:t + 1])
    wavefield.set_array(background[::-1] + u[::-1])
    return wavefield


anim = animation.FuncAnimation(
    fig, animate, frames=maxit / snapshots, interval=1)
mpl.show()

