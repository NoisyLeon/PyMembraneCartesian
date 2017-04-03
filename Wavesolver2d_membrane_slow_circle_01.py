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
myVecolity.CircleCosineAnomaly( Xc=Xc, Zc=Zc, R=a, va=2500)
velocity=myVecolity.velocity


fc = 0.1
sources = [wavefd.GaussSource(10 * ds, 400 * ds, area, shape,  1., fc)]
dt = wavefd.scalar_maxdt(area, shape, np.max(velocity))
duration = 1500.0
maxit = int(duration / dt)
snapshots = 5  # every 3 iterations plots one
# stations = [[400 * ds, 400 * ds], [285 * ds, 675 * ds], [700 * ds, 400 * ds], [697 * ds, 460 * ds], [658 * ds, 635 * ds]]
# stations = [[400 * ds, 400 * ds], [285 * ds, 675 * ds], [700 * ds, 400 * ds], [699 * ds, 424 * ds],[658 * ds, 635 * ds]]

stations=[[400 * ds, 400 * ds], [425 * ds, 400 * ds], [450 * ds, 400 * ds], [475 * ds, 400 * ds], [500 * ds, 400 * ds], [700 * ds, 400 * ds]]
simulation = wavefd.membrane(
    velocity, area, dt, maxit, sources, stations=None, snapshot=snapshots, padding=50)

# This part makes an animation using matplotlibs animation API
background = (velocity - 4000) * 5* 10 ** 2
fig = mpl.figure(figsize=(16, 12))
mpl.subplots_adjust(right=0.98, left=0.11, hspace=0.5, top=0.93)
wavefield = mpl.imshow(np.zeros_like(velocity), extent=area,
                       cmap=mpl.cm.seismic_r, vmin=-1500000, vmax=1500000)
mpl.points(stations, '^b', size=8)
mpl.ylim(area[2:][::-1])
mpl.xlabel('x (km)')
mpl.ylabel('z (km)')
mpl.m2km()



def animate(i):
    t, u, seismogram = simulation.next()
    # seismogram1.set_data(times[:t + 1], seismogram[0][:t + 1])
    wavefield.set_array(background[::-1] + u[::-1])
    return wavefield


anim = animation.FuncAnimation(
    fig, animate, frames=maxit / snapshots, interval=200)
anim.save('wavefield_circle_R_600_01.mp4', dpi=400)
mpl.show()

