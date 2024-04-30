#
# A quick and dirty work around to emulate axis('equal') in mplot3d
#

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure()
ax = Axes3D(fig)

Scale = 50
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

x = 10 * np.outer(np.cos(u), np.sin(v))
y = 10 * np.outer(np.sin(u), np.sin(v))
z = 13 * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b')



Xstart=(x.min())
Xend  =(x.max())

Ystart=(y.min())
Zstart=(z.min())

Scalex=(Xend-Xstart)


## ax.set_xlim3d(Xstart,Xstart+Scale)
## ax.set_ylim3d(Ystart,Ystart+Scale)
## ax.set_zlim3d(Zstart,Zstart+Scale)
ax.set_xlim3d(-10,10)
ax.set_ylim3d(-10,10)
ax.set_zlim3d(-10,10)

plt.show()
