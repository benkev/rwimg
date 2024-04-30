import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
#ax = fig.gca(projection='3d')
ax = fig.add_subplot(111, projection='3d')

theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
z = np.linspace(-2, 2, 100)
r = z**2 + 1
x = r * np.sin(theta)
y = r * np.cos(theta)
ax.plot(x, y, z, label='parametric curve')
ax.legend()

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()

print 'Axis limits:'
print 'ax.get_xlim3d() = ', ax.get_xlim3d()
print 'ax.get_ylim3d() = ', ax.get_ylim3d()
print 'ax.get_ylim3d() = ', ax.get_ylim3d()

xls = ax.get_xlim3d(); xl = xls[1] - xls[0]
yls = ax.get_ylim3d(); yl = yls[1] - yls[0]
zls = ax.get_zlim3d(); zl = zls[1] - zls[0]
rr = np.sqrt(xl**2 + yl**2 + zl**2)
r1 = (xl/rr, yl/rr, zl/rr)
print r1, np.dot(r1,r1)

