from pylab import *
#from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys

# from http://www.xsi-blog.com/archives/115

def fibonacci(N):
    inc = np.pi * (3 - np.sqrt(5))
    off = 2. / N
    r2d = 180./np.pi
    k = np.arange(0,N)
    y = k*off - 1. + 0.5*off
    r = np.sqrt(1 - y*y)
    phi = k * inc
    x = np.cos(phi)*r
    z = np.sin(phi)*r
    theta = np.arctan2(np.sqrt(x**2+y**2),z)
    phi = np.arctan2(y,x)
    lats = 90.-r2d*theta
    lons = r2d*phi
    return lats, lons

#npts = int(sys.argv[1])
#print sys.argv[1]
npts = 1000

lats, lons = fibonacci(npts)
lam = radians(lats)
phi = radians(lons)
r = 1;
x = r*cos(lam)*cos(phi)
y = r*cos(lam)*sin(phi)
z = r*sin(lam)

fig = plt.figure()
ax = Axes3D(fig)

ax.scatter3D(x, y, z, zdir='z', c = 'r')

plt.show()
