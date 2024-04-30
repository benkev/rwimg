#
# Exploring the Stokes V parameter behavior 
#
from pylab import *
import raytrace as rt; import sys; import os


a = rt.implane((40,1),[0., 0.5, 1., 0.5], rsph = 25.0, freq=200e6, cnu=3.,
               mode='TbrIQUV', msun=1.,
               trkrays=[[0,0],[1,0],[2,0],[3,0],[4,0],[5,0],
                        [6,0],[7,0],[8,0],[9,0],[10,0],[11,0],
                        [12,0],[13,0],[14,0],[15,0],[16,0],[17,0],
                        [18,0],[19,0]],
               trkparms=['Stokes', 'Arclen', 'dist', 'tbriquv', 'bmagn',
                         'bfield', 'angbd'], trknpmax=3000)

a.trace(3000)


abdm = nanmin(a.traj.angbd)
col = ['b', 'g', 'r', 'c', 'm', 'y', 'orange', 'brown', 'violet', 'pink', 'k']
nc = len(col)
iz = empty(20, dtype=int)
zx = empty(20, dtype=double)
zy = empty(20, dtype=double)
hpi = pi/2.

figure();
for i in xrange(20):
    plot(a.traj.arclen[i,:], a.traj.tbriquv[i,:,3], col[i%nc])
    plot(a.traj.arclen[i,:], 2e3*(a.traj.angbd[i,:]-abdm) - 15000, col[i%nc])
    x = copy(a.traj.angbd[i,:])
    ix = find(~isnan(x))
    x = x[ix]
    ix = find(x > hpi); x[:] = 0.; x[ix] = 1.
    iiz = find(diff(x) != 0.)
    if len(iiz) > 0:
        iz[i] = iiz[0]     # Zero crossing point index
        zx[i] = a.traj.arclen[i,iz[i]]  # Zero crossing point
        zy[i] = a.traj.tbriquv[i,iz[i],3]  # Zero crossing point of Tb V
        plot([zx[i],zx[i]], [-15000, 7000], 'k-')
        plot(zx[i], zy[i], col[i%nc], marker='o')
grid(1)
l90 = 2e3*(pi/2. - abdm) - 15000
plot([0.,nanmax(a.traj.arclen[0,:])], [l90, l90], 'k')
title('Stokes V')

figure();
for i in xrange(20):
    plot(a.traj.arclen[i,:], a.traj.tbriquv[i,:,2], col[i%nc])
    plot(a.traj.arclen[i,:], 2e3*(a.traj.angbd[i,:]-abdm) - 15000, col[i%nc])
    x = copy(a.traj.angbd[i,:])
    ix = find(~isnan(x))
    x = x[ix]
    ix = find(x > hpi); x[:] = 0.; x[ix] = 1.
    iiz = find(diff(x) != 0.)
    if len(iiz) > 0:
        iz[i] = iiz[0]     # Zero crossing point index
        zx[i] = a.traj.arclen[i,iz[i]]  # Zero crossing point
        zy[i] = a.traj.tbriquv[i,iz[i],2]  # Zero crossing point of Tb U
        plot([zx[i],zx[i]], [-15000, 7000], 'k-')
        plot(zx[i], zy[i], col[i%nc], marker='o')
grid(1)
l90 = 2e3*(pi/2. - abdm) - 15000
plot([0.,nanmax(a.traj.arclen[0,:])], [l90, l90], 'k')
title('Stokes U')

show()
