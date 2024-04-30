import time, sys
import raytrace_cuda as rt
import numpy as np
import pylab as pl
from pylab import *

t = rt.implane([1,10],
               #[-2.0, -0.01, 0.0, 0.01],
               [-2.0, -0.001, 0.0, 0.001],
               (215,0,0),
               25.0,
               140e6,
               mode='Tbr')

before = time.time()

t.trace(2000, tpts='All')

print
print 'elapsed time: %8.3f s' % (time.time()-before)

#np.save('solar_burst_tbr.npy', t.tbr)
#np.save('solar_burst_rays.npy', t.traj)


yleft = -1.8;  yright = -1.6
xfar  = -0.1;   xnear =  0.1

bur = []
rnum = []

for i in xrange(t.nrays):
    ix = np.where((t.traj[i,:,1] > yleft) & (t.traj[i,:,1] < yright) & \
                  (t.traj[i,:,0] > xfar) & (t.traj[i,:,1] < xnear))[0]
    if len(ix) > 0:
        rnum.append(i)
        bur.append((i, ix))

rnum = np.array(rnum, dtype=int)

figure(); grid(1); axis('equal')

al = linspace(-pi, pi, 100)
plot(cos(al), sin(al), 'orange', lw=2)
fill(array([xfar, xfar, xnear, xnear]),
     array([yright, yleft, yleft, yright]), 'r', alpha=0.3)

for i in rnum:
    plot(t.traj[i,:,0], t.traj[i,:,1])
for i in xrange(t.nrays):
    plot(t.traj[i,:,0], t.traj[i,:,1])



show()
