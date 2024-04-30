import raytrace_cuda as rt
from pylab import *

trace = rt.implane([50,50],
                   [-2.0,-2.0,2.0,2.0],
                   (215,0,0),
                   25.0,
                   80e6,
                   mode='Tbr',
                   traj=True)

traj = trace.trace(3000, \
        tracePoints=array([[20,25],[21,25],[22,25],[23,25],[24,25]], dtype=int))
f = open('out.out', 'w')
for i in range(0,len(trace.tbr)):
    f.write(str(trace.tbr[i]))
