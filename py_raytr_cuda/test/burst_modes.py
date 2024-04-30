#
# Type II solar bursts modes
#
import raytrace_cuda as rt
from pylab import *


t = rt.implane([5,5],
               #[-2.0,-2.0,2.0,2.0],  # rect=[xleft, ybottom, xright, ytop]
               [-1.2,-0.5,-0.8,-0.4],  # rect=[xleft, ybottom, xright, ytop]
               (215,0,0),
               25.0,
               100e6,
               mode='Tbr')  #,
               #traj=True)

t.trace(3000)

#tp = array([[2,3],[3,2]],dtype=int32)
#trj = t.trace(3000, tracePoints=tp)

#trj = trace.trace(3000, tracePoints=array([[23,23],[26,26]],dtype=int))
#trace.trace(3000)

## f = open('out.out', 'w')
## for i in range(0,len(trace.tbr)):
##     f.write(str(trace.tbr[i]))
## f.close()

