#!/usr/bin/python

import raytrace_cuda as rt


trace = rt.implane([500,500],
                   [-2.0,-2.0,2.0,2.0],
                   (215,0,0),
                   25.0,
                   200e6,
                   mode='Tbr')

trace.trace(3000)
#fig = figure()
#imshow(trace.tbr[:,:],cmap=cm.gist_heat)
#colorbar()
#show()
