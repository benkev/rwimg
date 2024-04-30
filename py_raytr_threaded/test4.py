import raytrace as rt
import numpy as np
import math
from pylab import *

trace = rt.implane([60,60],
                   [-2.0,-2.0,2.0,2.0],
                   (215,0,0),
                   25.0,
                   80e6,
                   mode='TbrIV')
trace.set_plfunc('plasma_parameters.c')
trace.remove_streamers()

trace.make_streamer(90,
                    90,
                    0)
                    

trace.trace(1500)
fig = figure()
imshow(trace.tbriv[:,:,0],cmap=cm.gist_heat)
colorbar()
show()
