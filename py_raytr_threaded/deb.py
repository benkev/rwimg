import pylab as pl
import numpy as np
import raytrace as rt
import sys

a = rt.implane((20,20),[-2.,-2.,2.,2.], rsph = 25.0, freq=275e6, cnu=3.,
               mode='TbrIQUV', trkrays=[[2,5],[3,3],[3,5],[5,5]],
               trkparms=['Stokes', 'Arclen', 'dist'], trknpmax=3000)

a.trace(3000)

