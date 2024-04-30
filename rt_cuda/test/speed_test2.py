#!/usr/bin/python

import raytrace_cuda as rt

import sys
import time


for size in [10, 20, 50, 100, 200]:
    trace = rt.implane(grid=[size,size],
                   rect=[-2.0,-2.0,2.0,2.0],
                   obs=(215,0,0),
                   rsph=25.0,
                   freq=80e6,
                   mode='TbrIQUV')

    before = time.time()


    trace.trace(1500)
    
    sys.stderr.write(str(size*size)+": "+str(time.time()-before)+"\n")
