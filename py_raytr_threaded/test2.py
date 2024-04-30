import raytrace as rt
from pylab import *
import numpy as np
import math


rt.run_freq_range(file='data/scattering_blah',grid=(5,5), mode='Tbr',steps=12,scattering=False)

'''
for a in [0,90]:
    for th in np.linspace(15,90,6):
        for pi in np.linspace(0,90,7):
            if(a==0):
                rt.run_freq_range(file='data/orientationX/theta'+str(th)+'_phi'+str(pi),theta=th,phi=pi,orientation=a,grid=(300,300), mode='Tbr',steps=12,scattering=True)
            else:
                rt.run_freq_range(file='data/orientationY/theta'+str(th)+'_phi'+str(pi),theta=th,phi=pi,orientation=a,grid=(300,300), mode='Tbr',steps=12,scattering=True)
'''            
'''
a = rt.implane(grid=(40,40),rect=(-2,-2,2,2),mode='TbrIV')
a.make_streamer(45,0)
#a.make_streamer(90,90)

a.trace(1500)
figure()
imshow(a.tbriv[:,:,1])
show()

'''
