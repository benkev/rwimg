import raytrace as rt
import numpy as np
import math


#rt.run_freq_range(file='data/goodheaders/scattering',grid=(300,300), mode='Tbr'#,steps=12,scattering=True)

th = 60
pi = 90
rt.run_freq_range(file='data/goodheaders/theta'+str(th)+'_phi'+str(pi),theta=th,phi=pi,orientation = 0,grid=(300,300), mode='Tbr',steps=12,scattering=True)


th = 60
pi = 60
rt.run_freq_range(file='data/goodheaders/theta'+str(th)+'_phi'+str(pi),theta=th,phi=pi,orientation = 0,grid=(300,300), mode='Tbr',steps=12,scattering=True)



th = 45
pi = 30
rt.run_freq_range(file='data/goodheaders/theta'+str(th)+'_phi'+str(pi),theta=th,phi=pi,orientation = 0,grid=(300,300), mode='Tbr',steps=12,scattering=True)

