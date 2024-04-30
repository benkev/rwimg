import raytrace_cuda as rt
import numpy as np

for a in [0,90]:
    for th in np.linspace(0,90,7):
        for pi in np.linspace(0,90,7):
            if(a==0):
                rt.run_freq_range(file='data/orientationX/theta'+str(th)+'_phi'+str(pi),theta=th,phi=pi,orientation=a,grid=(60,60), steps=12)
            else:
                rt.run_freq_range(file='data/orientationY/theta'+str(th)+'_phi'+str(pi),theta=th,phi=pi,orientation=a,grid=(60,60), steps=12)
