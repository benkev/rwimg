import raytrace as rt
from pylab import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

p = np.zeros((2,3))
p = np.zeros((2,2))
p[0,0]
p[0,0] =25
p[0,1] =15
p[1,0] =25
p[1,1] =35

a = rt.implane(grid=(50,50),rect=(-2,-2,2,2),freq=300e6,mode='TbrIV')
traj = a.trace(1500,p)

prof = a.plprofile(traj[0,:,:])


#plt.quiver(traj[0,:,0],traj[0,:,2],prof[2][:,0],prof[2][:,1])
#plt.show()


########################

n = 200
        
xr = linspace(-5, 5, n)
yr = linspace(-5, 5, n)

X,Y = meshgrid(xr,yr)

M = np.zeros(X.shape, dtype='bool')
print(M[n/2,n/2])
for i in range(0, n):
    for j in range(0, n):
        if(np.sqrt(xr[i]**2+yr[j]**2)<=1):
            #print (i,j,xr[i],yr[j])
            M[i,j] = True

print(M[n/2,n/2])

xpass = np.reshape(X,n**2)
ypass = np.reshape(Y,n**2)

array = np.empty((n**2,3))
array[:,0] = xpass
array[:,2] = ypass
array[:,1] = np.zeros(n**2)

#print(array)

prof = a.plprofile(array)

B = prof[2]
rho = prof[0]


B2 = np.reshape(B,(n,n,3))

plt.figure(2)


for i in range(0, n):
    for j in range(0, n):
        if(np.sqrt(xr[i]**2+yr[j]**2)<=1):
            #print (i,j,xr[i],yr[j])
            B2[i,j,0] = 0
            B2[i,j,2] = 0

plt.quiver(X,Y, B2[:,:,0], B2[:,:,2], units='inches', #'dots',
           pivot='middle', scale=10, minlength=.01)
axis('equal')
circle = linspace(-2*np.pi,2*np.pi,100)
plt.plot(np.cos(circle),np.sin(circle))

plt.plot(traj[0,:,0],traj[0,:,2])

plt.plot(traj[1,:,0],traj[1,:,2])

plt.show()
