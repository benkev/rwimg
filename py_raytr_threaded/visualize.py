
import raytrace as rt

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



from pylab import *
import numpy as np
from matplotlib.colors import LogNorm

size = 1000

a = rt.implane((50,50),(-2,-2,2,2),mode='TbrIV')
rt.make_streamer(90,90,orientation=0,density=5,baseStrength=10.0,stalkStrength=2.0,scaleX=.75,scaleY=.75,scaleZ=.75)
#a.trace(3000)
#imshow(a.tbriv[:,:,0])

## b = np.empty((size,3))
## b[:,0] = linspace(1.1,2,size)
## c = a.plprofile(b)
## plt.figure(1)
## plt.plot(b[:,0],c[0])
## plt.show()

####################################
# BIG X-Y

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
array[:,1] = ypass
array[:,2] = np.zeros(n**2)

#print(array)

prof = a.plprofile(array)

B = prof[2]
rho = prof[0]


B2 = np.reshape(B,(n,n,3))
rho2 = np.reshape(rho,(n,n))



plt.figure(2)


mag = np.empty((n,n))

mag[:,:] = np.sqrt(B2[:,:,0]**2+B2[:,:,2]**2+B2[:,:,1]**2)

mag = np.ma.masked_array(mag, mask=M)

plt.contour(X,Y, mag, 30,colors='k')

circle = linspace(-2*np.pi,2*np.pi,100)

plt.plot(np.cos(circle),np.sin(circle))

plt.show()

plt.figure(2)

B2[:,:,0] = np.ma.masked_array(B2[:,:,0], mask=M)
B2[:,:,1] = np.ma.masked_array(B2[:,:,1], mask=M)

plt.quiver(X,Y, B2[:,:,0], B2[:,:,1], units='inches', #'dots',
           pivot='middle', scale=10, minlength=.01)
axis('equal')
plt.plot(np.cos(circle),np.sin(circle))

plt.show()

plt.figure(3)

rho2 = np.ma.masked_array(rho2, mask=M)

plt.contour(X,Y, rho2, 200,colors='k')
axis('equal')
plt.plot(np.cos(circle),np.sin(circle))

plt.show()


####################################
# X-Y

n = 76
        
xr = linspace(-1.5, 1.5, n)
yr = linspace(.5,10, n)

X,Y = meshgrid(xr,yr)

xpass = np.reshape(X,n**2)
ypass = np.reshape(Y,n**2)

array = np.empty((n**2,3))
array[:,0] = xpass
array[:,1] = ypass
array[:,2] = np.zeros(n**2)

#print(array)


B = a.plprofile(array)[2]
B2 = np.reshape(B,(n,n,3))
#print(B2)



plt.figure(2)


mag = np.empty((n,n))

mag[:,:] = np.sqrt(B2[:,:,0]**2+B2[:,:,2]**2+B2[:,:,1]**2)
plt.contour(X,Y, mag, 30,colors='k')



plt.show()


plt.figure(2)

plt.quiver(X,Y, B2[:,:,0], B2[:,:,1], units='inches', #'dots',
           pivot='middle', scale=None, minlength=.01)
axis('equal')


plt.show()

################################################
# X-Z


n = 76
        
xr = linspace(-1.5, 1.5, n)
yr = linspace(-1.5,1.5, n)

X,Y = meshgrid(xr,yr)

xpass = np.reshape(X,n**2)
ypass = np.reshape(Y,n**2)

array = np.empty((n**2,3))
array[:,0] = xpass
array[:,1] = np.zeros(n**2)+2
array[:,2] = ypass

#print(array)


B = a.plprofile(array)[2]
B2 = np.reshape(B,(n,n,3))
#print(B2)



plt.figure(2)


mag = np.empty((n,n))

mag[:,:] = np.sqrt(B2[:,:,0]**2+B2[:,:,2]**2+B2[:,:,1]**2)
plt.contour(X,Y, mag, 30,colors='k')



plt.show()


plt.figure(2)

plt.quiver(X,Y, B2[:,:,0], B2[:,:,2], units='inches', #'dots',
           pivot='middle', scale=None, minlength=.01)
axis('equal')


plt.show()


##########################################################
# Y-Z

n = 76
        
xr = linspace(-1.5, 1.5, n)
yr = linspace(.5,10, n)

X,Y = meshgrid(xr,yr)

xpass = np.reshape(X,n**2)
ypass = np.reshape(Y,n**2)

array = np.empty((n**2,3))
array[:,0] = np.zeros(n**2)
array[:,1] = ypass
array[:,2] = xpass

#print(array)


B = a.plprofile(array)[2]
B2 = np.reshape(B,(n,n,3))
#print(B2)



plt.figure(2)


mag = np.empty((n,n))

mag[:,:] = np.sqrt(B2[:,:,0]**2+B2[:,:,2]**2+B2[:,:,1]**2)
plt.contour(X,Y, mag, 30,colors='k')



plt.show()


plt.figure(2)

plt.quiver(X,Y, B2[:,:,2], B2[:,:,1], units='inches', #'dots',
           pivot='middle', scale=None, minlength=.01)
axis('equal')


plt.show()


