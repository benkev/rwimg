
#
# Polinomial "stitching" for electron number density
# between the chromosphere (exponential) and corona (Saito's)
#

import pylab
from pylab import *

def g(r,th):
    g1 = 3.09e8; g2 = 1.58e8; g3 = 2.51e4
    return g1*pow(r,-16)*(1.0 - 0.5*abs(cos(th))) \
           + g2*pow(r,-6)*(1.0 - 0.95*abs(cos(th)))     \
           + g3*pow(r,-2.5)*(1.0 - sqrt(abs(cos(th))))



R0 = 6.950e5   # km, Solar Radius
R = R0*(1. + 5.*arange(100001)/100000.)   # Distance from the sun centre
rh = R/R0   # rh = R/R0
h = R0*(rh - 1)   # Height in chromosphere, 0 .. 10000 km
hc = 10000.0    # Chromosphere height
rhc = hc/R0 + 1.0   # rh at the chromosphere/corona boundary 
print 'rhc = ', rhc
r1 = (hc - 1000.0)/R0 + 1.0   # Merge point in chromosphere
r2 = (hc + 1000.0)/R0 + 1.0   # Merge point in corona

print 'r1, r2 = ', r1, r2

#
# Solve the 3-system of linear equations for the
# polinomial coefficients a[0:3]
#
g1 = 3.09e8; g2 = 1.58e8; g3 = 2.51e4

th = 10.0*pi/180.0

A = mat([[1, r1, r1**2],
         [1, r2, r2**2],
         [0., 1., 2.*r1]])

b = mat([[5.7e11*exp(-7.7e-4*(R0*(r1-1.0) - 500.0))],
         [0.],
         [-438.9e6*R0*exp(-7.7e-4*(R0*(r1-1.0) - 500.0))]])

a = solve(A, b)   # Get the polinomial coefficients
#save('poly.txt', a)

print 'a = ', a
print 'A = ', A

print 'A*a = ', A*a
print 'b = ', b
print 'A*a - b = ', A*a - b

print 'eigvals(A) =\n', eigvals(A)

aa = empty(3, dtype=float)
aa[2] = -(b[0]/(r2-r1) + b[2])/(r2-r1)
aa[1] = b[2] - 2.0*r1*aa[2]
aa[0] = b[0] -r1*aa[1] - r1*r1*aa[2]
print 'aa = ', aa

a = aa



#r12 = r1 + (r2 - r1)*arange(1001)/1000.
r12 = linspace(r1, r2, 1001)
#a = a.A   # Convert to array
p = empty(1001)
for i in xrange(1001):
    p[i] = a[0] + r12[i]*(a[1] + r12[i]*a[2]) \
           + (r12[i]-r1)/(r2-r1)*g(r12[i],th)
p1 = a[0] + r1*(a[1] + r1*a[2])
p2 = a[0] + r2*(a[1] + r2*a[2]) + g(r2,th)



N = where(h < 10000.,
          
          5.7e11*exp(-7.7e-4*(h-500)),
          
            g1*pow(rh,-16)*(1.0 - 0.5*cos(th))     \
          + g2*pow(rh,-6)*(1.0 - 0.95*cos(th))     \
          + g3*pow(rh,-2.5)*(1.0 - sqrt(cos(th)))
          )

#den = load('den.txt')


figure();
semilogy(rh, N); hold(True); grid(True);
#semilogy(den[:,0], den[:,1]); hold(True); grid(True);
semilogy(r12, p); semilogy([r1, r2], [p1, p2], 'o')

figure();
plot(rh, N); hold(True); grid(True);
#plot(den[:,0], den[:,1]); hold(True); grid(True);
plot(r12, p); plot([r1, r2], [p1, p2], 'o')


print 'N = ', N

print '\na/a[0] = '
print a/a[0]
