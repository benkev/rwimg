
#
# Polinomial "stitching" for electron number density
# between the chromosphere (exponential) and corona (Saito's)
#

import pylab
from pylab import *


R0 = 6.950e5   # km, Solar Radius
R = R0*(1. + 5.*arange(100001)/100000.)   # Distance from the sun centre
rh = R/R0   # rh = R/R0
h = R0*(rh - 1)   # Height in chromosphere, 0 .. 10000 km
hc = 10000.0    # Chromosphere height
rhc = hc/R0 + 1.0   # rh at the chromosphere/corona boundary 
print 'rhc = ', rhc
rh1 = (hc - 1000.0)/R0 + 1.0   # Merge point in chromosphere
rh2 = (hc + 1000.0)/R0 + 1.0   # Merge point in corona

print 'rh1, rh2 = ', rh1, rh2

#
# Solve the 4-system of linear equations for the
# polinomial coefficients a[0:3]
#
g1 = 3.09e8; g2 = 1.58e8; g3 = 2.5e4

th1 = th2 = 30.0*pi/180.0

A = mat([[1, rh1, rh1**2, rh1**3, th1, th1**2],
         [1, rh2, rh2**2, rh2**3, th2, th2**2],
         [0., 1., 2.*rh1, 3.*rh1**2, 0., 0.],
         [0., 1., 2.*rh2, 3.*rh2**2, 0., 0.],
         [0., 0., 0., 0., th1, 2.*th1],
         [0., 0., 0., 0., th2, 2.*th2]])

b = mat([[5.7e11*exp(-7.7e-4*(R0*(rh1-1.0) - 500.0))],
         [g1*pow(rh1,-16)*(1.0 - 0.5*cos(th1)) \
          + g2*pow(rh1,-6)*(1.0 - 0.95*cos(th1)) \
          + g3*pow(rh1,-2.5)*(1.0 - sqrt(cos(th1)))],
         [-438.9e6*R0*exp(-7.7e-4*(R0*(rh1-1.0) - 500.0))],
         [-16.*g1*pow(rh2,-17)*(1.0 - 0.5*cos(th2)) \
          - 6.*g2*pow(rh1,-7)*(1.0 - 0.95*cos(th2)) \
          - 2.5*g3*pow(rh1,-3.5)*(1.0 - sqrt(cos(th2)))],
         [0.],
         [g1*pow(rh1,-16)*0.5*sin(th2) \
          + g2*pow(rh1,-6)*0.95*cos(th2) \
          + g3*pow(rh1,-2.5)*sin(th2)/sqrt(cos(th2))]])

a = solve(A, b)   # Get the polinomial coefficients
#save('poly.txt', a)

print 'a = ', a
print 'A = ', A

print 'A*a = ', A*a
print 'b = ', b

#rh12 = rh1 + (rh2 - rh1)*arange(1001)/1000.
rh12 = linspace(rh1, rh2, 1001)
a = a.A   # Convert to array
p = a[0] + rh12*(a[1] + rh12*(a[2] + rh12*a[3])) + th1*(a[4] + th1*a[5])
p1 = a[0] + rh1*(a[1] + rh1*(a[2] + rh1*a[3])) + th1*(a[4] + th1*a[5])
p2 = a[0] + rh2*(a[1] + rh2*(a[2] + rh2*a[3])) + th2*(a[4] + th2*a[5])



N = where(h < 10000.,
          5.7e11*exp(-7.7e-4*(h-500)),
          1e8*(1.55/rh**6 + 2.99/rh**16)
          )

den = load('den.txt')


figure();
semilogy(rh, N); hold(True); grid(True);
semilogy(den[:,0], den[:,1]); hold(True); grid(True);
semilogy(rh12, p); semilogy([rh1, rh2], [p1, p2], 'o')

figure();
plot(rh, N); hold(True); grid(True);
plot(den[:,0], den[:,1]); hold(True); grid(True);
plot(rh12, p); plot([rh1, rh2], [p1, p2], 'o')


print 'N = ', N

print '\na/a[0] = '
print a/a[0]
