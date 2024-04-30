
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
#R = R0*(1. + 5.*arange(100001)/100000.)   # Distance from the sun centre
R = R0*logspace(0., 0.3, 100001)   # Distance from the sun centre
rh = R/R0   # rh = R/R0
h = R0*(rh - 1)   # Height in chromosphere, 0 .. 10000 km
hc = 10000.0    # Chromosphere height
rhc = hc/R0 + 1.0   # rh at the chromosphere/corona boundary 
print 'rhc = ', rhc
r1 = (hc - 1000.0)/R0 + 1.0   # Merge point in chromosphere
r2 = (hc + 1000.0)/R0 + 1.0   # Merge point in corona
rdisrup = hc/R0 + 1.0   # Disruption point
print 'r1, r2 = ', r1, r2

#
# Solve the 4-systems of linear equations for the
# polinomial coefficients a[0:3] and b[0:3] 
#
g1 = 3.09e8; g2 = 1.58e8; g3 = 2.51e4

th = 0.0*pi/180.0

A = mat([[1., r1, r1**2, r1**3],
         [1., r2, r2**2, r2**3],
         [0., 1., 2.*r1, 3.*r1**2],
         [0., 1., 2.*r2, 3.*r2**2]]);

print 'A = ', A

rhsa = mat([[5.7e11*exp(-7.7e-4*(R0*(r1-1.0) - 500.0))],
            [0.],
            [0.],
            [0.]]);
print 'rhsa = ', rhsa

#[-438.9e6*R0*exp(-7.7e-4*(R0*(r1-1.0) - 500.0))]])
am = solve(A, rhsa)   # Get the a(r) polinomial coefficients
a = empty(4,float) 
for i in xrange(4): a[i] = am[i,0]
print 'a = ', a

gr2th =   g1*pow(r2,-16)*(1.0 - 0.5*cos(th))  \
        + g2*pow(r2,-6)*(1.0 - 0.95*cos(th))  \
        + g3*pow(r2,-2.5)*(1.0 - sqrt(cos(th)));


dfdr1 = -438.9e6*R0*exp(-7.7e-4*(R0*(r1-1.0) - 500.0));   # df(r2)/dr

dgdr2 =   - 16.*g1*pow(r2,-17)*(1.0 - 0.5*cos(th))     \
          -  6.*g2*pow(r2, -7)*(1.0 - 0.95*cos(th))     \
          - 2.5*g3*pow(r2,-3.5)*(1.0 - sqrt(cos(th)))


rhsb = mat([[0.],
            [1.],
            [dfdr1/gr2th],
            [dgdr2/gr2th]]);

print 'rhsb = ', rhsb

bm = solve(A, rhsb)   # Get the a(r) polinomial coefficients
b = empty(4,float) 
for i in xrange(4): b[i] = bm[i,0]
print 'b = ', b

r12 = linspace(r1, r2, 1001)
#a = a.A   # Convert to array
p = empty(1001)
for i in xrange(1001):
    r = r12[i];
    p[i] = a[0] + r*(a[1] + r*(a[2] + r*a[3])) \
                     + (b[0] + r*(b[1] + r*(b[2] + r*b[3])))*g(r12[i],th)

#p1 =    a[0] + r1*(a[1] + r1*(a[2] + r1*a[3])) is enough!
p1 =    a[0] + r1*(a[1] + r1*(a[2] + r1*a[3])) \
     + (b[0] + r1*(b[1] + r1*(b[2] + r1*b[3])))*g(r1,th)


#p2 =    a[0] + r2*(a[1] + r2*(a[2] + r2*a[3])) + g(r2,th) is enough!
p2 =    a[0] + r2*(a[1] + r2*(a[2] + r2*a[3])) \
     + (b[0] + r2*(b[1] + r2*(b[2] + r2*b[3])))*g(r2,th)


N = where(h < 10000.,
          
          5.7e11*exp(-7.7e-4*(h-500)),
          
            g1*pow(rh,-16)*(1.0 - 0.5*cos(th))     \
          + g2*pow(rh,-6)*(1.0 - 0.95*cos(th))     \
          + g3*pow(rh,-2.5)*(1.0 - sqrt(cos(th)))
          )

#den = load('den.txt')


#figure();
#semilogy(rh, N); hold(True); grid(True);
#semilogy(den[:,0], den[:,1]); hold(True); grid(True);
#semilogy(r12, p); semilogy([r1, r2], [p1, p2], 'o')

#figure();
plot(rh, N, 'b', linewidth=2); hold(True); grid(True);

plot(r12, p, 'g', linewidth=4); plot([r1, r2], [p1, p2], 'ro')

xlabel(r'Distance from the Sun center, $r$, in units of $R_{\odot}$',
       fontsize=14)
ylabel(r'Electron number density, $N_e$, cm$^{-3}$', fontsize=14)
text(1.0157, 2e7, r'$r_2$', fontsize=20)
text(1.0129, 2e7, r'$r_1$', fontsize=20)
text(1.012,.9e9, r'$f(r)$', fontsize=20)
text(1.0162, .27e9, r'$g(r,\theta)$', fontsize=20)
text(1.0144, .6e9, r'$p(r,\theta)$', fontsize=20)
text(1.0168, 1.25e8, r'$\theta=0^{\circ}$', fontsize=20)
text(1.0168, 2.0e8, r'$\theta=45^{\circ}$', fontsize=20)
text(1.0168, 3.8e8, r'$\theta=90^{\circ}$', fontsize=20)

print 'N = ', N

print '\na/a[0] = '
print a/a[0]


#dgdth =   0.50*g1*pow(r2,-16)*sin(th)  \                  # dg(r2,th)/dth
#        + 0.95*g2*pow(r2,-6)*sin(th)) \
#        + 0.50*g3*pow(r2,-2.5)*sin(th)/sqrt(cos(th)));


