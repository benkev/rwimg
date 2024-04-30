
#
# Polinomial merging for electron number density
# between the chromosphere and corona
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
A = mat([[1, rh1, rh1**2, rh1**3],
         [0., 1., 2.*rh1, 3.*rh1**2],
         [1, rh2, rh2**2, rh2**3],
         [0., 1., 2.*rh2, 3.*rh2**2]])

b = mat([[5.7e11*exp(-7.7e-4*(R0*(rh1-1.0) - 500.0))],
         [-438.9e6*R0*exp(-7.7e-4*(R0*(rh1-1.0) - 500.0))],
         [1e8*(1.55/rh2**6 + 2.99/rh2**16)],
         [-1e8*(9.3/rh2**7 + 47.84/rh2**17)]])

a = solve(A, b)   # Get the polinomial coefficients
save('poly.txt', a)

print 'a = ', a
print 'A = ', A

print 'A*a = ', A*a
print 'b = ', b

#rh12 = rh1 + (rh2 - rh1)*arange(1001)/1000.
rh12 = linspace(rh1, rh2, 1001)
a = a.A   # Convert to array
p = a[0] + rh12*(a[1] + rh12*(a[2] + rh12*a[3]))
p1 = a[0] + rh1*(a[1] + rh1*(a[2] + rh1*a[3]))
p2 = a[0] + rh2*(a[1] + rh2*(a[2] + rh2*a[3]))



N = where(h < 10000.,
          5.7e11*exp(-7.7e-4*(h-500)),
          1e8*(1.55/rh**6 + 2.99/rh**16)
          )

#den = load('den.txt')


#figure();
#semilogy(rh, N); hold(True); grid(True);
#semilogy(den[:,0], den[:,1]); hold(True); grid(True);
#semilogy(rh12, p); semilogy([rh1, rh2], [p1, p2], 'o')

figure();
plot(rh, N, linewidth=2); hold(True); grid(True);
#plot(den[:,0], den[:,1]); hold(True); grid(True);
plot(rh12, p, linewidth=4); plot([rh1, rh2], [p1, p2], 'o')
xlabel(r'Distance from the Sun center, $r$, in units of $R_{\odot}$',
       fontsize=14)
ylabel(r'Electron number density, $N_e$, cm$^{-3}$', fontsize=14)
text(1.0157,.2e9, r'$r_2$', fontsize=20)
text(1.0129,.2e9, r'$r_1$', fontsize=20)
text(1.012,.9e9, r'$f(r)$', fontsize=20)
text(1.0164, .27e9, r'$g(r)$', fontsize=20)
text(1.0144, .5e9, r'$p(r)$', fontsize=20)

print 'N = ', N

print '\na/a[0] = '
print a/a[0]
