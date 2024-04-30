#
# Table the dependance of the solar critical surface radius vs frequency
# for high frequencies, where the Ne distribution is exponential (chromosphere)
#
import numpy
from numpy import *

R0 = 6.950e5            # solar radius, km
q = 4.8e-10;            # Proton charge, StatCoulombs, SGSe
mp = 1.672621636E-24;   # Proton mass, g
me = 9.10938215E-28;    # Electron mass, g

print "f, GHz     rho_cr, g/cm^3   R_crit, sol.rad."
print "--------------------------------------------"

#for i in (1., 2., 5., 10., 20., 50., 100., 200., 500, 1000., 2000., 5000):
for i in xrange(1000):
    f = i*100e6   # Freq, Hz
    rh_cr = pi*mp*me*(f/q)**2
    r = 1.0 + (500.0 - (log(rh_cr) + 27.6787477967)/7.4e-4)/R0
    print "%5.1f        %10.4e       %7f" % (f/1e9, rh_cr, r)

print "--------------------------------------------"

