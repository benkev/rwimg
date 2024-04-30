#
# dens.py
#
# Plot the curves of solar plasma electron number density distribution
# for several heliographic latitudes
# The formula for Ne is taken from the article:
# Kuniji Sato, "A non-spherical asymmetric model of the solar K corona
#   of the minimum type", 1970
#
import pylab
from pylab import *

RhoRsun = 3.3452E-16          # g/cm^3, plasma density at the solar surface
print 'RhoRsun = ', RhoRsun
ProtonMass = 1.6726E-24       # g, Proton mass
print 'ProtonMass = ', ProtonMass
NRsun = RhoRsun/ProtonMass    # cm^(-3), electron # density at solar surface
print 'NRsun = ', NRsun
figure(); hold('on')
plot([1,1], [1,5e10], '--', color=(1,0.8,0), lw=3)
for ph_deg in  (0, 30, 60, 90):
    print 'Lat = ', ph_deg
    r = (100 + arange(1000))/100.0
    ph = ph_deg*pi/180.0   # Heliographic latitude
    Ne = 3.09E8/pow(r,16)*(1 - 0.5*sin(ph)) \
         + 1.58E8/pow(r,6)*(1 - 0.95*sin(ph)) \
         + 0.0251E8/pow(r,2.5)*(1 - sqrt(sin(ph)))
    semilogy(r, Ne, 'g', lw=2)
    
#
# Baumbach's profile
#
Ne = 1e8*(1.55*pow(r,-6) + 2.99*pow(r,-16))
semilogy(r, Ne, 'b', lw=2)
#
#
#
#N1 = NRsun/pow(r,2) # Plasma number density vs solar distance as Ne = NRsun/R^2
#semilogy(r, N1, 'b', lw=2)  
grid(True)
title('Electron Number Density Distributions in Solar Corona', fontsize=18)
xlabel(r'Solar distance, $R_{\odot}$', fontsize=18)
ylabel(r'Electron Number Density, $N_e$', fontsize=18) 
text(11, 2.5, r'$\bf {90^o}$', fontsize=18)
text(11, 250, r'$\bf {60^o}$', fontsize=18)
text(11, 1300, r'$\bf {30^o}$', fontsize=18)
text(11, 4000, r'$\bf {0^o}$', fontsize=18)
text(9.5, 1e4, r'Saito profile', fontsize=18)
text(9.7, 45, r'Baumbach', fontsize=18)
text(9.7, 15, r'profile', fontsize=18)
#text(9, 7e5, r'inverse squares', fontsize=18)

fh = open('Newkirk_density_profile.txt')

fh.readline()
a = loadtxt(fh)
semilogy(a[:,0], a[:,1], 'r', lw=2)
fh.close()
text(2, 2e8, r'Newkirk profile', fontsize=18)


#
# Ne at the chromospheric height of 10000 km
#
R0 = 6.95e5   # Solar radius
r = 1.0 + 10000.0/R0  # Chromosphere height of 10000 km


for ph_deg in  (0, 30, 60, 90):
    ph = ph_deg*pi/180.0   # Heliographic latitude
    
    kn = (3.09E8/pow(r,16)*(1 - 0.5*sin(ph)) \
         + 1.58E8/pow(r,6)*(1 - 0.95*sin(ph)) \
         + 0.0251E8/pow(r,2.5)*(1 - sqrt(sin(ph)))) \
         /(exp(-7.7e-4*(R0*(r-1) - 500.0)))
    print 'kn = ', kn

    Ne_chr = kn*exp(-7.7e-4*(R0*(r-1) - 500.0))
    Ne0 = kn*exp(-7.7e-4*(-500.0))
    print 'Lat = ', ph_deg, ', Ne_chr = ', Ne_chr, ', Ne0 = ', Ne0
    
    Ne_cor = 3.09E8/pow(r,16)*(1 - 0.5*sin(ph)) \
             + 1.58E8/pow(r,6)*(1 - 0.95*sin(ph)) \
             + 0.0251E8/pow(r,2.5)*(1 - sqrt(sin(ph)))
    print 'Lat = ', ph_deg, ', Ne_cor = ', Ne_cor

show();