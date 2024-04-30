#
# Brightness temperature of the solar disk calculation
# for the central ray by integration
#

import pylab
from pylab import *
import sys


def densNe(r):
    """
    Coronal and chromospheric plasma electron number density, cm^-3
    by Baumbach-Allen's formulae
    r: distance from the centre of the sun in solar radii
    """
    R0 = 6.950e5;   # km, Solar radius
    #
    # Polynomial coefs to smoothly merge exponential and power laws of
    # Ne distributions
    #
    a0 = 1.686846091695834800e+16
    a1 = -4.981085198383761600e+16
    a2 =  4.902880054853674400e+16 
    a3 = -1.608634376506040800e+16

    # r_chromo = 1 + (h_chromo - 1000)/R0,
    # r_corona = 1 + (h_chromo + 1000)/R0,
    r_chromo = 1.01294964029;  # *R0, km
    r_corona = 1.01582733813;  # *R0, km
    
    if r < r_chromo:
        expo = exp(-7.7e-4*(R0*(r-1) - 500.0));
        Ne = 5.7e11*expo;
        #dNdr = -438.9e6*R0*expo/r;
        #print 'r, expo, Ne =', r, expo, Ne
        
        
    elif r < r_corona: #  but r >= r_chromo
        Ne = a0 + r*(a1 + r*(a2 + r*a3));
        # dNdr = (1/r)*dNe/dr
        #dNdr = a1/r + 2.0*a2 + 3*a3*r;
        
    else: # r >= r_corona:
        r2 = pow(r,2)
        rm2 = 1/r2;                  # r^(-2)
        rm6 = pow(r,-6);             # r^(-6)
        rm8 = rm6*rm2;               # r^(-8)
        rm16 = pow(r,-16);           # r^(-16)
        rm18 = rm16*rm2;             # r^(-18)
        Ne = 1e8*(1.55*rm6 + 2.99*rm16);
        # dNdr = (1/r)*dNe/dr
        #dNdr = -1e8*(9.3*rm8 + 47.84*rm18);

    return Ne



np = 101
Freq = 18e6
R0 = 6.950e5;   # km, Solar radius
R0cm = 6.95E10;   # cm
r_chromo = 1.01294964029;  # *R0, km
r_corona = 1.01582733813; # *R0, km
Te_corona = 1.0e6; # K
Te_chromo = 3.0e4; # K
ProtonChargeSGSe = 4.8e-10;      # StatCoulomb, SGSe
cProtonMass = 1.672621636E-24;   # g
cElectronMass = 9.10938215E-28;  # g


#DensityCr = cPi*cProtonMass*cElectronMass*pow(Freq/ProtonChargeSGSe,2);
Necr = pi*cElectronMass*pow(Freq/ProtonChargeSGSe,2) # Critical Ne
#r = 1.0 + arange(np,dtype=float)/float(np-1)
#r = 1.0 + logspace(-3, 1, np)
b2 = 5.0 - 1.0 - r_corona   # Distance 2..25 with DS = 0.1
b1 = 1.0                     # Distance r_corona..2 with DS = 0.01
b0 = r_corona - 1.0         # Distance 1..r_corona with DS = 0.001
n2 = int(b2/0.01)
n1 = int(b1/0.0001)
n0 = int(b0/0.000001)
np = n0 + n1 + n2
r = empty(np, dtype=float)
Ne = empty(np, dtype=float)
print 'b0, b1, b2 = ',  b0, b1, b2
print 'n0, n1, n2 = ',  n0, n1, n2
print 'np = ', np
print len(r[0:n0]), len(arange(n0,dtype=float)/n0*b0)
r[:n0] = 1.0 + arange(n0,dtype=float)/n0*b0
print len(r[n0:n0+n1]), len(arange(n1,dtype=float)/n0*b0)
r[n0:n0+n1] = r_corona + arange(n1,dtype=float)/n1*b1
print len(r[n0+n1:]), len(arange(n2,dtype=float)/n2*b2)
r[n0+n1:] = r_corona + 1.0 + arange(n2,dtype=float)/n2*b2

#plot(r)
#sys.exit(0)

CAbsorp = 0.0
OpDepth = 0.0
Tb = 0.0
r_pred = 25.0

#for i in xrange(np-1, -1, -1):
#    DS = r_pred - r[i]
for i in xrange(np):
    DS = r[i] - r_pred
    Ne[i] = densNe(r[i])
    Eps = 1.0 - Ne[i]/Necr

    
    if r[i] > r_chromo:
	Te = Te_corona;
	xi = 2.0; # cm^6 K^1.5 Hz^2
    else:
	Te = Te_chromo;
	xi = 1.4; # cm^6 K^1.5 Hz^2
        
    if Eps > 0.0:
        if i%100 == 0:
            print 'r, DS, CAbsorp, OpDepth, Tb = ', \
                  "%5.3f %6.4f %9.3e %10.7f" % (r[i], DS, CAbsorp, OpDepth), Tb
        CAbsorp = xi*pow(Ne[i],2)/(pow(Te,1.5)*pow(Freq,2)*sqrt(Eps))
        dOpDepth = CAbsorp*R0cm*DS
        OpDepth += dOpDepth
        dTb =  Te*exp(OpDepth)*dOpDepth
        Tb += dTb
        
    r_pred = r[i]

print '\nOpDepth = ', OpDepth, ', Tb_before = ', Tb
Tb = exp(-OpDepth)*Tb

print '\nOpDepth = ', OpDepth, ', Tb = ', Tb
print 'exp(-OpDepth) = ', exp(-OpDepth)
#grid(True)
#semilogy(r,Ne)
#print 'r, Ne = \n'
#for i in xrange(np): print r[i], Ne[i]





