#
# dist.py
#
# Plot the curves of solar plasma electron number density distribution
# for several heliographic latitudes
# The formula for Ne is taken from the article:
# Kuniji Sato, "A non-spherical asymmetric model of the solar K corona of the minimum type", 1970
#
RhoRsun = 3.3452E-16          # g/cm^3, plasma density at the solar surface
print 'RhoRsun = ', RhoRsun
ProtonMass = 1.6726E-24       # g, Proton mass
print 'ProtonMass = ', ProtonMass
NRsun = RhoRsun/ProtonMass    # cm^(-3), electron number density at the solar surface
print 'NRsun = ', NRsun
figure(); hold('on')
plot([1,1], [1,5e10], '--', color=(1,0.8,0), lw=3)
for ph_deg in  (0, 30, 60, 90):
    print 'Lat = ', ph_deg
    r = (100 + arange(1000))/100.0
    ph = ph_deg*pi/180.0   # Heliographic latitude
    Ne = 3.09E8/pow(r,16)*(1 - 0.5*sin(ph)) + 1.58E8/pow(r,6)*(1 - 0.95*sin(ph)) \
         + 0.0251E8/pow(r,2.5)*(1 - sqrt(sin(ph)))
    semilogy(r, Ne, 'g', lw=2)  

#
#
#
N1 = NRsun/pow(r,2) # Plasma number density vs solar distance as Ne = NRsun/R^2
semilogy(r, N1, 'b', lw=2)  
grid()
title('Electron Number Density Distribution in Corona by Kuniji Saito')
xlabel(r'Solar distance, $R_{\odot}$')
ylabel(r'Electron Number Density, $N_e$') 
text(11, 2.5, r'$\bf {90^o}$')
text(11, 250, r'$\bf {60^o}$')
text(11, 1300, r'$\bf {30^o}$')
text(11, 4000, r'$\bf {0^o}$')
text(7, 2000, r'Saito profile')
text(9, 7e5, r'inverse squares')

fh = open('Newkirk_density_profile.txt')

fh.readline()
a = load(fh)
semilogy(a[:,0], a[:,1], 'r', lw=2)
fh.close()
text(2, 2e8, r'Newkirk profile')
