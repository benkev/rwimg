#
# traj.py
#
# Plot the ray trajectories
#
from pylab import *

cProtonChargeSGSe = 4.8e-10;          # StatCoulombs, SGSe
cProtonMass = 1.672621636E-24;        #g
cElectronMass = 9.10938215E-28;       #g
DensityAtSolarSurface = 3.3452E-16;   #g/cm^3


fh = open("traj.txt", mode='r')

for i in range(4): fh.readline()    # Skip header
RadioFrequency = float(fh.readline())
#
# Calculate the critical density and radius of critical surface
# from the frequency
#
DensityCr = pi*cProtonMass*cElectronMass* \
            pow(RadioFrequency/cProtonChargeSGSe,2);
RadiusCr = sqrt(DensityAtSolarSurface/DensityCr);
print 'RadioFrequency = ', RadioFrequency, ', RadiusCr = ', RadiusCr

for i in xrange(2): fh.readline()    # Skip header
l = (fh.readline()).split()   # Lengths of the trajectories to read
for i in xrange(2): fh.readline()    # Skip header
llt = len(l)
lenTraj = []
for i in xrange(llt): lenTraj.append(int(l[i]))  # Convert lengths to integer

figure()
p = 0; q = lenTraj[0];
for lt in lenTraj:
    tr = empty((lt, 3), dtype=float)
    for i in range(lt):
        l = (fh.readline()).split()
        for j in range(3):
            tr[i,j] = float(l[j])
    plot(tr[:,0], tr[:,2], color=(0., 0., 0.), linewidth=1); hold('on')
grid('on')
axis('equal')

fh.close()

#
# Circle to outline the photosphere
#
ph = 2*pi*arange(101)/100.0
cx = cos(ph); cy = sin(ph);
#
# Circle to outline the chromosphere
#
R0 = 6.950e5   # km, the solar radius
hc = 10000.     # km, the chromosphere height
rchrom = 1.0 + hc/R0   # The chromosphere radius relative to R0
ccx = rchrom*cos(ph); ccy = rchrom*sin(ph);

fill(cx, cy, edgecolor=(1., 0.7, 0.), facecolor=(1., 0.9, 0.))
plot(ccx, ccy, color=(1.0, 0.3, 0.0)); # Circle to outline the chromosphere
#plot(RadiusCr*cx, RadiusCr*cy, 'r', linewidth=1)
xlabel(r'X, Solar Radii, $R_{\odot}$')        
ylabel(r'Z, Solar Radii, $R_{\odot}$')
#
# Annotation
#
title('Electromagnetic Rays Refracting Near Sun at '+ \
      str(RadioFrequency/1e6) +' MHz')     #, size=30)
text(26, -0.2, r'to Earth')
text(-4, -1, r'Sun')
#text(-15, -3.5, r'Critical surface, $R_{cr}=1.05R_{\odot}$')
