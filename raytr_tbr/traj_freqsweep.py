#!/usr/bin/python
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

flist = ['traj_10.txt', 'traj_18.txt', 'traj_40.txt', 'traj_80.txt', \
         'traj_200.txt', 'traj_3000.txt']
#         'traj_400.txt'] #, 'traj_800.txt', 'traj_3000.txt']
figure()
ax = subplot(111)

for fl in flist:
    fh = open('/home/benkev/raytr_tbr/Tb_Saito/'+fl, mode='r')

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
    llt = len(l)  # How many trajectories
    lenTraj = []
    for i in xrange(llt): lenTraj.append(int(l[i]))  # Convert lengths to int

    tcnt = 0
    for lt in lenTraj:
        tcnt += 1
        tr = empty((lt, 3), dtype=float)
        for i in range(lt):
            l = (fh.readline()).split()
            for j in range(3):
                tr[i,j] = float(l[j])
        if tcnt == 70:
            plot(tr[:,0], tr[:,2], color=(1., 0., 0.), linewidth=1);
            #break
        if tcnt == 95:    
            plot(tr[:,0], tr[:,2], color=(0., 0., 1.), linewidth=1);
            break
    hold('on')
    fh.close()
    
grid('on')
#xlim(0.925, 2.430); ylim(-0.015, 1.151)
#axis('equal')
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
xlim(0.9, 2.430); ylim(-0.015, 1.151)
axis('equal')
#plot(RadiusCr*cx, RadiusCr*cy, 'r', linewidth=1)
xlabel(r'X, Solar Radii, $R_{\odot}$', fontsize=18)        
ylabel(r'Z, Solar Radii, $R_{\odot}$', fontsize=18)
#
# Annotation
#
title('Refraction Near Sun for Several Frequencies', fontsize=18)
text(26, -0.2, r'to Earth')
text(0.65, 0.1, r'Sun', fontsize=24)
#text(-15, -3.5, r'Critical surface, $R_{cr}=1.05R_{\odot}$')
text(2.25, 0.32, '10 MHz')
text(1.85, 0.32, '18 MHz')
text(1.45, 0.26, '40 MHz')
text(1.20, 0.22, '80 MHz')
text(1.66, 0.59, '200', rotation=25)
text(1.62, 0.46, '3GHz', rotation=20)

text(1.98, 1.52, '10 MHz')
text(1.45, 1.52, '18 MHz')
text(0.86, 1.52, '40 MHz')
text(0.61, 1.52, '80')
text(0.54, 1.03, '200MHz', rotation=100)
text(0.71, 1.25, '3GHz', rotation=87)

ax.annotate("", xy=(2.76, 0.68), xytext=(2.76, 0.0),
            arrowprops=dict(arrowstyle="<->", connectionstyle="arc3"))
text(2.535, 0.03, '0.18', rotation=90)

ax.annotate("", xy=(2.60, 0.19), xytext=(2.60, 0.0),
            arrowprops=dict(arrowstyle="<->", connectionstyle="arc3"))
text(2.705, 0.3, '0.67', rotation=90)

show();


