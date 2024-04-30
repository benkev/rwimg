#
# Plot the solar brightness temperature profiles across the solar disk
# for several frequencies
#

import pylab
from pylab import *
import glob
import os

cdir = os.getcwd()
#os.chdir('Tb_Saito')
os.chdir('Tb_Smerd')
flist = glob.glob('image_*.txt')

#
# Sort filenames in frequency ascending orders
#
nfm1 = len(flist) - 1
sorted = True
f1 = ''
while True:
    sorted = True
    for i in xrange(nfm1):
        freqi = int((flist[i].split('_')[1]).split('.')[0])
        freqi1 = int((flist[i+1].split('_')[1]).split('.')[0])
        if freqi1 < freqi:
            sorted = False
            f1 = flist[i+1]
            flist[i+1] = flist[i]
            flist[i] = f1
    if sorted: break

for f in flist: print f

freqs = [60, 100, 200, 400, 600, 1200, 3000, 30000];
leg = []

#for fname in flist: print fname
figure(1,figsize=(8,8))
figure(2,figsize=(8,8))

for fname in flist:
    print fname
    freq = fname.split('.')
    freq = freq[0].split('_')
    freq = int(freq[1])
    print freq
    if freq not in freqs: continue
    print 'freq = ', freq, ' MHz'
    fh = open(fname, mode='r')
    for i in range(4): fh.readline()    # Skip header
    Freq = float(fh.readline())/1e6
    print 'RadioFrequency = ', Freq
    if Freq < 1000.0:
        leg.append(str(Freq)+' MHz')
    else:
        leg.append(str(Freq/1e3)+' GHz')
    for i in xrange(4): fh.readline()    # Skip header
    (x0, y0, x1, y1) = (fh.readline()).split()  # Im. range as (x0, y0, x1, y1)
    (x0, y0, x1, y1) = map(float, (x0, y0, x1, y1))
    for i in xrange(2): fh.readline()    # Skip header
    (Ny, Nx) = (fh.readline()).split()  # Image Y and X sizes in pixels
    (Ny, Nx) = map(int, (Ny, Nx))
    for i in xrange(2): fh.readline()    # Skip header
    pic = empty((Ny, Nx), dtype=float32)
    for i in xrange(Ny):
        pic[i,:] = map(float, (fh.readline()).split())
    Tby = pic[Ny/2:,Nx/2]/1e6
    Tbx = pic[Ny/2,Nx/2:]/1e6
    axx = x0 + (x1-x0)*arange(Nx, dtype=float)/float(Nx-1) # Plot x-axix
    axx = axx[Ny/2:]
    axy = y0 + (y1-y0)*arange(Ny, dtype=float)/float(Ny-1) # Plot x-axix
    axy = axy[Ny/2:]
    figure(1);
    plot(axy, Tby, linewidth=2); grid(True); hold(True) 
    figure(2);
    plot(axx, Tbx, linewidth=2); grid(True); hold(True) 
    #xax = pic[:,0]
    #yax = pic[:,1]/1e6
    #plot(xax, yax, linewidth=2); grid(True); hold(True) 
    #fig1 = figure();
    #imshow(pic, extent=(x0, x1, y0, y1)); gray(); hold('on');

os.chdir(cdir)


figure(1);
title('Soilar Brightness Temperature in South-North Direction', fontsize=18) 
xlabel(r'Distance from center of disk, ($R_{\odot}$)', fontsize=18)        
ylabel(r'Brightness temperature $T_e (10^6 K)$', fontsize=18)
legend(leg)

figure(2);    
title('Soilar Brightness Temperature in West-East Direction', fontsize=18) 
xlabel(r'Distance from center of disk, ($R_{\odot}$)', fontsize=18)        
ylabel(r'Brightness temperature $T_b (10^6 K)$', fontsize=18)
legend(leg)

