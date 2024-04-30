#
# Plot the solar brightness temperature images
# for several frequencies
#

import pylab
from pylab import *
import glob
#import os

#
# Circle to outline the photosphere
#
ph = 2*pi*arange(101)/100.0
cx = cos(ph); cy = sin(ph);

#cdir = os.getcwd()
#os.chdir('Tb_Smerd')
flist = glob.glob('/home/benkev/raytr_tbr/image_*.txt')
#freqs = [40] #, 100, 200, 400, 600, 1200, 3000, 30000];
freqs = [100];
#leg = []

#
# Leave only filenames
#
nfm = len(flist)
nfm1 = len(flist) - 1
for i in xrange(nfm):
    flist[i] = flist[i].split('/')[-1]
print flist
#
# Sort filenames in frequency ascending orders
#

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

#for fname in flist: print fname
#figure(1,figsize=(8,8))
#figure(2,figsize=(8,8))
fn = 1

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
        cFreq = str(Freq)+' MHz'
    else:
        cFreq = str(Freq/1e3)+' GHz'
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
    pic = pic/1e6
    #figure(fn); fn = fn + 1;
    figure(); fn = fn + 1;
    #subplot(121);
    print 'PLOTTING ', cFreq
    #imshow(pic, extent=(x0, x1, y0, y1), cmap=cm.gist_heat); hold('on');
    contourf(flipud(pic), 25, extent=(x0, x1, y0, y1), cmap=cm.gist_heat);
    hold('on');
    #plot(cx, cy, color=(1.0, 1.0, 0.0)); # Circle to outline the sun disk
    plot(cx, cy, color='y'); # Circle to outline the sun disk
    #subplot(122);
    #imshow(pic, extent=(x0, x1, y0, y1), cmap=cm.hot); hold('on');
    colorbar() 
    title(r'Helmet Streamers in $T_b$, MK at '
    #title(r'Sun in Brightness Temperature, $T_b (10^6 K)$ at '
          + cFreq, fontsize=16)
    #      +str(int(Freq))+' MHz', fontsize=16)
    xlabel(r'Solar Radii, $R_{\odot}$', fontsize=16)
    ylabel(r'Solar Radii, $R_{\odot}$', fontsize=16) 
    hold('off');

#os.chdir(cdir)
#os.chdir('..')

text(1.2, 1.3, r'$+4 N_e$', fontsize=16, color='w')
text(1.2, -1.3, r'$+2 N_e$', fontsize=16, color='w')
text(-1.7, -1.3, r'$+1 N_e$', fontsize=16, color='w')
text(-1.7, 1.3, r'$+0.5 N_e$', fontsize=16, color='w')

show()
