#
# Plot the solar brightness temperature profiles across the solar disk
# for several frequencies
#

import pylab
from pylab import *
import glob

flist = glob.glob('tb_*.txt')

figure();

for fname in flist:
    freq = fname.split('.')
    freq = freq[0].split('_')
    freq = freq[1]
    print freq
    dat = load(fname)
    xax = dat[:,0]
    yax = dat[:,1]/1e6
    plot(xax, yax, linewidth=2); grid(True); hold(True) 

xlabel(r'Distance from center of disk, ($R_{\odot}$)', fontsize=22)        
ylabel(r'Brightness temperature $T_e (10^6 K)$', fontsize=22)

plot([1.3, 1.6], [0.7, 0.7], 'r', linewidth=2);
grid(True); hold(True) 
plot([1.3, 1.6], [0.74, 0.74], 'g', linewidth=2);
grid(True); hold(True) 
plot([1.3, 1.6], [0.78, 0.78], 'b', linewidth=2);
grid(True); hold(True) 

text(1.65, 0.69, '231 MHz', fontsize=22)
text(1.65, 0.73, '202 MHz', fontsize=22)
text(1.65, 0.77, '114 MHz', fontsize=22)
