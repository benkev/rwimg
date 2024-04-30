#
# Read brightness temperature data
#

from pylab import *

fh = open("tbr.txt", mode='r')

for i in range(4): fh.readline()    # Skip header
RadioFrequency = float(fh.readline())

for i in xrange(2): fh.readline()    # Skip header
l = (fh.readline()).split()   # Lengths of the trajectories to read
for i in xrange(2): fh.readline()    # Skip header
llt = len(l)
lenTrajTb = []
for i in xrange(llt): lenTrajTb.append(int(l[i]))  # Convert lengths to integer

#nr = int(raw_input("Ray number?: "))
figure()
p = 0; q = lenTrajTb[0];
for lt in lenTrajTb:
    tb = empty(lt, dtype=float)
    for i in range(lt):
        l = fh.readline()
        tb[i] = float(l)
    semilogy(tb, color=(0., 0., 0.), linewidth=2); hold('on')
grid('on')

fh.close()


