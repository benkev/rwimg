#
# Newkirk's streamer model
#

from pylab import *

a = 1.6498
b = -14.322
c = -0.6925
d = 1.0
e = 41.4
f = 52.24
g = 25.06

r = linspace(1,2.5,1000)

#s = (c*r + d)/((((a*r + b)*r + e)*r + f)*r + g)
s = empty(1000, float)
#for i in xrange(1000): s1[i] = (a*r[i] + b)/polyval([c,d,e,f,g], r[i])
s = 0.15 + 0.15*exp(-(r-1.5)**2/.22)
#figure(1);
plot(r, s, linewidth=2); grid(True); hold(True);

r2 = [1., 1.125, 1.275, 1.375, 1.5, 1.75, 2.]
s2 = [.2, .235, .269, .31, .27, .256, .197]

plot(r2, s2, 'or'); grid(True)
