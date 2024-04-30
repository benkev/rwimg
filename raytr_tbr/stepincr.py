import pylab 
from pylab import *

x = .10;
y = empty(26);
y[0] = .25;
for i in range(25):
	y[i+1] = (2 - y[i]/x)*y[i]
print y
f1 = figure();
plot(y, linewidth=2); hold(True); grid(True)
title(r'Solutions to $y_{i+1} = (2 - y_{i}/x_i) y_i$')
xlabel(r'i')
ylabel(r'$y_i$') 

#y[0] = .01;
#for i in range(25):
#	y[i+1] = (2 - y[i]/x)*y[i]
#print y
#plot(y, color='k', linewidth=2); #hold(True);
#plot(y, linewidth=2); #hold(True);

#y[0] = .0001;
#for i in range(25):
#	y[i+1] = (2 - y[i]/x)*y[i]



#print y
#plot(y, color='k', linewidth=2); #hold('on');
#plot(y, linewidth=2); #hold('on');

#y[0] = .000001;
#for i in range(25):
#	y[i+1] = (2 - y[i]/x)*y[i]
#print y
#plot(y, color='k', linewidth=2)
#plot(y, linewidth=2)
