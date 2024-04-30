from pylab import *   # fsolve, trig funcs etc.
import scipy
from scipy.optimize import *   # fsolve, trig funcs etc.

def saito(r, t, v):
    """
    Kuniji Saito electron number density distribution
    (r,t) are spherical coordinates (radius, theta)
    v is a specified value of the Saito function.
    """
    from scipy import pi, cos, sqrt   # fsolve, trig funcs etc.
    g1 = 3.09e8; g2 = 1.58e8; g3 = 2.51e6;
    #r = x[0]; t = x[1];
    s = g1*pow(r,-16)*(1.0 - 0.5*cos(t)) \
           + g2*pow(r,-6)*(1.0 - 0.95*cos(t)) \
           + g3*pow(r,-2.5)*(1.0 - sqrt(cos(t)))
    return s - v

#print saito(2, 0);
#print saito(2, pi/6);
#print saito(2, pi/4);
#print saito(2, pi/3);
#print saito(2, pi/2);

#t = 0 #pi/3.0
#r = fsolve(saito, 1.0, args=(t, 114.7e-2))

#print 'At r = ', r, ' saito is ', saito(r, t, 814114.7)
#
# Circle to fill
#
ph = 2*pi*arange(101)/100.0
cx = cos(ph); cy = sin(ph);
x0, x1, y0, y1 = -10, 10, -10, 10

np = 101
ths = linspace(0., pi/2., np);
r = empty(np, dtype=float)
#Nes = [1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7, 1e8]
Nes = [1e4, 1e5, 1e6, 1e7, 1e8]


#figure(figsize=(12.5,6.5)); hold(True); axis('equal');
figure(figsize=(12.5,6.5)); hold(True); axis('equal');
for v in Nes:
    for i in xrange(np):
        r[i] = fsolve(saito, 1.0, args=(ths[i], v))
        #print 90. - ths[i]*180/pi, r[i], saito(r[i], ths[i], v)
    plot(r*sin(ths), r*cos(ths), 'b', linewidth=2);  
    plot(-r*sin(ths), r*cos(ths), 'b', linewidth=2);  
    plot(r*sin(ths), -r*cos(ths), 'b', linewidth=2);  
    plot(-r*sin(ths), -r*cos(ths), 'b', linewidth=2);  #grid(True)  

#r1 = r[::-1]  # Reverse order
#r[:np] = r1[:np]

#fill(cx, cy, color='#ffff00', facecolor='#ffff00')
fill(cx, cy, ec='#ff0000', fc='#ffff00', lw=2);

plot((-12,-10), (0,0), 'k', linewidth=2)
plot((12,10), (0,0), 'k', linewidth=2)
plot((0,0), (5,6.5), 'k', linewidth=2)
plot((0,0), (-4.5,-6), 'k', linewidth=2)

title(r'Coronal Electron Number Density, cm$^{-3}$', fontsize=16)
xlabel(r'Solar Radii, $R_{\odot}$', fontsize=16)
ylabel(r'Solar Radii, $R_{\odot}$', fontsize=16) 

text(9.9, 0.2, 'equator', fontsize=16)
text(-0.5, 4, 'pole', fontsize=16)
text(-0.5, -0.25, 'SUN', fontsize=16)
text(7, 1.6, r'$10^4$', fontsize=16)
text(3.5, 1.1, r'$10^5$', fontsize=16)
text(2.1, 0.6, r'$10^6$', fontsize=16)
text(1.5, 0.2, r'$10^7$', fontsize=8)
text(-1.5, 6.7, r'Saito Model', fontsize=18)
