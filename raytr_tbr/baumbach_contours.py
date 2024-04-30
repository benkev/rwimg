from pylab import *   # fsolve, trig funcs etc.
import scipy
from scipy.optimize import *   # fsolve, trig funcs etc.

def baumbach(r, v):
    """
    Baumbach electron number density distribution
    """
    b = 1e8*(1.55*pow(r,-6) + 2.99*pow(r,-16))
    return b - v
    
ph = 2*pi*arange(501)/500.0
cx = cos(ph); cy = sin(ph);
#Nes = [1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7, 1e8]
#Nes = [5e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7, 1e8]
Nes = [1e4, 1e5, 1e6, 1e7, 1e8]

figure(figsize=(6.5,6.5)); axis('equal');  hold(True);

for v in Nes:
        r = fsolve(baumbach, 1.0, args=(v))
        plot(r*cx, r*cy, 'b', linewidth=2)

fill(cx, cy, ec='#ff0000', fc='#ffff00', lw=2);


plot((-5.5,-6.5), (0,0), 'k', linewidth=2)
plot((5.5,6.5), (0,0), 'k', linewidth=2)
plot((0,0), (5.5,6.5), 'k', linewidth=2)
plot((0,0), (-5.5,-6.4), 'k', linewidth=2)

title(r'Coronal Electron Number Density, cm$^{-3}$', fontsize=16)
xlabel(r'Solar Radii, $R_{\odot}$', fontsize=16)
ylabel(r'Solar Radii, $R_{\odot}$', fontsize=16) 


text(4.5, 0.2, 'equator', fontsize=16)
text(-0.6, 4, 'pole', fontsize=16)
text(-0.5, -0.25, 'SUN', fontsize=16)
text(4.5, 2.05, r'$10^4$', fontsize=16)
text(3.15, 1.25, r'$10^5$', fontsize=16)
text(2.2, 0.7, r'$10^6$', fontsize=14)
text(1.55, 0.25, r'$10^7$', fontsize=8)
text(-3, 6.7, r'Baumbach Model', fontsize=18)
