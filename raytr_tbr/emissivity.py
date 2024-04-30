#
#emissivity.py
#
# Plot the curve of solar plasma emissivity versus rho/rho_critical
#

rh = arange(1.0, 0.0, -0.01)
w = 4*pow(rh,2)*pow((0.5 - rh),2)

plot(rh, w, linewidth=4, color='r')
axis([1, 0, 0, 1])

plot([.25], [w[74]], 'ob') # local maximum at 0.75
grid('on');
title('Plasma Volume Emissivity');
xlabel(r'$\rho / \rho_{cr}$')
ylabel(r'Emissivity, rel. units') 
