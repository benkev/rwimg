#
# rwimage.py
#
# Plot the result of write_rwimage.c program, the radio wave image of the sun
#
from pylab import *

#N = 200    # Size of the picture in both dimensions

#
# Circle to outline the photosphere
#
ph = 2*pi*arange(101)/100.0
cx = cos(ph); cy = sin(ph);

fh = open("image.txt", mode='r')

for i in range(4): fh.readline()    # Skip header
Freq = float(fh.readline())
print 'RadioFrequency = ', Freq
for i in xrange(4): fh.readline()    # Skip header
(x0, y0, x1, y1) = (fh.readline()).split()  # Image range as (x0, y0, x1, y1)
(x0, y0, x1, y1) = map(float, (x0, y0, x1, y1))
for i in xrange(2): fh.readline()    # Skip header
(Ny, Nx) = (fh.readline()).split()  # Image Y and X sizes in pixels
(Ny, Nx) = map(int, (Ny, Nx))
for i in xrange(2): fh.readline()    # Skip header
pic = empty((Ny, Nx), dtype=float32)
for i in xrange(Ny):
    pic[i,:] = map(float, (fh.readline()).split())
#imshow(pic, extent=(x0, x1, y0, y1)); gray(); hold('on');

amax = max(pic[Ny/2,:])
print '\nmax value = ', amax, '\n'
#save('xsec.txt', pic[N/2,:])
#pic = pic/amax            # scale to 1
pic = pic/1e6           # scale to 1 = 1e6 K

ax = y0 + (y1-y0)*arange(Ny, dtype=float)/float(Ny-1) # Plot x-axix

fig1 = figure();
imshow(pic, extent=(x0, x1, y0, y1)); gray(); hold('on');
plot(cx, cy, color=(1.0, 1.0, 0.0));    # Yellow circle to outline the sun disk
plot(cx, cy, color=(1.0, 0.7, 0.0));    # Yellow circle to outline the sun disk
title('Solar Image at '+str(int(Freq/1e6))+' MHz')
xlabel(r'Solar Radii, $R_{\odot}$')
ylabel(r'Solar Radii, $R_{\odot}$') 
hold('off');

#
# 32T model
#
#dm66 = load('model066_slice.txt')
#dm181 = load('model181_slice.txt')

#dm181[:,2] /= max(dm181[:,2])   # Normalize so that positive max = 1


fig2 = figure();
plot([-1,-1], [-0.1, 1.1], color=(1.0, 0.7, 0.0), \
     linewidth=3, linestyle='--'); hold('on');
plot([1,1], [-0.1, 1.1], color=(1.0, 0.7, 0.0), \
     linewidth=3, linestyle='--');
#plot(dm181[:,0]/8.-64.3, dm181[:,2], linewidth=3); grid(True)
#plot(ax, pic[N/2,:], linewidth=3, color='blue'); grid(True); #hold('on');
plot(ax, pic[:,Ny/2], linewidth=3, color='green'); grid(True); #hold('on');
#plot(ax, pic[N/2,:]+0.1, 'o'); grid(True); #hold('on');
#axis([-10, 10, 0, 1.05]);
title('Solar Brightness Temperature Profile at '+str(int(Freq/1e6))+' MHz')
xlabel(r'Solar Radii, $R_{\odot}$')
#ylabel(r'Intensity') 
ylabel(r'Brightness Temperature, $T_e (10^6 K)$') 
hold('off');

for i in xrange(Nx):
    print pic[i,Ny/2]


show()
