#
# rwimage.py
#
# Plot the result of write_rwimage.c program, the radio wave image of the sun
#

#
# Circle to outline the sun
#
ph = 2*pi*arange(101)/100.0
cx = cos(ph); cy = sin(ph);


pic = load('image.txt');
#pic = load('image_40.0MHz_400x400.txt');
#pic = load('image_80.0MHz_400x400.txt');
#pic = load('image_120.0MHz_400x400.txt');

amax = max(pic[200,:]) 
pic = pic/amax            # scale to 1

ax = (arange(400) - 199.5)/20.0  


fig1 = figure();
imshow(pic, extent=(-10, 10, -10, 10)); gray(); hold('on');
#plot(cx, cy, color=(1.0, 1.0, 0.0));    # Yellow circle to outline the sun disk
plot(cx, cy, color=(1.0, 0.7, 0.0));    # Yellow circle to outline the sun disk
#title('Solar Image at 120.0 MHz')
#title('Solar Image at 80.0 MHz')
title('Solar Image at 40.0 MHz')
xlabel(r'Solar Radii, $R_{\odot}$')
ylabel(r'Solar Radii, $R_{\odot}$') 
hold('off');


fig2 = figure();
#plot([-1,-1], [0, 1], color=(1.0, 0.8, 0.0), linestyle='--');
#plot([1,1], [0, 1], color=(1.0, 0.8, 0.0), linestyle='--');
axes([.2125, .2, .6, .7])
plot([-1,-1], [0, 1], color=(1.0, 0.7, 0.0), linestyle='--'); hold('on');
plot([1,1], [0, 1], color=(1.0, 0.7, 0.0), linestyle='--');
plot(ax, pic[200,:], color='blue'); grid(); #hold('on');
axis([-10, 10, 0, 1.05]);
xlabel(r'Solar Radii, $R_{\odot}$')
ylabel(r'Intensity') 
hold('off');
