from pylab import *

L = 100.0;


#a = loadtxt('ray0.txt'); plot(a[:,0], a[:,1]); hold('on');

#al = 0.1*pi/180;
#x0 = 0; y0 = 0;
#y = arange(10001)/10000.*a[size(a[:,0])-1,1];
#x = x0 + L*pow(cos(al),2) - L*pow((cos(al) - 1./(2.0*L*sin(al))*(y - y0)), 2);
#plot(x,y, 'k'); grid('on');

a = loadtxt('ray1.txt'); 

al = 1.0*pi/180; #plot(a[:,0], a[:,1]); hold('on');
x0 = 0; y0 = 0;
y = arange(10001)/10000.*a[size(a[:,0])-1,1];
x = x0 + L*pow(cos(al),2) - L*pow((cos(al) - 1./(2.0*L*sin(al))*(y - y0)), 2);
plot(x,y, 'k', linewidth=2); grid('on');
plot(a[:,0], a[:,1], linewidth=2); hold('on');

a = loadtxt('ray2.txt'); #plot(a[:,0], a[:,1]); hold('on');

al = 10*pi/180;
x0 = 0; y0 = 0;
y = arange(10001)/10000.*a[size(a[:,0])-1,1];
x = x0 + L*pow(cos(al),2) - L*pow((cos(al) - 1./(2.0*L*sin(al))*(y - y0)), 2);
plot(x,y, 'k', linewidth=2); grid('on');
plot(a[:,0], a[:,1], linewidth=2); hold('on');

a = loadtxt('ray3.txt'); #plot(a[:,0], a[:,1]); hold('on');

al = 30*pi/180;
x0 = 0; y0 = 0;
y = arange(10001)/10000.*a[size(a[:,0])-1,1];
x = x0 + L*pow(cos(al),2) - L*pow((cos(al) - 1./(2.0*L*sin(al))*(y - y0)), 2);
plot(x,y, 'k', linewidth=2); grid('on');
plot(a[:,0], a[:,1], linewidth=2); hold('on');

a = loadtxt('ray4.txt'); #plot(a[:,0], a[:,1]); hold('on');

al = 60*pi/180;
x0 = 0; y0 = 0;
y = arange(10001)/10000.*a[size(a[:,0])-1,1];
x = x0 + L*pow(cos(al),2) - L*pow((cos(al) - 1./(2.0*L*sin(al))*(y - y0)), 2);
plot(x,y, 'k', linewidth=2); grid('on');
plot(a[:,0], a[:,1], 'brown',  linewidth=2); hold('on');
grid('on');

show()
