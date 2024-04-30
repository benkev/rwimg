x = 1.0;
y = empty(26);
y[0] = .000001;
for i in range(25):
	y[i+1] = (2 - y[i]/x)*y[i]
print y
plot(y, color='k', linewidth=2); hold(True);
title(r'Solutions to $y_{i+1} = (2 - y_{i}/x_i) y_i$')
xlabel(r'$i$')
ylabel(r'$y_i$') 
