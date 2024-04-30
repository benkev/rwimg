import numpy as np

mu, sigma = 0, 0.1 # mean and standard deviation
s = np.random.normal(mu, sigma/2.0, 100000)*2


import matplotlib.pyplot as plt
count, bins, ignored = plt.hist(s, 50, normed=True)
plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
         np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
         linewidth=2, color='r')
plt.show()
