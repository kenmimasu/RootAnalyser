from matplotlib import pyplot as plt
import numpy as np

with open('Results/subset_metdist_test_DoubleHist.dat', 'r') as dat:
    xx, yy, zz, dz = np.loadtxt(dat, unpack=True)
    
nx, ny = len(set(xx)), len(set(yy))

XX = xx.reshape(nx,ny)
YY = yy.reshape(nx,ny)
ZZ =  np.flipud(zz.reshape(nx,ny))

fig, ax = plt.subplots()

xmin, xmax, ymin, ymax = min(xx), max(xx), min(yy), max(yy)

ax.imshow(ZZ, interpolation='none', extent = (min(xx), max(xx), min(yy), max(yy)) ,aspect=(xmax-xmin)/(ymax-ymin))

plt.show()