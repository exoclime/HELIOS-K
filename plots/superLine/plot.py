import os

import matplotlib
#matplotlib.use('PS')
matplotlib.use('agg')

import pylab as pl
import numpy as np
from matplotlib.patches import Rectangle
matplotlib.rcParams['agg.path.chunksize'] = 10000

pl.rc('font', size=12)
params = {'legend.fontsize': 12}
pl.rcParams.update(params)



pl.figure(figsize=(8, 6))


dirname = os.path.abspath(os.path.dirname(__file__))
print(dirname)
filename1 = os.path.join(dirname, 'Out_all.dat')
filename2 = os.path.join(dirname, 'Out_super.dat')

ax1=pl.subplot(211)

nu1, k1 = np.loadtxt(filename1, unpack=True)
nu2, k2 = np.loadtxt(filename2, unpack=True)


pl.plot(nu1, k1, lw = 1.5, label= r'all lines')
pl.plot(nu2, k2, lw = 0.5, label= r'super lines')

#pl.xlabel(r'$\nu$ [cm$^{-1}$]')
pl.ylabel(r'$\kappa$ [cm$^2$ / g]')


ax1.set_yscale('log')
pl.xlim([0, 43000])
pl.legend(loc='upper right')

ax1=pl.subplot(212)

nu1, k1 = np.loadtxt(filename1, unpack=True)
nu2, k2 = np.loadtxt(filename2, unpack=True)

dk = (k1 - k2) / k1

print(np.max(dk))

pl.plot(nu2, dk, lw = 1.0, c='g',label= r'(all - super) / all')

pl.xlabel(r'$\nu$ [cm$^{-1}$]')
pl.ylabel(r'fractional difference')


ax1.set_yscale('log')
pl.xlim([0, 43000])
pl.ylim([0.001,1])
pl.legend(loc='upper right')

name = 'SuperLine.png'

pl.savefig(name, dpi=300)
pl.clf()
	

