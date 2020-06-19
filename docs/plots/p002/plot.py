import os

import matplotlib
#matplotlib.use('PS')
matplotlib.use('agg')

import pylab as pl
import numpy as np

pl.rc('font', size=12)
params = {'legend.fontsize': 12}
pl.rcParams.update(params)


#yt = np.arange(1.0)*0.1

pl.figure(figsize=(8, 6))


dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, 'Out_i.dat')
filenameb0 = os.path.join(dirname, 'Out_i_bin0000.dat')
filenameb1 = os.path.join(dirname, 'Out_i_bin0001.dat')
filenameb2 = os.path.join(dirname, 'Out_i_bin0002.dat')
filenameb3 = os.path.join(dirname, 'Out_i_bin0003.dat')
filenameb4 = os.path.join(dirname, 'Out_i_bin0004.dat')


ax1=pl.subplot(111)

nu, k = np.loadtxt(filename, unpack=True)
y0, k0 = np.loadtxt(filenameb0, unpack=True)
y1, k1 = np.loadtxt(filenameb1, unpack=True)
y2, k2 = np.loadtxt(filenameb2, unpack=True)
y3, k3 = np.loadtxt(filenameb3, unpack=True)
y4, k4 = np.loadtxt(filenameb4, unpack=True)


pl.plot(nu, k, lw = 0.5, label= r'full Opacity function $\kappa(\nu)$')

pl.xlabel(r'$\nu$ [cm$^{-1}$]')
pl.ylabel(r'$\kappa$ [cm$^2$ / g]')

pl.vlines(2000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(6000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(12000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(15000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(25000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(30000, 1E-10, 1E6, linestyles = ':', lw = 0.8)


pl.plot(y0*4000 + 2000, k0, lw = 2, label ='bin 0 $\kappa(y)$')
pl.plot(y1*6000 + 6000, k1, lw = 2, label ='bin 1 $\kappa(y)$')
pl.plot(y2*3000 + 12000, k2, lw = 2, label ='bin 2 $\kappa(y)$')
pl.plot(y3*10000 + 15000, k3, lw = 2, label ='bin 3 $\kappa(y)$')
pl.plot(y4*5000 + 25000, k4, lw = 2, label ='bin 4 $\kappa(y)$')

ax1.set_yscale('log')
pl.ylim([1E-8,1E6])


pl.legend()


ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())

ax2.set_xticks([2000, 4000, 6000, 9000, 12000, 13500, 15000, 20000, 25000, 27500, 30000])
ax2.set_xticklabels([0, 0.5, 0, 0.5, 0, 0.5, 0, 0.5, 0, 0.5, 1.0])
pl.xlabel('y')


name = 'plot001.png'

pl.savefig(name, dpi=300)
pl.clf()
	

