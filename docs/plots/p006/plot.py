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

dirname = os.path.abspath(os.path.dirname(__file__))
filename = os.path.join(dirname, 'Out_i.dat')
filenameb0 = os.path.join(dirname, 'Out_i_tr0000.dat')
filenameb1 = os.path.join(dirname, 'Out_i_tr0001.dat')
filenameb2 = os.path.join(dirname, 'Out_i_tr0002.dat')
filenameb3 = os.path.join(dirname, 'Out_i_tr0003.dat')
filenameb4 = os.path.join(dirname, 'Out_i_tr0004.dat')

ax1=pl.subplot(211)
pl.subplots_adjust(wspace=0.3, hspace=0.3)

nu, k = np.loadtxt(filename, unpack=True)
y0, k0 = np.loadtxt(filenameb0, unpack=True)
y1, k1 = np.loadtxt(filenameb1, unpack=True)
y2, k2 = np.loadtxt(filenameb2, unpack=True)
y3, k3 = np.loadtxt(filenameb3, unpack=True)
y4, k4 = np.loadtxt(filenameb4, unpack=True)


pl.plot(nu, k, lw = 0.5, label= r'full Opacity function $\kappa(\nu)$')

pl.xlabel(r'$\nu$ [cm$^{-1}$]')
pl.ylabel(r'$\kappa$ [cm$^2$ / g]')

pl.vlines(0, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(6000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(12000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(18000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(24000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(30000, 1E-10, 1E6, linestyles = ':', lw = 0.8)

pl.legend()

ax1.set_yscale('log')

pl.xlim([-500,30500])

ax1=pl.subplot(212)

pl.plot(y0,  k0, lw = 2, label =r'bin 0 $\tau$')
pl.plot((y1) * 1e23, k1, lw = 2, label =r'bin 1 $\tau$')
pl.plot((y2) * 1e46, k2, lw = 2, label =r'bin 2 $\tau$')
pl.plot((y3) * 1e69, k3, lw = 2, label =r'bin 3 $\tau$')
pl.plot((y4) * 1e92, k4, lw = 2, label =r'bin 4 $\tau$')

pl.legend(loc='lower right')


ax1.set_xscale('log')


ax1.set_xticks([1e-11, 1.0, 1e11, 3e22, 1e34, 3e45, 1e57, 3e68, 1e80, 3e91, 1e103])
ax1.set_xticklabels([1e-12, 1, 1e-12, 1, 1e-12, 1, 1e-12, 1, 1e-12, 1, '1e12'])

pl.xlabel('m[g cm−2]')
pl.ylabel(r'$\tau$')

pl.xlim([1e-13,1e105])


name = 'plot001.png'

pl.savefig(name, dpi=300)
pl.clf()
	

