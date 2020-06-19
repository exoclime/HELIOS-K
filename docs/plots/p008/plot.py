import os

import matplotlib
#matplotlib.use('PS')
matplotlib.use('agg')

import pylab as pl
import numpy as np
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)


pl.rc('font', size=12)
params = {'legend.fontsize': 12}
pl.rcParams.update(params)



pl.figure(figsize=(8, 6))

dirname = os.path.abspath(os.path.dirname(__file__))

filename1 = os.path.join(dirname, 'Out_c0_25.dat')
filename2 = os.path.join(dirname, 'Out_c0_500.dat')
filename3 = os.path.join(dirname, 'Out_c0_0.dat')

ax1=pl.subplot(111)

nu1, k1 = np.loadtxt(filename1, unpack=True)
nu2, k2 = np.loadtxt(filename2, unpack=True)
nu3, k3 = np.loadtxt(filename3, unpack=True)


pl.plot(nu1, k1, lw = 0.5, label= r'cut = 25 cm$^{-1}$')
pl.plot(nu2, k2, lw = 0.5, label= r'cut = 500 cm$^{-1}$')
pl.plot(nu3, k3, lw = 0.5, label= r'cut = $\infty$ cm$^{-1}$')


pl.xlabel(r'$\nu$ [cm$^{-1}$]')
pl.ylabel(r'$\kappa$ [cm$^2$ / g]')


ax1.set_yscale('log')
pl.ylim([1E-9,1E6])


pl.legend(loc='lower left')



ax2 = pl.axes([0.0,0.0,1.0,1.0])
pl.xlim([22500,23300])
pl.ylim([1e-9,1e-2])

ip = InsetPosition(ax1, [0.7,0.5,0.25, 0.45])
ax2.set_axes_locator(ip)
ax2.set_yscale('log')

pl.plot(nu1, k1, lw = 0.5, label= r'cut = 25 cm$^{-1}$')
pl.plot(nu2, k2, lw = 0.5, label= r'cut = 500 cm$^{-1}$')
pl.plot(nu3, k3, lw = 0.5, label= r'cut = $\infty$ cm$^{-1}$')

mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')
pl.xticks([])
pl.yticks([])


name = 'plot001.png'

pl.savefig(name, dpi=300)
pl.clf()
	

