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


#yt = np.arange(1.0)*0.1

pl.figure(figsize=(8, 6))

dirname = os.path.abspath(os.path.dirname(__file__))
filenameb0 = os.path.join(dirname, 'Out_i_bin.dat')
filenameb1 = os.path.join(dirname, 'Out_r_bin.dat')

fig, ax1=pl.subplots()


y0, k0 = np.loadtxt(filenameb0, unpack=True)
y1, k1 = np.loadtxt(filenameb1, unpack=True)



pl.xlabel(r'y')
pl.ylabel(r'$\kappa$ [cm$^2$ / g]')


pl.plot(y0, k0, lw = 2, label ='original')
pl.plot(y1[y1>0.8537], k1[y1>0.8537], lw = 1.5, label ='resampled')

pl.annotate('', xy=(1.0, 1e-6), xytext=(0.8537, 1e-6), 
            arrowprops=dict(facecolor='black', arrowstyle='|-|'),)
pl.text(0.91, 2e-6, 'Chebyshev support')

pl.annotate('', xy=(0.8537, 1e-6), xytext=(0.8, 1e-6), 
            arrowprops=dict(facecolor='black', arrowstyle='|-|'),)
pl.text(0.818, 2e-6, 'kmin')


ax1.set_yscale('log')
#pl.ylim([1E-8,1E6])
pl.xlim([0.8, 1.0])
pl.legend(loc='upper left')




name = 'plot001.png'

pl.savefig(name, dpi=300)
pl.clf()
	

