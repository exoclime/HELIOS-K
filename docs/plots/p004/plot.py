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


filenameb0 = 'Out_i_bin.dat'
filenameb1 = 'Out_r_bin.dat'

fig, ax1=pl.subplots()


y0, k0 = np.loadtxt(filenameb0, unpack=True)
y1, k1 = np.loadtxt(filenameb1, unpack=True)



pl.xlabel(r'y')
pl.ylabel(r'$\kappa$ [cm$^2$ / g]')


pl.plot(y0, k0, lw = 2, label ='original')
pl.plot(y1, k1, lw = 2, label ='resampled')

ax1.set_yscale('log')
pl.ylim([1E-8,1E6])
pl.legend(loc='lower right')

ax2 = pl.axes([0.0,0.0,1.0,1.0])
pl.xlim([0.21,0.23])
pl.ylim([2e-6,5e-6])

ip = InsetPosition(ax1, [0.2,0.4,0.5, 0.5])
ax2.set_axes_locator(ip)
ax2.set_yscale('log')


mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

pl.plot(y0, k0, lw = 2, label ='i')
pl.plot(y1, k1, lw = 2, label ='r')

pl.xticks()
pl.yticks()




name = 'plot001.png'

pl.savefig(name, dpi=300)
pl.clf()
	

