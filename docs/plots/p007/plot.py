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
filenamem = os.path.join(dirname, 'Out_i_mean.dat')

ax1=pl.subplot(111)

nu, k = np.loadtxt(filename, unpack=True)

mean = np.loadtxt(filenamem, unpack=True)

print(mean)

pl.plot(nu, k, lw = 0.5, label= r'full Opacity function $\kappa(\nu)$')

pl.hlines(mean[0], 0, 30000, colors ='g', lw = 2.0, label= 'Planck Mean')
pl.hlines(mean[1], 0, 30000, colors = 'orange', lw = 2.0, label= 'Rosseland Mean')

pl.xlabel(r'$\nu$ [cm$^{-1}$]')
pl.ylabel(r'$\kappa$ [cm$^2$ / g]')


ax1.set_yscale('log')
pl.ylim([1E-8,1E6])


pl.legend()


name = 'plot001.png'

pl.savefig(name, dpi=300)
pl.clf()
	

