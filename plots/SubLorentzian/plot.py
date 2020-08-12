import os

import matplotlib
#matplotlib.use('PS')
matplotlib.use('agg')

import pylab as pl
import numpy as np
from matplotlib.patches import Rectangle

pl.rc('font', size=12)
params = {'legend.fontsize': 12}
pl.rcParams.update(params)



pl.figure(figsize=(8, 6))

pl.subplots_adjust(left=0.17, bottom=None, right=None, top=None, wspace=None, hspace=None)

dirname = os.path.abspath(os.path.dirname(__file__))
print(dirname)
filename1 = os.path.join(dirname, 'Out_i.dat')
filename2 = os.path.join(dirname, 'Out_j.dat')
filename3 = os.path.join(dirname, 'Out_k.dat')

ax1=pl.subplot(211)

nu1, k1 = np.loadtxt(filename1, unpack=True)
nu2, k2 = np.loadtxt(filename2, unpack=True)


pl.plot(nu2, k2, lw = 1.0, color = 'orange', label= r'2: sub-Lorentzian wings, 2100-2600 cm$^{-1}$')
pl.plot(nu1, k1, lw = 1.5, color = 'blue', label= r'1: full Voigt')

pl.ylabel(r'$\kappa$ [cm$^2$ / g]')


ax1.set_yscale('log')
#pl.xlim([20318.8, 20320.2])
pl.ylim([1.0e-18,3.0e5])
pl.legend(loc='lower right')

ax1=pl.subplot(212)

nu1, k1 = np.loadtxt(filename1, unpack=True)
nu2, k2 = np.loadtxt(filename3, unpack=True)


pl.plot(nu2, k2, lw = 1.0, color = 'orange', label= r'2: sub-Lorentzian wings, full range')
pl.plot(nu1, k1, lw = 1.5, color = 'blue', label= r'1: full Voigt')

pl.xlabel(r'$\nu$ [cm$^{-1}$]')
pl.ylabel(r'$\kappa$ [cm$^2$ / g]')


ax1.set_yscale('log')
#pl.xlim([20318.8, 20320.2])
pl.ylim([1.0e-18,3.0e5])

pl.legend(loc='lower right')

name = 'plot001.png'

pl.savefig(name, dpi=300)
pl.clf()
	

