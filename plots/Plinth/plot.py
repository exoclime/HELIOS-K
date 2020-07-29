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
filename3 = os.path.join(dirname, 'plinth.dat')

ax1=pl.subplot(121)

nu1, k1 = np.loadtxt(filename1, unpack=True)
nu2, k2 = np.loadtxt(filename2, unpack=True)
il, nu3, k3 = np.loadtxt(filename3, unpack=True)


pl.plot(nu1, k1, lw = 1.5, label= r'1: full Voigt')
pl.plot(nu2, k2, lw = 1.0, label= r'2: Voight - plinth')
#pl.plot(nu3, 2.51241e-6, lw = 1.0, label= r'3: plinth')
ax1.add_patch(Rectangle((20318.99, 0), 1.0, 2.51241e-6,alpha=0.3, color='g'))
ax1.text(20319.4,1.0e-6, 'Plinth', color='g')

pl.xlabel(r'$\nu$ [cm$^{-1}$]')
pl.ylabel(r'$\kappa$ [cm$^2$ / g]')


#ax1.set_yscale('log')
pl.xlim([20318.8, 20320.2])
pl.ylim([0.0,0.00002])
pl.legend(loc='upper right')

ax1=pl.subplot(122)

nu1, k1 = np.loadtxt(filename1, unpack=True)
nu2, k2 = np.loadtxt(filename2, unpack=True)
il, nu3, k3 = np.loadtxt(filename3, unpack=True)


pl.plot(nu1, k1, lw = 1.5, label= r'1: full Voigt')
pl.plot(nu2, k2, lw = 1.0, label= r'2: Voight - plinth')
ax1.add_patch(Rectangle((20318.99, 3e-8), 1.0, 2.51241e-6-3e-8,alpha=0.3, color='g'))
ax1.text(20319.4,3.0e-7, 'Plinth', color='g')


pl.xlabel(r'$\nu$ [cm$^{-1}$]')


ax1.set_yscale('log')
pl.xlim([20318.8, 20320.2])
pl.ylim([3E-8,5E-5])

pl.legend(loc='upper right')

name = 'plot001.png'

pl.savefig(name, dpi=300)
pl.clf()
	

