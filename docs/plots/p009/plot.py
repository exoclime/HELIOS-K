import os

import matplotlib
#matplotlib.use('PS')
matplotlib.use('agg')

import pylab as pl
import numpy as np


pl.rc('font', size=12)
params = {'legend.fontsize': 12}
pl.rcParams.update(params)



pl.figure(figsize=(8, 6))

dirname = os.path.abspath(os.path.dirname(__file__))
print(dirname)
filename1 = os.path.join(dirname, 'Out_p1.dat')
filename2 = os.path.join(dirname, 'Out_p2.dat')
filename3 = os.path.join(dirname, 'Out_p3.dat')
filename4 = os.path.join(dirname, 'Out_p4.dat')

ax1=pl.subplot(111)

nu1, k1 = np.loadtxt(filename1, unpack=True)
nu2, k2 = np.loadtxt(filename2, unpack=True)
nu3, k3 = np.loadtxt(filename3, unpack=True)
nu4, k4 = np.loadtxt(filename4, unpack=True)


pl.plot(nu1, k1, lw = 1.5, label= r'1: Voigt')
pl.plot(nu2, k2, lw = 1.0, label= r'2: Lorentz')
pl.plot(nu3, k3, lw = 1.0, label= r'3: Doppler')
pl.plot(nu3, k4, lw = 0.5, label= r'4: Integrated Binned Gaussian')


pl.xlabel(r'$\nu$ [cm$^{-1}$]')
pl.ylabel(r'$\kappa$ [cm$^2$ / g]')


ax1.set_yscale('log')
pl.xlim([24455.9, 24458.4])
pl.ylim([1E-7,2E-4])


pl.legend(loc='upper right',ncol=2)

name = 'plot001.png'

pl.savefig(name, dpi=300)
pl.clf()
	

