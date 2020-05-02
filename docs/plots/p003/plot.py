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


filename = 'Out_i.dat'
filenameb0 = 'Out_i_bin0000.dat'
filenameb1 = 'Out_i_bin0001.dat'
filenameb2 = 'Out_i_bin0002.dat'
filenameb3 = 'Out_i_bin0003.dat'
filenameb4 = 'Out_i_bin0004.dat'

filenamee0 = 'Out_e_bin0000.dat'
filenamee1 = 'Out_e_bin0001.dat'
filenamee2 = 'Out_e_bin0002.dat'
filenamee3 = 'Out_e_bin0003.dat'
filenamee4 = 'Out_e_bin0004.dat'

ax1=pl.subplot(111)

nu, k = np.loadtxt(filename, unpack=True)
y0, k0 = np.loadtxt(filenameb0, unpack=True)
y1, k1 = np.loadtxt(filenameb1, unpack=True)
y2, k2 = np.loadtxt(filenameb2, unpack=True)
y3, k3 = np.loadtxt(filenameb3, unpack=True)
y4, k4 = np.loadtxt(filenameb4, unpack=True)

ye0, ke0 = np.loadtxt(filenamee0, unpack=True)
ye1, ke1 = np.loadtxt(filenamee1, unpack=True)
ye2, ke2 = np.loadtxt(filenamee2, unpack=True)
ye3, ke3 = np.loadtxt(filenamee3, unpack=True)
ye4, ke4 = np.loadtxt(filenamee4, unpack=True)


pl.plot(nu, k, lw = 0.5, label= r'full Opacity function $\kappa(\nu)$')

pl.xlabel(r'$\nu$ [cm$^{-1}$]')
pl.ylabel(r'$\kappa$ [cm$^2$ / g]')

pl.vlines(0, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(6000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(12000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(18000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(24000, 1E-10, 1E6, linestyles = ':', lw = 0.8)
pl.vlines(30000, 1E-10, 1E6, linestyles = ':', lw = 0.8)


pl.plot(y0*30000/5.0, k0, lw = 0.8, label ='bin 0 $\kappa(y)$')
pl.plot(y1*30000/5.0 + 1 * 30000 / 5.0, k1, lw = 0.8, label ='bin 1 $\kappa(y)$')
pl.plot(y2*30000/5.0 + 2 * 30000 / 5.0, k2, lw = 0.8, label ='bin 2 $\kappa(y)$')
pl.plot(y3*30000/5.0 + 3 * 30000 / 5.0, k3, lw = 0.8, label ='bin 3 $\kappa(y)$')
pl.plot(y4*30000/5.0 + 4 * 30000 / 5.0, k4, lw = 0.8, label ='bin 4 $\kappa(y)$')

pl.scatter(ye0*30000/5.0 + 0 * 30000 / 5.0, ke0, color='orange')
pl.scatter(ye1*30000/5.0 + 1 * 30000 / 5.0, ke1, color='green')
pl.scatter(ye2*30000/5.0 + 2 * 30000 / 5.0, ke2, color='red')
pl.scatter(ye3*30000/5.0 + 3 * 30000 / 5.0, ke3, color='magenta')
pl.scatter(ye4*30000/5.0 + 4 * 30000 / 5.0, ke4, color='brown')

ax1.set_yscale('log')
pl.ylim([1E-8,1E6])


pl.legend()


ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())

ax2.set_xticks([0, 3000, 6000, 9000, 12000, 15000, 18000, 21000, 24000, 27000, 30000])
ax2.set_xticklabels([0, 0.5, 0, 0.5, 0, 0.5, 0, 0.5, 0, 0.5, 1.0])
pl.xlabel('y')


name = 'plot001.png'

pl.savefig(name, dpi=300)
pl.clf()
	

