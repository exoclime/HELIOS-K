import numpy as np
import argparse
import matplotlib.pylab as plt
from numpy.polynomial.chebyshev import chebval, chebfit

'''
This script reconstructs sampled cbin files, to check the sampling

Date: January 2020
Author: Simon Grimm
'''




#change here the bin size:
fig = plt.figure()

def main(M):
	data = np.loadtxt('%s.dat' % M)
	data_c = np.loadtxt('%s_cbin.dat' % M)

	print('%s_cbin.dat' % M)
	print('%s_r.dat' % M)

	#data = np.loadtxt('1H2-16O__POKAZATEL_e0DAT/Out_00000_42000_01500_n800.dat')
	#data_c = np.loadtxt('1H2-16O__POKAZATEL_cbin/Out_01500_n800_cbin.dat')
	#data = np.loadtxt('12C2-1H2__HITRAN2016_e0DAT/Out_00000_10000_00100_n800.dat')
	#data_c = np.loadtxt('12C2-1H2__HITRAN2016_cbin/Out_00100_n800_cbin.dat')
	#data = np.loadtxt('14N-1H__Yueqi_e0DAT/Out_00000_17000_01500_n800.dat')
	#data_c = np.loadtxt('14N-1H__Yueqi_cbin/Out_01500_n800_cbin.dat')

	kStart = data[0,0]	
	Nbins = len(data_c)
	binSize0 = int(len(data) / Nbins)
	binSize = int(binSize0 * 1.0)	#here the Chebyshev output can be reduced by a factor e.g. 0.1

	print("number of nu points:",len(data), Nbins, binSize, kStart) 

	for binIndex in range(0, Nbins):

		print(binIndex)
		k = data[(binIndex * binSize0):(binIndex + 1) * binSize0, 1]
		x0 = data[(binIndex * binSize0):(binIndex + 1) * binSize0, 0]
		ks = np.sort(k)


		x = np.linspace(0, 1.0, num=binSize, endpoint=True)

		#extract Chebyshev coefficients
		c = data_c[binIndex,2:]
		#extract starting point in x of opacity function	
		xs = data_c[binIndex,1]

		#print(x)
		print(c)

		print(c[0], xs)

		#rescale x to the standard Chebychev polynomial range [-1:1]
		x1 = x * 2.0 - 1.0
		k_res = chebval(x1,c,tensor=False)
		x2 = x * (1.0 - xs) + xs

		k_res = np.exp(k_res)


		plt.plot(x0, k, c="blue", label='Full')
		plt.plot(x0, ks, c="red", label='Sorted')
		plt.plot(x2*(binSize0 - 1) + binIndex * binSize0 + kStart, k_res, c="green", label='Cheb Poly')

		if(binIndex == 0):
			plt.legend()

			plt.xlabel(r'$\nu$')
			plt.ylabel(r'k($\nu$)')

		plt.yscale('log')
		#plt.xscale('log')
	plt.show()

if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	parser.add_argument('-M', '--M', type=str,
		help='Species', default = '')

	args = parser.parse_args()

	M = args.M

	if(M == ''):
		print("Error, no species specified, run python exomol.py -M <id>")

	main(M)
