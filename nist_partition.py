#https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19700011342.pdf


import numpy as np
import sys
import argparse

def partition(z, I):

	file = "test.dat"


	outfile = f = open('NIST%02d%02d.pf' % (z, I),'w')

	Term, g, E = np.loadtxt(file, dtype='str', skiprows=1, usecols=(1, 3,4), delimiter='\t',unpack=True)
	print('NIST%02d%02d.pf' % (z, I), len(Term))


	for t in range(10, 10000, 10):
		T = float(t)

		Z = 0.0

		for i in range(len(g)):

			try:
				gf = float(g[i])
				Ef = float(E[i])
				z = gf * np.exp(-1.4387770 * Ef/T)
				if(t == 0):
					print(i, gf, Ef, Z, z)

				Z += gf * np.exp(-1.4387770 * Ef/T)
			except:
				continue


		print(T, Z, file = outfile)

	outfile.close()

if __name__ == "__main__":


	parser = argparse.ArgumentParser()

	parser.add_argument('-Z', '--Z', type=int,
		help='Z', default = 1)
	parser.add_argument('-I', '--I', type=int,
		help='I', default = 0)
	
	args = parser.parse_args()

	Z = args.Z
	I = args.I

	partition(Z,I)


