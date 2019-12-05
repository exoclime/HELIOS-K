#https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19700011342.pdf


'''
This script calculates the partition function from the Energy levels.

Date: November 2019
Author: Simon Grimm

'''


import numpy as np
import sys
import argparse

def gTot():
	#This function calculates the total statistical weight for a quantum number n
	for n in range(0, 15):
		gtot = 0
		for l in range(0, n):
			g = 2 * l + 1
			gtot += g * 2
			#print(n, l, g)
		print(n, gtot)




def partition(z, I):

	file = ("NIST_ELevels%02d%02d.dat" % (z, I))


	outfile = f = open('NIST%02d%02d.pf' % (z, I),'w')

	Conf, Term, J, g, E = np.loadtxt(file, dtype='str', skiprows=1, usecols=(0, 1, 2, 3,4), delimiter='\t',unpack=True)
	print('NIST%02d%02d.pf' % (z, I), len(Term))

	#check if the Term is given, Note that for H I Energy levels are duplicated
	#Duplicated levels can be filtered out by this check
	noTerm = 0

	for i in range(len(Term)):
		#print(Conf[i], Term[i], J[i], g[i], E[i])
		if Term[i] == '':
			noTerm = 1		
	if(noTerm == 1):
		print("no Term")


	for t in range(10, 10000, 10):
		T = float(t)

		Z = 0.0

		for i in range(len(g)):

			try:
				if(Term[i] != '' or z > 1):
					gf = float(g[i])
					Ef = float(E[i])
					#exp(-h c E /(k_B T))
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

	#gTot()
	partition(Z,I)


