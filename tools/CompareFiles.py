import sys
import numpy as np
import argparse

'''
This script compares two output files

Date: January 2020
Author: Simon Grimm
'''
def main(t0, t1, out):
	file0 = 'Out_%s.dat' % t0
	file1 = 'Out_%s.dat' % t1

	nu0, k0 = np.loadtxt(file0, usecols=(0,1), unpack = True)
	nu1, k1 = np.loadtxt(file1, usecols=(0,1), unpack = True)

	check = 1

	for i in range(len(k0)):
		if(nu0[i] != nu1[i]):
			print("Error, nu0 and nu1 not the same")
			break
		d = k0[i] - k1[i]

		if(d > 1.0e-10):
			print('%.20g %.20g %.20g %.20g' % (nu0[i], d, k0[i], k1[i]))
			check = 0
			
		if(out == 1):
			print(nu0[i], d)

	if(check == 0):
		sys.exit(100)


if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	parser.add_argument('-f0', '--f0', type=str,
                    help='File name 0', default = '')
	parser.add_argument('-f1', '--f1', type=str,
                    help='File name 1', default = '')
	parser.add_argument('-out', '--out', type=int,
                    help='print full output', default = 0)

	args = parser.parse_args()

	if(args.f0 == ''):
		print("Error, no f0 specified")
	if(args.f1 == ''):
		print("Error, no f1 specified")
		print("Error, no f1 specified")


	main(args.f0, args.f1, args.out)

