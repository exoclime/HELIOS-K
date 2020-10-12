#this script converts a HELIOS-K opacity *.dat file to a binary opacity *.bin file

#Date: October 2020
#Author: Simon Grimm

import numpy as np
import struct
import argparse


def main(filename):


	datFilename = '%s.dat' % filename
	binFilename = '%s.bin' % filename

	binFile = open(binFilename, "wb")

	nu, K = np.loadtxt(datFilename, unpack=True)

	for i in range(len(nu)):
		kb = struct.pack('f', K[i])
		binFile.write(kb)

		#print(nu[i], K[i])



	binFile.close()


if __name__ == '__main__':


	parser = argparse.ArgumentParser()

	parser.add_argument('-n', '--name', type=str,
		help='file name', default = 'Out_i')


	args = parser.parse_args()

	print("Convert %s.dat to %s.bin" % (args.name, args.name))

	main(args.name)

