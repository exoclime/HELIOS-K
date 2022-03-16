import struct
import argparse


#filename = "/home/sigrimm/scratch/EXOMOL/1H2-16O__POKAZATEL__00000-00100.bin"


def main(filename):
	file = open(filename, "rb")

	for i in range(100000000000):
		try:
			nu = struct.unpack('d', file.read(8))[0]
			S = struct.unpack('d', file.read(8))[0]
			EL = struct.unpack('d', file.read(8))[0]
			A = struct.unpack('d', file.read(8))[0]

			print(nu, S, EL, A)
		except:
			return


if __name__ == '__main__':

	parser = argparse.ArgumentParser()

	parser.add_argument('-name', '--name', type=str,
		help='name', default = '')

	args = parser.parse_args()
	name = args.name

	main(name)
