import struct

filename = "/home/sigrimm/scratch/EXOMOL/1H2-16O__POKAZATEL__00000-00100.bin"

file = open(filename, "rb")

for i in range(10):
	nu = struct.unpack('d', file.read(8))[0]
	S = struct.unpack('d', file.read(8))[0]
	EL = struct.unpack('d', file.read(8))[0]
	A = struct.unpack('d', file.read(8))[0]

	print(nu, S, EL, A)
