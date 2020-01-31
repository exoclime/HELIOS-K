import struct
  
filename = "Out_00000_09000_02900_p300.bin"

file = open(filename, "rb")

for i in range(20000):
        K = struct.unpack('f', file.read(4))[0]

        print(K)

