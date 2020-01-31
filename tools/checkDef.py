import exomol 
import numpy as np


M = np.loadtxt("Exomol_species.dat", usecols=(2,), unpack = True, dtype='str')


#diff -x '*.txt' -x '*.bin' -x '*.param' -x '*.py' -x '*.trans' -x '*.states' -x '*.par' EXOMOL/ EXOMOL_test/ > test
for m in M:
	print("********", m, "**************")

	exomol.main(m, 2, 0) 
