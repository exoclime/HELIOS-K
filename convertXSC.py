import numpy as np


def main(infile, outfile, unitScale):

	#filename = 'O3_293.0K-0.0Torr_28901.0-40999.0_118.xsc'

	outf = open(outfile,'w')


	with open(infile) as f:
		lines = f.readlines()
		header = lines[0]
		
		numin = float(header[20:30])
		numax = float(header[30:40])
		Np = int(header[40:47])
		T = float(header[47:54])
		p = float(header[54:60])
		
		dnu = (numax - numin) / (Np - 1)


		#print(numin, numax, Np, T, p, dnu)

		count = 0

		for i in range(1,len(lines)):
		#for i in range(1,5):
			line = lines[i]
			words = line.split()
			
			for w in words:
				nu = numin + dnu * count

				w = float(w) / unitScale

				print(nu, w)
				print(nu, w, file = outf)
				count = count + 1

	outf.close()

if __name__ == '__main__':


	#open param.dat file and extract name, Species name and unit
	with open('param.dat', 'r') as pfile:

		lines = pfile.readlines()
		for line in lines:
			if line.find("name") != -1:
				name = line.split("=")[1]
				#remove blank spaces
				name = name.strip()
				print("name =", name)
			if line.find("Species Name") != -1:
				M = line.split("=")[1]
				#remove blank spaces
				M = M.strip()
				print("Species Name =", M)

			if line.find("Units") != -1:
				unit = line.split("=")[1]
				#remove blank spaces
				unit = int(unit.strip())
				print("unit =", unit)
	#open species param file and extract mass
	with open('%s.param' % M, 'r') as pfile:

		lines = pfile.readlines()
		print(lines[5])
		mass = float(lines[5].split()[4])
		print("mass= ",mass)


	infile = ("%s.xsc" % name)
	outfile = ("Out_%s.dat" % name)
	
	print(infile, outfile)

	unitScale = 1.0
	def_NA = 6.0221412927e23  #Avogadro Constant  1/mol
	if(unit == 0):
		unitScale = 1.0 / def_NA * mass;


	main(infile, outfile, unitScale)


