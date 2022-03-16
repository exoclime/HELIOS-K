import numpy as np
import math
import struct
import os
import argparse

# This script calculates the molecular line intensities
# It writes data files to prepare the HELIOS-K binary files


# Date: January 2022
# Author: Simon Grimm


#choose if the file contains wavenumber or not
Wavenumber = 2			#1: Wavenumber, 2: vacuum wavelength, 3: air based wavelenght


def main(name, mass, pfname, printA):
	filename= '%s.asc' % name

	outname = "%s.bin" % name
	
	print(name, outname, mass)


	output_file = open(outname,"wb")

	numax = 0.0
	nl = 0	

	if(printA == 1):
		Afile = open("%s_A.dat" % name, "w")

	with open(filename) as f:
		line = f.readlines()


		for ii in range(len(line)):
			l = line[ii]

			#E in cm^-1
			#atomic line list format
			if(Wavenumber == 1):
				wn = float(l[0:10])
				wl = 1.0E7/wn		#wavelenght in nm
			if(Wavenumber == 2):
				wl = float(l[0:10])	#wavelenght in nm
				wn = 1.0E7/wl
			if(Wavenumber == 3):
				wlAir = float(l[0:10])	#wavelenght in nm


			loggf = float(l[10:17])
			JLow = float(l[17:22])
			ELow = float(l[22:32])
			JUP = float(l[32:37])
			EUP = float(l[37:48])
			code = l[48:70]
			

			e = 4.80320425E-10      #electron charge in cgs units [statcoulomb = cm^(3/2) g^(1/2) s^-1]
			c = 2.99792458E10       #Speed of light cm/s
			me = 9.1093835611E-28   #mass of electron in g
			NA = 6.0221412927e23	#Avogadro Constant  1/mol


			ELow = abs(ELow)
			EUP = abs(EUP)


			#somethimes ELOW is larger than EUP
			
			if(ELow > EUP):
				t = EUP
				EUP = ELow
				ELow = t

				t = JUP
				JUP = JLow
				JLow = t
				
			#convert air wavelength to vacuum wavelength
			#http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion

			if(Wavenumber == 3):
				if(wlAir > 200):
					wlAir = wlAir * 10  #convert nm to Angstrom
					s = 10000.0 / wlAir
					n = 1.0 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s*s) + 0.0001599740894897 / (38.92568793293 - s*s)
					wl = wlAir * n
					wl = wl * 0.1  #convert Angstrom to nm
				else:
					wl = wlAir
				wn = 1.0E7/wl		#wavelenght in nm
				

			gUP = 2 * JUP + 1
			gLow = 2 * JLow + 1

			A = 8.0 * math.pi * wn * wn * (10.0**loggf) / gUP * math.pi * e * e / (me * c)

			nl = nl + 1

			numax = max(numax, wn)

			S = math.pi * e * e * 10.0**loggf * NA / (c * c * me * mass) 

			A = 8.0 * math.pi * wn * wn * 10.0**loggf / gUP * math.pi * e * e / (me * c)
		
			if(printA == 1):
				print(i, wn, A, ELow, gUP, mass, file = Afile)	
			#print(wn, 1.0E7/wn, loggf, ELow, EUP, JLow, JUP, mass, A)

			#print(wn, S, A, ELow, EUP, gLow, gUP)

			s = struct.pack('d', wn)
			output_file.write(s)
			s = struct.pack('d', S)
			output_file.write(s)
			s = struct.pack('d', ELow)
			output_file.write(s)
			s = struct.pack('d', A)
			output_file.write(s)



	print(" Lines:",nl, end='')

	output_file.close()
	if(printA == 1):
		Afile.close()


	if(nl > 0):
		f = open("%s.param" % name,'w')

		print("Database = 2", file = f)
		print("Molecule number = 1", file = f)
		print("Name = %s" % name, file = f)
		print("Number of Isotopologues = 1", file = f)
		print("#Id Abundance      Q(296K)   g     Molar Mass(g)  partition file :", file = f)
		print("0 1.0             0.0       0      %s        %s" % (mass, pfname), file = f)
		print("Number of columns in partition File = 2", file = f)
		print("Number of line/transition files = 1", file = f)
		print("Number of lines per file :", file = f)
		print("%d" % nl, file = f)
		print("Line file limits :", file = f)
		print("0", file = f)
		print("%d" % (int(numax)+1), file = f)
		print("#ExoMol :", file = f)
		print("Number of states = 0", file = f)
		print("Number of columns in transition files = 0", file = f)
		print("Default value of Lorentzian half-width for all lines = 0.07", file = f)
		print("Default value of temperature exponent for all lines = 0.5", file = f)
		print("Version = %s" % filename, file = f)

		f.close()



if __name__ == '__main__':

	parser = argparse.ArgumentParser()

	parser.add_argument('-name', '--name', type=str,
		#help='name', default = 'coax')
		help='name', default = 'tiototo')

	parser.add_argument('-mass', '--mass', type=float,
		#help='Isotopologue mass g/mol', default = 28.0101)
		help='Isotopologue mass g/mol', default = 61.947545)

	parser.add_argument('-pfFile', '--pfFile', type=str,
		#help='partition function file name', default = '12C-16O__Li2015.pf')
		help='partition function file name', default = '46Ti-16O__Toto.pf')

	parser.add_argument('-printA', '--printA', type=int,
		help='print A to file', default = 0)


	args = parser.parse_args()
	name = args.name
	mass = args.mass
	pfname = args.pfFile
	printA = args.printA

	print("name: %s, mass: %g; pfFile: %s; printA: %d" % (name, mass, pfname, printA))

	main(name, mass, pfname, printA)

