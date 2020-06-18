import numpy as np
import math
import struct
import os
import argparse

# This script downloads the Kurucz gfnew files and calculates the line intensities
# It writes data files to prepare the HELIOS-K binary files


# Date: June 2020
# Author: Simon Grimm


elt0=[
[  100, "H"  , 1.00794],
[  200, "He" , 4.002602],
[  300, "Li" , 6.941],
[  400, "Be" , 9.012182],
[  500, "B"  , 10.811],
[  600, "C"  , 12.011],
[  700, "N"  , 14.00674],
[  800, "O"  , 15.9994],
[  900, "F"  , 18.9984032],
[ 1000, "Ne" , 20.1797],
[ 1100, "Na" , 22.989768],
[ 1200, "Mg" , 24.3050],
[ 1300, "Al" , 26.981539],
[ 1400, "Si" , 28.0855],
[ 1500, "P"  , 30.973762],
[ 1600, "S"  , 32.066],
[ 1700, "Cl" , 35.4527],
[ 1800, "Ar" , 39.948],
[ 1900, "K"  , 39.0983],
[ 2000, "Ca" , 40.078],
[ 2100, "Sc" , 44.955910],
[ 2200, "Ti" , 47.88],
[ 2300, "V"  , 50.9415],
[ 2400, "Cr" , 51.9961],
[ 2500, "Mn" , 54.93805],
[ 2600, "Fe" , 55.847],
[ 2700, "Co" , 58.93320],
[ 2800, "Ni" , 58.6934],
[ 2900, "Cu" , 63.546],
[ 3000, "Zn" , 65.39],
[ 3100, "Ga" , 69.723],
[ 3200, "Ge" , 72.61],
[ 3300, "As" , 74.92159],
[ 3400, "Se" , 78.96],
[ 3500, "Br" , 79.904],
[ 3600, "Kr" , 83.80],
[ 3700, "Rb" , 85.4678],
[ 3800, "Sr" , 87.62],
[ 3900, "Y"  , 88.90585],
[ 4000, "Zr" , 91.224],
[ 4100, "Nb" , 92.90638],
[ 4200, "Mo" , 95.94],
[ 4300, "Tc" , 97.9072],
[ 4400, "Ru" ,101.07],
[ 4500, "Rh" ,102.90550],
[ 4600, "Pd" ,106.42],
[ 4700, "Ag" ,107.8682],
[ 4800, "Cd" ,112.411],
[ 4900, "In" ,114.818],
[ 5000, "Sn" ,118.710],
[ 5100, "Sb" ,121.757],
[ 5200, "Te" ,127.60],
[ 5300, "I"  ,126.90447],
[ 5400, "Xe" ,131.29],
[ 5500, "Cs" ,132.90543],
[ 5600, "Ba" ,137.327],
[ 5700, "La" ,138.9055],
[ 5800, "Ce" ,140.115],
[ 5900, "Pr" ,140.90765],
[ 6000, "Nd" ,144.24],
[ 6100, "Pm" ,144.9127],
[ 6200, "Sm" ,150.36],
[ 6300, "Eu" ,151.965],
[ 6400, "Gd" ,157.25],
[ 6500, "Tb" ,158.92534],
[ 6600, "Dy" ,162.50],
[ 6700, "Ho" ,164.93032],
[ 6800, "Er" ,167.26],
[ 6900, "Tm" ,168.93421],
[ 7000, "Yb" ,173.04],
[ 7100, "Lu" ,174.967],
[ 7200, "Hf" ,178.49],
[ 7300, "Ta" ,180.9479],
[ 7400, "W"  ,183.84],
[ 7500, "Re" ,186.207],
[ 7600, "Os" ,190.23],
[ 7700, "Ir" ,192.22],
[ 7800, "Pt" ,195.08],
[ 7900, "Au" ,196.96654],
[ 8000, "Hg" ,200.59],
[ 8100, "Tl" ,204.3833],
[ 8200, "Pb" ,207.2],
[ 8300, "Bi" ,208.98037],
[ 8400, "Po" ,208.9824],
[ 8500, "At" ,209.9871],
[ 8600, "Rn" ,222.0176],
[ 8700, "Fr" ,223.0197],
[ 8800, "Ra" ,226.0254],
[ 8900, "Ac" ,227.0278],
[ 9000, "Th" ,232.0381],
[ 9100, "Pa" ,231.03588],
[ 9200, "U"  ,238.0289],
[ 9300, "Np" ,237.0482],
[ 9400, "Pu" ,244.0642],
[ 9500, "Am" ,243.0614],
[ 9600, "Cu" ,247.0703],
[ 9700, "Bk" ,247.0703],
[ 9800, "Cf" ,251.0796],
[ 9900, "Es" ,252.0830],
[10000, "Fm" ,257.0951],
[10100, "Md" ,258.0984],
[10200, "No" ,259.1011],
[10300, "Lr" ,262.1098],
[10400, "Rf" ,261.1089],
[10500, "Db" ,262.1144],
[10600, "Sg" ,263.1186],
[10700, "Bh" ,264.12],
[10800, "Hs" ,265.1306],
[10900, "Mt" ,268.00],
[11000, "Ds" ,268.00],
[11100, "Rg" ,272.00],
[11200, "Cn" ,277.00],
[11300, "Uut" ,0.00],
[11400, "Fl" ,289.00],
[11500, "Uup" ,0.00],
[11600, "Lv" ,289.00],
[11700, "Uus" ,294.00],
[11800, "Uuo" ,293.00]
]




def main(Z, I, printA):
	# i molecule id 0 to 100
	# j ion id 0 to 3

	el=elt0[Z - 1]

	if(I == 0):
		el[1] = el[1] + " 1"
	if(I == 1):
		el[1] = el[1] + " 2"
	if(I == 2):
		el[1] = el[1] + " 3"

	els= "% 6.2f" % (el[0] / 100.0)

	name ="vald%02d%02d" % (Z, I)
	NISTname ="NIST%02d%02d" % (Z, I)
	inname = "VALD%02d%02d.dat" % (Z, I)
	outname = "%s.bin" % name
	mass = el[2]

	exists = os.path.isfile(inname)
	numax = 0.0
	nl = 0	

	if(exists != 0):
		
		print(el[0], els, el[1], inname, outname, mass)


		output_file = open(outname,"wb")



		if(printA == 1):
			Afile = open("vald_A%02d%02d.dat" % (Z, I), "w")
	
		wn_a = []
		S_a = []
		ELow_a = []
		gLow_a = []
		EUp_a = []
		gUp_a = []
		GammaRad_a = []
		A_a = []


		with open(inname) as f:
			line = f.readlines()


			for ii in range(len(line)):
			#for ii in range(50):
				l = line[ii]

				l1 = l.find(el[1])
				if(l1 != 1):
					continue

				l2 = l.split(",");

				#print(l2)

				wn = float(l2[1])
				loggf = float(l2[2])
				ELow = float(l2[3])
				JLow = float(l2[4])
				EUp = float(l2[5])
				JUp = float(l2[6])
				GammaRad = float(l2[10])
			
				wl = 1.0E7/wn		#wavelenght in nm

				e = 4.80320425E-10      #electron charge in cgs units [statcoulomb = cm^(3/2) g^(1/2) s^-1]
				c = 2.99792458E10       #Speed of light cm/s
				me = 9.1093835611E-28   #mass of electron in g
				NA = 6.0221412927e23	#Avogadro Constant  1/mol


				gUp = 2 * JUp + 1
				gLow = 2 * JLow + 1
				
				A = 8.0 * math.pi * wn * wn * (10.0**loggf) / gUp * math.pi * e * e / (me * c)

			
				S = math.pi * e * e * 10.0**loggf * NA / (c * c * me * mass)


				wn_a.append(wn)
				S_a.append(S)
				ELow_a.append(ELow)
				gLow_a.append(gLow)
				EUp_a.append(EUp)
				gUp_a.append(gUp)
				GammaRad_a.append(GammaRad)
				A_a.append(A)

				nl = nl + 1

				numax = max(numax, wn)

		for i in range(len(wn_a)):
		
			wn = wn_a[i]
			S = S_a[i]
			ELow = ELow_a[i]
			EUp = EUp_a[i]
			gLow = gLow_a[i]
			gUp = gUp_a[i]
			GammaRad = GammaRad_a[i]
			A = A_a[i]

	
			#Compute GammaRad
			GammaR = 0.0
			#print(" ", wn, EUp_a[i], gUp_a[i], "|", ELow_a[i], gLow_a[i], A_a[i])
			for j in range(len(wn_a)):

				if(EUp_a[i] == EUp_a[j] and gUp_a[i] == gUp_a[j] and ELow_a[j] < EUp_a[i]):
					#print("-", EUp_a[j], gUp_a[j], "|", ELow_a[j], gLow_a[j], A_a[j])
					GammaR += A_a[j]

			for j in range(len(wn_a)):
				if(ELow_a[i] == EUp_a[j] and gLow_a[i] == gUp_a[j] and ELow_a[j] < ELow_a[i]):
					#print("+", EUp_a[j], gUp_a[j], "|", ELow_a[j], gLow_a[j], A_a[j])
					GammaR += A_a[j]
			#print(GammaR, 10**GammaRad)		

			if(printA == 1):
				print(ii, wn, S, A, ELow, EUp, gLow, gUp, GammaR, 10**GammaRad, file = Afile)
			#print(ii, wn, S, A, ELow, EUp, gLow, gUp, 10**GammaRad)

			if(i % 10000 == 0):
				print("reached line %d from %d" % (i, len(wn_a)))

			G = 10**GammaRad
			if(G == 1.0):
				G = GammaR


			s = struct.pack('d', wn)
			output_file.write(s)
			s = struct.pack('d', S)
			output_file.write(s)
			s = struct.pack('d', ELow)
			output_file.write(s)
			s = struct.pack('d', 0.0)
			output_file.write(s)
			s = struct.pack('d', G)
			output_file.write(s)



		print("                                                     Lines: ",nl)

		output_file.close()
		if(printA == 1):
			Afile.close()



	if(nl > 0):
		f = open("%s.param" % name,'w')

		print("Database = 32", file = f)
		print("Molecule number = %d" % el[0], file = f)
		print("Name = %s" % name, file = f)
		print("Number of Isotopes = 1", file = f)
		print("#Id Abundance      Q(296K)   g     Molar Mass(g)  partition file :", file = f)
		print("0 1.0             0.0       0      %s        %s.pf" % (mass, NISTname), file = f)
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
		print("Default value of Lorentzian half-width for all lines = 0.0", file = f)
		print("Default value of temperature exponent for all lines = 0.0", file = f)
		print("Version = 0",file = f)

		f.close()




if __name__ == '__main__':

	parser = argparse.ArgumentParser()

	parser.add_argument('-Z', '--Z', type=int,
		help='Z', default = -1)
	parser.add_argument('-I', '--I', type=int,
		help='I', default = -1)
	parser.add_argument('-printA', '--printA', type=int,
		help='print A to file', default = 0)


	args = parser.parse_args()
	Z = args.Z
	I = args.I
	printA = args.printA

	print("Z:%d, I: %d, printA: %d" % (Z, I, printA))

	main(Z, I, printA)

