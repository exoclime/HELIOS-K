import numpy as np
import math
import struct
import os
import argparse

# This script corrects missing isotope fractions in hyperfine level splittings


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




def main(Z, I):

	if(Z == -1 and I == -1):
		# all species in Z and I
		for i in range(0,100):
			for j in range(0,3):
				processLineList(i, j)

	if(Z == -1 and I > -1):
		# all species in Z
		for i in range(0,100):
			processLineList(i, I)

	if(Z > -1 and I == -1):
		# all species in I
		for j in range(0,3):
			processLineList(Z - 1, j)

	if(Z > -1 and I > -1):
		processLineList(Z - 1, I)
	

def processLineList(i, j):
	# i molecule id 0 to 100
	# j ion id 0 to 3

	name ="gfnew%02d%02d" % (i + 1, j)
	hyperfile = open("gfnewCorr%02d%02d.dat" % (i + 1, j), "w")

	inname = "%s.dat" % name
	el=elt0[i]


	numax = 0.0

	filesize = os.path.getsize(inname)
	if(filesize != 0):
		#ELOWAverage and EUPAverage do not include hyperfine energy shifts
		wn, S, A, ELow, EUP, ELowAverage, EUPAverage, gLow, gUP, GammaR, isotope, hyperFineFraction, ISOFraction, sameLabel = np.loadtxt(inname, unpack=True)
	else:
		exit()

	for i in range(len(wn)):


		if(sameLabel[i] == 0.0):
			HFISO = hyperFineFraction[i] * ISOFraction[i]
		else:
			HFISO += hyperFineFraction[i] * ISOFraction[i]

		#print("-", wn[i], S[i], A[i], ELow[i], EUP[i], gLow[i], gUP[i], GammaR[i], isotope[i], hyperFineFraction[i], ISOFraction[i], sameLabel[i], HFISO)
		
		if(sameLabel[i] == 1.0 and ISOFraction[i] == 1.0 and hyperFineFraction[i] == 1.0 and hyperFineFraction[i - 1] != 1.0):
			HFISOj = HFISO
			for j in range(1,20):

				if(sameLabel[i + j] == 0.0 or i + j >= len(wn)):
					break
	
				if(ISOFraction[i + j] == 1.0 and hyperFineFraction[i + j == 1.0]):
					print("Error a, two lines with IsoFraction 1.0", wn[i])

				HFISOj += hyperFineFraction[i + j] * ISOFraction[i + j]
			print("ISOa", HFISOj, 2.0 - HFISOj)
			ISOFraction[i] = 2.0 - HFISOj
			HFISO -= 1.0
			HFISO += hyperFineFraction[i] * ISOFraction[i]
				

		if(i + 1 < len(wn) and sameLabel[i] == 0.0 and sameLabel[i + 1] == 1.0 and ISOFraction[i] == 1.0 and hyperFineFraction[i] == 1.0 and hyperFineFraction[i + 1] != 1.0):
			HFISOj = HFISO
			for j in range(1,20):

				if(sameLabel[i + j] == 0.0 or i + j >= len(wn)):
					break
	
				if(ISOFraction[i + j] == 1.0 and hyperFineFraction[i + j] == 1.0):
					print("Error b, two lines with IsoFraction 1.0", wn[i])

				HFISOj += hyperFineFraction[i + j] * ISOFraction[i + j]
			print("ISOb", HFISOj, 2.0 - HFISOj)
			ISOFraction[i] = 2.0 - HFISOj
			HFISO -= 1.0
			HFISO += hyperFineFraction[i] * ISOFraction[i]

		#print("+", wn[i], S[i], A[i], ELow[i], EUP[i], gLow[i], gUP[i], GammaR[i], isotope[i], hyperFineFraction[i], ISOFraction[i], sameLabel[i], HFISO)
		
		

		if(HFISO >= 1.1):
			print("Warning HFISO > 1.001:", HFISO, wn[i])

		print(wn[i], S[i], A[i], ELow[i], EUP[i], ELowAverage[i], EUPAverage[i], gLow[i], gUP[i], GammaR[i], isotope[i], hyperFineFraction[i], ISOFraction[i], sameLabel[i], file = hyperfile)


	hyperfile.close()

if __name__ == '__main__':

	parser = argparse.ArgumentParser()

	parser.add_argument('-Z', '--Z', type=int,
		help='Z', default = -1)
	parser.add_argument('-I', '--I', type=int,
		help='I', default = -1)


	args = parser.parse_args()
	Z = args.Z
	I = args.I

	print("Z:%d, I: %d" % (Z, I))

	main(Z, I)

