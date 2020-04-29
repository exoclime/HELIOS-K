#https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19700011342.pdf

'''
This script generates the HELIOS-K binary files and <species>.param files 
for NIST line lists

Date: November 2019
Author: Simon Grimm

'''



import numpy as np
import struct
import math
import csv
import pandas
import argparse


def Lines2(Z, I, printA):

	datafile = ("NIST_Lines%02d%02d.dat" % (Z, I))

	def_c = 2.99792458e10
	def_NA =  6.0221412927e23

	masses = np.zeros(120)

	with open("masses.txt") as mf:
		lines = mf.readlines()

		j = 0
		z = ""
		for line in lines:
			j += 1
			i = 0
			w1 = ""
			w2 = ""
			w3 = ""
			w4 = ""
			w5 = ""
			for word in line.split():
				i += 1
				if(i == 1):
					w1 = word
				if(i == 2):
					w2 = word
				if(i == 4):
					w4 = word
				if(i == 5):
					w5 = word
			if(w1 == "Atomic" and w2 == "Number"):		
			
				z = int(w4)
				#print(w1, w2, w4, z)
			if(w1 == "Standard" and w2 == "Atomic" and w5 != ""):		

				w6 = w5.replace('[', '')
				w7 = w6.replace(']', '')
				w6 = w7.replace('(', '')
				w7 = w6.replace(')', '')
				w8 = w7.split(",")	
				#print(w1, w2, w4, z, w5, w8, w8[0])
				try:
					masses[z] = float(w8[0])
				except:
					masses[z] = 0.0


	#for i in range(len(masses)):
	#	print(i, masses[i])

	try:
		data = pandas.read_csv(datafile,quotechar='"',quoting=0)
		wn = data['wncm-1']
		A = data['Akis^-1']
		gUP = data['g_k']
		ELow = data['Eicm-1']
	except:
		A =[]
		wn=[]
		gUP=[]
		ELow=[]

	name = "NIST%.2d%.2d" % (Z, I)
	outname = "%s.bin" % name
	paramname = "%s.param" % name
	output_file = open(outname,"wb")

	print(paramname)
	nl = 0
	numax = 0

	if(printA == 1):
		Afile = open("NIST_A%02d%02d.dat" % (Z, I), "w")

	for i in range(len(A)):

		a = A[i].replace('"', '')
		a = a.replace('a', '')
		a = a.replace('b', '')
		if(a == ""):
			continue
		e = ELow[i].replace('"', '')
		e = e.replace('+x', '')
		e = e.replace('&dagger;', '')
		if(e == ""):
			continue
		w = wn[i].replace('"', '')
		if(w == ""):
			continue
		#print(a)
		Af = float(a)
		gUPf = float(gUP[i])
		ELowf = float(e)
		wnf = float(w)

		if(printA == 1):
			print(i, wnf, Af, ELowf, gUPf, Z, masses[Z], file = Afile)

		S = gUPf * Af /(8.0 * math.pi * def_c * wnf * wnf * masses[Z] / def_NA);

		s = struct.pack('d', wnf)
		output_file.write(s)
		s = struct.pack('d', S)
		output_file.write(s)
		s = struct.pack('d', ELowf)
		output_file.write(s)
		s = struct.pack('d', 0.0)
		output_file.write(s)
		s = struct.pack('d', 0.0)
		output_file.write(s)

		nl += 1
		numax = max(numax, wnf)


	output_file.close()
	print("number of lines", nl)

	if(nl > 0):
		with open(paramname,"w") as f:

			print("Database = 31", file = f)
			print("Molecule number = %.2d%.2d" % (Z, I), file = f)
			print("Name = %s" % name, file = f)
			print("Number of Isotopes = 1", file = f)
			print("#Id Abundance      Q(296K)   g     Molar Mass(g)  partition file :", file = f)
			print("0 1.0             0.0       0      %s        %s.pf" % (masses[Z], name), file = f)
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
			print("Version = ", file = f)

		f.close()

if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	parser.add_argument('-Z', '--Z', type=int,
		help='Z', default = 1)
	parser.add_argument('-I', '--I', type=int,
		help='I', default = 0)
	parser.add_argument('-printA', '--printA', type=int,
		help='print A to file', default = 0)
	
	args = parser.parse_args()

	Z = args.Z
	I = args.I
	printA = args.printA

	Lines2(Z,I, printA)

