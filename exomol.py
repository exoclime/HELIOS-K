# This script downloads and unpacks  the *.states, *.pf, *.trans
# and *.def files from wwww.exomol.com.
# And it generates the < species name >.param files

# Date: May 2019
# Author: Simon Grimm

#run with "python3 exomol.py -M id" where id is the name of the molecule

import sys
import os
import subprocess
import numpy as np
import argparse
import math
import exomol2


def getTime(url):

	com = "wget -S --spider %s 2>&1 | grep -i Last-Modified" % url

	date = str(subprocess.check_output(com, shell=True))
	#cut string
	date = date.split('Last-Modified:')[1]
	date = date.split('\\n')[0]
	#convert to timestamp
	com = "date -d \'%s\' +%%s" % date
	date1 = str(subprocess.check_output(com, shell=True))
	date1 = date1.split('b\'')[1]
	date1 = date1.split('\\n')[0]

	print(date, date1)

	return int(date1)


# PrintISO: 1: create the <species>.param files, 0: no <species>.param files are created
# DownloadFiles: 0: no download
#		 1: download all files (*.def, *.pf, *.states and *.trans files)
#		 2: download only *.def and *.pf files
#		 3: download all file for super lines (*.def, *.pf, *.super)
def main(M, DownloadFiles, PrintISO, getTimeStamp, Temp):


	if(getTimeStamp == 1 and DownloadFiles == 0):
		print("Error, whe using getTimeStamp, then DownloadFiles must be 1 or 2")
		exit()

	print("Molecule = %s" % M)

	P = ""
	s = 0		#file range
	ntcol = 0	#columns in transition files
	npfcol = 0	#columns in partition function file

	dL = ""
	dn = ""
	dg = 0		#number of digits in .trans files ranges

	exfile = "Exomol_species.dat"
	date = 0
	date1 = 0


	#check if Exomol file exists and create it otherwise
	exists = os.path.isfile(exfile)
	if(exists == 0):
		exomol2.main()

	M0, P0, s0, nn0, dg0 = np.loadtxt(exfile, usecols=(2,3,4,5,6), unpack=True, dtype=str)

	if(M == "SH_X-X_A-X"):
		P = M
		s = 0
		nn = 1
		dg = 0
	else:
		#M = "1H2-16O__BT2"
		P = P0[M0 == M][0]
		s = int(s0[M0 == M][0])
		nn = int(nn0[M0 == M][0])
		dg = int(dg0[M0 == M][0])

	if(M == '75As-1H3__CYT18'):
		dg = 5


	print(M, P, s, nn, dg)

	if(DownloadFiles == 2):
		#download *.def and *.pf files, don't exit when files are missing
		com = "wget http://exomol.com/db/%s/%s.def" % (P, M)
		er=os.system(com)
		if(er != 0):
			print("Error in download .def file")

		com = "wget http://exomol.com/db/%s/%s.pf" % (P, M)
		er=os.system(com)
		if(er != 0):
			print("Error in download .pf file")

		#return 0

	if(DownloadFiles == 1):
		#downlaod *.def, *.pf and *.states files, exit when files are missing
		com = "wget http://exomol.com/db/%s/%s.def" % (P, M)
		er=os.system(com)
		if(er != 0):
			print("Error in download .def file")
			exit()
		com = "wget http://exomol.com/db/%s/%s.pf" % (P, M)
		er=os.system(com)
		if(er != 0):
			print("Error in download .pf file")
			exit()
	 
		print("er", er)
		com = "wget http://exomol.com/db/%s/%s.states.bz2" % (P, M)
		er=os.system(com)
		if(er == 0):
			com = "bzip2 -d %s.states.bz2" % M
			er=os.system(com)
		else:
			com = "wget http://exomol.com/db/%s/%s.states" % (P, M)
			er=os.system(com)
			if(er != 0):
				print("Error in download .states file")
				exit()
	if(DownloadFiles == 3):
		#downlaod *.def and *.pf files, exit when files are missing
		com = "wget http://exomol.com/db/%s/%s.def" % (P, M)
		er=os.system(com)
		if(er != 0):
			print("Error in download .def file")
			exit()
		com = "wget http://exomol.com/db/%s/%s.pf" % (P, M)
		er=os.system(com)
		if(er != 0):
			print("Error in download .pf file")
			exit()

	if(getTimeStamp == 1):
		url = "http://exomol.com/db/%s/%s.states.bz2" % (P, M)

		#first check if file exists
		com = "wget -S --spider %s" % url
		er=os.system(com)
		#print("er", er)
		if(er == 0):
			date = getTime(url)
			date1 = max(date, date1)
		else:
			url = "http://exomol.com/db/%s/%s.states" % (P, M)

			#first check if file exists
			com = "wget -S --spider %s" % url
			er=os.system(com)
			#print("er", er)
			if(er == 0):

				date = getTime(url)
				date1 = max(date, date1)

	#check if .pf file exists
	exists = os.path.isfile("%s.pf" % M)
	if(exists == 0):
		print("Error, partition file not found: %s.pf" % M)
		return 0
	#determine the number of columns in the partition file
	with open("%s.pf" % M) as pfFile:
		line = pfFile.readline()
		npfcol = len(line.split())
		print("number of columns in pf file", npfcol)

	
	#check if .def file exists
	exists = os.path.isfile("%s.def" % M)
	if(exists == 0):
		print("Error, definition file not found: %s.def" % M)

		if(M == "SH_X-X_A-X"):
			n = 1	
			mass = 32.9798960322
			dL = 0.07
			dn = 0.5
			version = 0
			nuMax = 39000.0 
		else:
			return 0
	else:
		n = 1
		nuMax = 0.0
		with open("%s.def" % M) as defFile:
			for line in defFile:
				if not "No. of transition files" in line:
					continue
				n = int(line.split()[0])
				print(line) 
				print(n)
		with open("%s.def" % M) as defFile:
			for line in defFile:
				if not "Isotopologue mass" in line:
					continue
				mass = line.split()[0]
				print(line) 
				print(mass)
		with open("%s.def" % M) as defFile:
			for line in defFile:
				if not "Default value of Lorentzian half-width for all lines" in line:
					continue
				dL = line.split()[0]
				print(line) 
				print(dL)
		with open("%s.def" % M) as defFile:
			for line in defFile:
				if not "Default value of temperature exponent for all lines" in line:
					continue
				dn = line.split()[0]
				print(line) 
				print(dn)
		with open("%s.def" % M) as defFile:
			for line in defFile:
				if not "Maximum wavenumber" in line:
					continue
				nuMax = line.split()[0]
				print(line) 
				print(nuMax)
		with open("%s.def" % M) as defFile:
			for line in defFile:
				if not "Version number with format" in line:
					continue
				version = line.split()[0]
				print(line) 
				print(version)


	if(n == 1):
		s = int(math.ceil(float(nuMax)))

	warning1 = 0
	if(len(dL) == 0):
		warning1 = 1
		dL = 0.07

	warning2 = 0
	if(len(dn) == 0):
		warning2 = 1
		dn = 0.5

	#correct now wrong number of files
	#if(M == "12C-1H3-37Cl__OYT"):
	#	n=64

	if(n != nn):
		print("Error, number of .trans files do not agree %d %d" % (n, nn))
		return 0

	if(nuMax == 0.0):
		print("Error, Maximum wavenumber in def file not found")
		return 0

	l=np.zeros(n, dtype=int)
	jarray=np.zeros(n+1, dtype=int)
	for nu in range(n + 1):
		jarray[nu] = nu * s

	transFile = ''

	if(DownloadFiles == 3):
		#used for super lines
		n = 1

	for nu in range(n):
		print("------Download file %d from %d -----" % (nu, n))
		if(M == "1H2-16O__BT2"):
			jarray = [0, 250, 500, 750, 1000, 1500, 2000, 2250, 2750, 3500, 4500, 5500, 7000, 9000, 14000, 20000, 30000] 
			transFile = "%s__%05d-%05d.trans" % (M, jarray[nu], jarray[nu+1])
			if(DownloadFiles == 1):
				com = "wget http://www.exomol.com/db/%s/%s.bz2" % (P, transFile)
				er=os.system(com)
				if(er != 0):
					print("Error in download .trans file")
					exit()
				else:
					com = "bzip2 -d %s.bz2" % (transFile)
					er=os.system(com)
			if(getTimeStamp == 1):
				url = "http://www.exomol.com/db/%s/%s.bz2" % (P, transFile)

				#first check if file exists
				com = "wget -S --spider %s" % url
				er=os.system(com)
				#print("er", er)
				if(er == 0):
					date = getTime(url)
					date1 = max(date, date1)


		elif(M == "1H-2H-16O__VTT"):
			jarray = [0, 250, 500, 750, 1000, 1500, 2000, 2250, 2750, 3500, 4500, 5500, 7000, 9000, 14000, 20000, 26000] 
			transFile = "%s__%05d-%05d.trans" % (M, jarray[nu], jarray[nu+1])
			if(DownloadFiles == 1):
				com = "wget http://www.exomol.com/db/%s/%s.bz2" % (P, transFile)
				er=os.system(com)
				if(er != 0):
					print("Error in download .trans file")
					exit()
				else:
					com = "bzip2 -d %s.bz2" % (transFile)
					er=os.system(com)
			if(getTimeStamp == 1):
				url = "http://www.exomol.com/db/%s/%s.bz2" % (P, transFile)

				#first check if file exists
				com = "wget -S --spider %s" % url
				er=os.system(com)
				#print("er", er)
				if(er == 0):
					date = getTime(url)
					date1 = max(date, date1)


		else:

			if(DownloadFiles == 3):
				#super lines
				nuMaxs = int(math.ceil(float(nuMax)))
				transFile = "%s__%05d-%05d__%dK.super" % (M, 0, nuMaxs, Temp)
				print('Download superlines', transFile)
				com = "wget http://www.exomol.com/db/%s/%s.bz2" % (P, transFile)
				url = "http://www.exomol.com/db/%s/%s.bz2" % (P, transFile)

				er=os.system(com)
				if(er == 0):
					com = "bzip2 -d %s.bz2" % transFile
					er=os.system(com)
				if(er != 0):
					print("Error in download .super file")
					exit()
				jarray[1] = nuMaxs
				l[0]=int(subprocess.check_output(['wc', '-l', "%s" % transFile]).split()[0])
			else:
				if(n > 1 or M == "75As-1H3__CYT18"):
					#multiple transition files
					transFile = "%s__%05d-%05d.trans" % (M, nu * s, nu * s + s)
					transFile4 = "%s__%04d-%04d.trans" % (M, nu * s, nu * s + s)
					if(DownloadFiles == 1):
						#trans files with only 4 digits
						if(dg == 4):
							com = "wget http://www.exomol.com/db/%s/%s.bz2" % (P, transFile4)
							url = "http://www.exomol.com/db/%s/%s.bz2" % (P, transFile4)
						elif(dg == 5):
							com = "wget http://www.exomol.com/db/%s/%s.bz2" % (P, transFile)
							url = "http://www.exomol.com/db/%s/%s.bz2" % (P, transFile)
						else:
							print("Error, number of digits not valid")
							return 0

						er=os.system(com)
						if(er != 0):
							print("Error in download .trans file")
							exit()

						if(dg == 4):
							com = "mv %s.bz2 %s.bz2" % (transFile4, transFile)
							er=os.system(com)

						com = "bzip2 -d %s.bz2" % (transFile)
						er=os.system(com)


						if(getTimeStamp == 1):
							#first check if file exists
							com = "wget -S --spider %s" % url
							er=os.system(com)
							#print("er", er)
							if(er == 0):
								date = getTime(url)
								date1 = max(date, date1)
							else:
								url = "http://www.exomol.com/db/%s/%s" % (P, transFile)

								#first check if file exists
								com = "wget -S --spider %s" % url
								er=os.system(com)
								#print("er", er)
								if(er == 0):

									date = getTime(url)
									date1 = max(date, date1)

				else:
					#only 1 transition file
					transFile = "%s.trans" % (M)
					if(DownloadFiles == 1):
						com = "wget http://www.exomol.com/db/%s/%s.bz2" % (P, transFile)
						er=os.system(com)
						if(er == 0):
							com = "bzip2 -d %s.bz2" % transFile
							er=os.system(com)
						else:
							com = "wget http://www.exomol.com/db/%s/%s" % (P, transFile)
							er=os.system(com)
							if(er != 0):
								print("Error in download .trans file")
								exit()
					if(getTimeStamp == 1):
						url = "http://www.exomol.com/db/%s/%s.bz2" % (P, transFile)

						#first check if file exists
						com = "wget -S --spider %s" % url
						er=os.system(com)
						#print("er", er)
						if(er == 0):
							date = getTime(url)
							date1 = max(date, date1)
						else:
							url = "http://www.exomol.com/db/%s/%s" % (P, transFile)

							#first check if file exists
							com = "wget -S --spider %s" % url
							er=os.system(com)
							#print("er", er)
							if(er == 0):

								date = getTime(url)
								date1 = max(date, date1)



		if(PrintISO == 1 and DownloadFiles != 3):
			#check if trans file exists
			exists = os.path.isfile("%s" % transFile)
			if(exists == 0):
				print("Error, transition file not found: %s" % transFile)
				return 0

			l[nu]=int(subprocess.check_output(['wc', '-l', "%s" % transFile]).split()[0])
			print("nu", nu, l[nu])


			#determine the number of columns in the transition files
			with open("%s" % transFile) as tFile:
				line = tFile.readline()
				ntcol = len(line.split())
				print("number of columns in trans file", ntcol)

			#check the maximum wavenumber. Its sometimes wrong in the def files
			if(n == 1 and ntcol == 4):
				nnuMax = 0.0
				with open(transFile) as infile:
					il = 0
					for line in infile:
						#print(line)
						#print("x", line.split()[3])
						if(len(line) > 1):
							nnuMax = max(nnuMax, float(line.split()[ntcol - 1]))
							if(il % 10000 == 0):
								print("nuMax", il, nnuMax)
							il += 1
				print("nuMax:", nuMax, nnuMax)
				jarray[1] = int(math.ceil(float(nnuMax)))

	print("download finished")

	if(warning1 == 1):
		print("Waring, <Default value of Lorentzian half-width for all lines> not found, used dL = %g" % dL)
	if(warning2 == 1):
		print("Waring, <Default value of temperature exponent for all lines> not found, used dn = %g" % dn)

	if(getTimeStamp == 1):
		f = open(("%s.time" % M),'w')
		print(date1, file = f)

	if(PrintISO == 1):
		if(DownloadFiles != 3):
			f = open(("%s.param" % M),'w')
		else:
			f = open(("%s__%05d-%05d__%dK_super.param" % (M, 0, nuMaxs, Temp)),'w')

		if(DownloadFiles != 3):
			lStates=int(subprocess.check_output(['wc', '-l', "%s.states" % M]).split()[0])
		else:
			lStates = 0

		if(DownloadFiles != 3):
			print("Database = 2", file = f)
		else:
			print("Database = 20", file = f)

		print("Molecule number = %d" % 1, file = f)
		if(DownloadFiles != 3):
			print("Name = %s" % M, file = f)
		else:
			print("Name = %s__%05d-%05d__%dK" % (M, 0, nuMaxs, Temp), file = f)

		print("Number of Isotopologues = 1", file = f)
		print("#Id Abundance      Q(296K)   g     Molar Mass(g)  partition file :", file = f)
		print("0 1.0             0.0       0      %s        %s.pf" % (mass, M), file = f)
		print("Number of columns in partition File = %s" % npfcol, file = f)
		print("Number of line/transition files = %d" % n, file = f)
		print("Number of lines per file :", file = f)
		for nu in range(n):
			print("%d" % l[nu], file = f)
		print("Line file limits :", file = f)
		for nu in range(n + 1):
			print("%d" % (jarray[nu]), file = f)
		print("#ExoMol :", file = f)
		print("Number of states = %d" % lStates, file = f)	
		print("Number of columns in transition files = %s" % ntcol, file = f)
		print("Default value of Lorentzian half-width for all lines = %s" % dL, file = f)
		print("Default value of temperature exponent for all lines = %s" % dn, file = f)
		print("Version = %s" % version, file = f)

if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	parser.add_argument('-M', '--M', type=str,
                    help='Species', default = '')
	parser.add_argument('-D', '--DownloadFiles', type=int,
                    help='Download Files', default = 1)
	parser.add_argument('-P', '--PrintISO', type=int,
                    help='Print ISO', default = 1)
	parser.add_argument('-Ts', '--getTimeStamp', type=int,
                    help='get time stamp', default = 0)
	parser.add_argument('-Temp', '--Temperature', type=int,
                    help='Temperature', default = 100)

	args = parser.parse_args()

	M = args.M
	DownloadFiles = args.DownloadFiles
	PrintISO = args.PrintISO
	getTimeStamp = args.getTimeStamp
	Temp = args.Temperature
	if(M == ''):
		print("Error, no species specified, run python exomol.py -M <id>")



	main(M, DownloadFiles, PrintISO, getTimeStamp, Temp)
