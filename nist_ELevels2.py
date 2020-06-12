#https://towardsdatascience.com/controlling-the-web-with-python-6fceb22c5f08
'''
This scripts downloads the energy levels from the NIST database.
The url is copied from a manual search, and Z and I are modified
The script stores the Energy levels in a file NIST_Elevles<Z>>I>.dat

Date: June 2020
Author: Simon Grimm

'''

import sys
import argparse
import csv
import requests

def nist(Z, I):

	print(Z, I)

	CSV_URL = 'https://physics.nist.gov/cgi-bin/ASD/energy1.pl?encodedlist=XXT2&de=0&spectrum=Z+%%3D+%d+%d&submit=Retrieve+Data&units=0&format=3&output=0&page_size=15&multiplet_ordered=0&conf_out=on&term_out=on&level_out=on&unc_out=1&j_out=on&g_out=on&lande_out=on&perc_out=on&biblio=on&temp=' % (Z, I)


	with requests.Session() as s:
		download = s.get(CSV_URL)
		data = download.content.decode('utf-8')


		check = data.find('Invalid charge')
		if(check > -1):
			print('Invalid charge')

		else:

			s1 = data.replace('"', '')
			with open("NIST_ELevels%02d%02d.dat" % (Z, I), "w") as f:
				f.write(s1)

if __name__ == "__main__":


	parser = argparse.ArgumentParser()

	parser.add_argument('-Z', '--Z', type=int,
		help='Z', default = 1)
	parser.add_argument('-I', '--I', type=int,
		help='I', default = 0)
	
	args = parser.parse_args()

	Z = args.Z
	I = args.I

	nist(Z,I)

