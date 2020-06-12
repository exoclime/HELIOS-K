'''
This scripts downloads the line lists from the NIST database.
The url is copied from a manual search, and Z and I are modified
The script stores the Energy levels in a file NIST_Lines<Z>>I>.dat

Date: June 2020
Author: Simon Grimm
'''


import sys
import argparse
import csv
import requests


def Lines(Z, I):

	print(Z, I)

	CSV_URL = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra=Z+%%3D+%d+I+%d&limits_type=0&low_w=&upp_w=&unit=1&de=0&format=2&line_out=1&remove_js=on&en_unit=0&output=0&bibrefs=1&page_size=15&show_obs_wl=1&show_calc_wl=1&show_wn=1&unc_out=1&order_out=0&max_low_enrg=&show_av=3&max_upp_enrg=&tsb_value=0&min_str=&A_out=0&intens_out=on&max_str=&allowed_out=1&forbid_out=1&min_accur=&min_intens=&conf_out=on&term_out=on&enrg_out=on&J_out=on&g_out=on&submit=Retrieve+Data' % (Z, I)




	#ascii
	#CSV_URL = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra=Z+%3D+41&limits_type=0&low_w=&upp_w=&unit=1&de=0&format=1&line_out=0&remove_js=on&en_unit=0&output=0&bibrefs=1&page_size=15&show_obs_wl=1&show_calc_wl=1&show_wn=1&unc_out=1&order_out=0&max_low_enrg=&show_av=5&max_upp_enrg=&tsb_value=0&min_str=&A_out=0&intens_out=on&max_str=&allowed_out=1&forbid_out=1&min_accur=&min_intens=&conf_out=on&term_out=on&enrg_out=on&J_out=on&g_out=on&submit=Retrieve+Data'


	with requests.Session() as s:
		download = s.get(CSV_URL)
		data = download.content.decode('utf-8')


		check = data.find('No lines are available')
		if(check > -1):
			print('No data found')
			
		else:

			#get the number of lines from the data
			#data1 = csv.reader(data.splitlines(), delimiter=',')
			#lines = list(data1)
			#print(len(lines))
			
			#print(data)

			s1 = data.replace('?', '')
			s2 = s1.replace('=', '')
			s3 = s2.replace('[', '')
			s4 = s3.replace(']', '')
			s5 = s4.replace('(', '')
			s6 = s5.replace(')', '')

			with open("NIST_Lines%02d%02d.dat" % (Z, I), "w") as f:
            			f.write(s6)


if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	parser.add_argument('-Z', '--Z', type=int,
		help='Z', default = 1)
	parser.add_argument('-I', '--I', type=int,
		help='I', default = 0)

	args = parser.parse_args()

	Z = args.Z
	I = args.I

	Lines(Z,I)

