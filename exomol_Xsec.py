# **********************************************************************************************
# This script downloads all available xsec files for a given temperature T and resolution dnu, and
# includes the linelist name in the filename.
# For some molecules, the temperature has to be changed, because the maximal xsec temperature on 
# the EXOMOL website is not the same as the maximal temperature of the line list data.
# 
# Parameters are: 
#    -T : temperature in K, default 1500
#    -dnu : wavelength resolution in cm^-1, default 0.1
#    -S : species name (optional), if not specified then all files are downloaded
#
# The scripts uses the Firefox browser for the webpage navigation
#
# Date: May 2019
# Author: Simon Grimm
#
# *********************************************************************************************

import argparse
import numpy as np
from bs4 import BeautifulSoup
import requests
import pandas as pd

from selenium import webdriver
from selenium.webdriver.common.keys import Keys

import wget
import os

import exomol2

def main(T, dnu, S):

	url="http://exomol.com/data/molecules/"
	page = requests.get(url).text
	soup = BeautifulSoup(page, "html.parser")

	List = soup.find_all('a', attrs={"class" : "list-group-item link-list-group-item"})

	exfile = "Exomol_xsec_species.dat"


	#check if Exomol file exists and create it otherwise
	exists = os.path.isfile(exfile)
	if(exists == 0):
		exomol2.main()

	xList, xmList = np.loadtxt(exfile, usecols=(1,3), unpack=True, dtype=np.str)

	#for i in range(len(xList)):
	#	print(xList[i], xmList[i]) 

	for i in range(0, len(xList)):

		print(i)
		x = xList[i]
		xm = xmList[i]

		#Download ony one specific molecule:
		if(len(S) > 0):
			if(x != S):
				continue
		
		print("xsec", x, xm)

		
		url = "http://exomol.com/xsec/" + x + "/"


		#scan for the numin and numax values
		page = requests.get(url).text
		soup = BeautifulSoup(page, "html.parser")

		List = soup.find_all('label', attrs={"for" : "id_numin"})

		if(len(List) <= 0):
			print("********Error xsec not available **********")
			continue

		l = str(List[0])
		l = l.split("<sub>min</sub> (")[1]
		l = l.split("cm<sup>-1</sup>")[0]
		numin = int(l.split(" - ")[0])
		numax = int(l.split(" - ")[1])
		print(l, numin, numax)

		#scan for the temperature values
		List = soup.find_all('label', attrs={"for" : "id_T"})

		l = str(List[0])
		l = l.split("<em>T</em> (")[1]
		l = l.split("K)")[0]
		Tmin = int(l.split(" - ")[0])
		Tmax = int(l.split(" - ")[1])

		TT = T
		if(Tmin > T):
			TT = Tmin
		if(Tmax < T):
			TT = Tmax

		#modify here the Temperature for molecules where the maximal temperature of the line list
		#is less than the maximal temperature of the xsec data
		if(x == "1H-14N-16O3"):
			TT = 500
			print("***************** modify T for 1H-14N-16O3 ***************")

		print(l, Tmin, Tmax, TT)


		driver = webdriver.Firefox()
		driver.get(url)

		spectra = driver.find_element_by_name("dnu")
		spectra.send_keys(str(dnu))

		spectra = driver.find_element_by_name("numin")
		spectra.send_keys(numin)

		spectra = driver.find_element_by_name("numax")
		spectra.send_keys(numax - 1)

		spectra = driver.find_element_by_name("T")
		spectra.send_keys(TT)

		spectra = driver.find_element_by_name("spoon_feed")
		spectra.click()

		spectra = driver.find_element_by_xpath("//input[@type='submit' and @value='Submit']")
		spectra.click()

		#filename on the EXOMOL webiste
		filename = ("%s_%d-%d_%dK_%.6f.sigma" % (x, numin, numax - 1, TT, dnu))
		#filename containing the line list name
		xm2 = xm.replace("HITEMP", "HITEMP2010") 
		filename2 = ("%s_%d-%d_%dK_%.6f.sigma" % (xm2, numin, numax - 1, TT, dnu))
		url = "http://exomol.com/results/" + filename
		print(filename)

		file = wget.download(url, filename2)
		print("")

		driver.close()

if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	parser.add_argument('-T', '--T', type=int,
                    help='Temperature', default = 1500)

	parser.add_argument('-dnu', '--dnu', type=float,
                    help='wavelength resolution', default = 0.1)

	parser.add_argument('-S', '--S', type=str,
                    help='Species, e.g. 1H-14N-16O3', default = '')

	args = parser.parse_args()

	T = args.T
	dnu = args.dnu
	S = args.S

	print(T)
	print(dnu)
	print(S)

	main(T, dnu, S)
