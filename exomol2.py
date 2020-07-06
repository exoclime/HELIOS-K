# **********************************************************************************************
# This code scans the EXOMOL webiste, and lists all available molecules and linelists.
# It writes to files: "Exomol_species.dat" and "Exomol_xsec_species.dat", and lists the
# corresponding file names

# Date: May 2019
# Author: Simon Grimm
#
# *********************************************************************************************

from bs4 import BeautifulSoup
import requests
import sys

#this function extracts the ranges of the .trans files
#it returns the ranges, or -1 when the ranges are not equal
#it returns the number of transition files
#it returns the number of digits of the ranges
def transitionRanges(url):
	#url="http://exomol.com/data/molecules/H2O/1H2-16O/BT2/"

	page = requests.get(url).text
	soup = BeautifulSoup(page, "html.parser")
	List = soup.find_all('li', attrs={"class" : "list-group-item link-list-group-item"})

	#print(List[4].a)
	#print(List[4].a.get('href'))

	transList = []

	#write a list with all transition file ranges
	for i in range(len(List)):
		el = List[i].a.get('href')
		#print(el)
		el1 = el.split('__')[-1]        #split at __ and take right part of it
		el2 = el1.split('.trans')[0]
		if(len(el1.split('.trans')) > 1):
			#print(el2)
			transList.append(el2)

	rangesList = []
	if(len(transList) > 1):
		#check range of files
		for x in transList:
			x0 = float(x.split('-')[0])
			x1 = float(x.split('-')[1])
			dg = len(x.split('-')[0])
			#print(x1-x0)
			rangesList.append(x1-x0)

		s = rangesList[0]
		for r in rangesList:
			if(r != s):
				s=-1
		n = len(rangesList)
	else:
		s = 0
		n = 1
		dg = 0
	#print(s, n)
	return(s, n, dg)


def main():

	print("Scan Exomol webiste for file names")

	url="http://exomol.com/data/molecules/"
	page = requests.get(url).text
	soup = BeautifulSoup(page, "html.parser")

	List = soup.find_all('a', attrs={"class" : "list-group-item link-list-group-item molecule_link"})

	efile = open("Exomol_species.dat", "w", buffering=1)
	exfile = open("Exomol_xsec_species.dat", "w", buffering=1)

	if(len(List) == 0):
		print("Error, no molecules found, maybe the Exomol homepage has changed")
		sys.exit(100)

	#Molecule
	for i in range(len(List)):
	#for i in range(20):
		el = List[i].get('href')
		print(el)

		url1 = url + el + "/"
		page1 = requests.get(url1).text
		soup1 = BeautifulSoup(page1, "html.parser")

		List1 = soup1.find_all('a', attrs={"class" : "list-group-item link-list-group-item"})

		#Isotopologue
		for j in range(len(List1)):
			el1 = List1[j].get('href')

			print("    ", el1)

			url2 = url1 + el1 + "/"
			page2 = requests.get(url2).text
			soup2 = BeautifulSoup(page2, "html.parser")

			List2 = soup2.find_all('a', attrs={"class" : "list-group-item link-list-group-item "})
			List2 += soup2.find_all('a', attrs={"class" : "list-group-item link-list-group-item recommended"})

			#Line list
			for k in range(len(List2)):
				el2 = List2[k].get('href')
				
				el3 = el2.replace("HITEMP", "HITEMP2010") 
				
				if(el2.find("xsec-") >= 0):

					#change HITEMP to HITEMP2010, but only in the prints
					
					name = el2.split("xsec-")[1]
					p3 = el1 + "__" + name
					print("%-16s %-24s %-32s %-32s" % (el, el1, el3, p3), file=exfile)

				else:
					print("        ", el1 + "__" + el3, el + "/" + el1 + "/" + el2 )

					p2 = el1 + "__" + el3
					p3 = el + "/" + el1 + "/" + el2

					url3 = url2 + el2
					#print(url3)
					s, n, dg = transitionRanges(url3)
					print(s, n, dg)
					print("%-16s %-24s %-32s %-40s %8g %8g %8g" % (el, el1, p2, p3, s, n, dg), file=efile)


	efile.close()
	print("Scan complete")

if __name__ == '__main__':
	main()
