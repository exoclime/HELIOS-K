import numpy as np
from bs4 import BeautifulSoup
import requests

#head -n -1 *.ref > test


exfile = "Exomol_species.dat"
M0, P0 = np.loadtxt(exfile, usecols=(2,3), unpack=True, dtype=str)

#M = "1H2-16O__BT2"
#P = P0[M0 == M][0]
#name = "BT2"

for m in range(len(M0)):
#for m in range(81, 82):
	M = M0[m]
	P = P0[m]
	name = M.split("__")[1]

	print("%s %s" % (M, P))

	filename = "%s.ref" % M
	f = open(filename,'w')


	url = "http://exomol.com/data/molecules/%s" % P

	page = requests.get(url).text
	soup = BeautifulSoup(page, "html.parser")

	#print(soup.text)

	for i in range(3):

		if(i == 0): 
			List = soup.find_all(text='%s: line list' % name)
			if(len(List)> 0):
				print("*** Line list ***", file = f)
				print("-----------------", file = f)
			else:
				List = soup.find_all(text='%s: line list (external)' % name)
				if(len(List)> 0):
					print("*** Line list ***", file = f)
					print("-----------------", file = f)
		if(i == 1): 
			List = soup.find_all(text='%s: energy levels' % name)
			if(len(List)> 0):
				print("*** Energy levels ***", file = f)
				print("---------------------", file = f)
			else:
				List = soup.find_all(text='%s: energy levels (external)' % name)
				if(len(List)> 0):
					print("*** Energy levels ***", file = f)
					print("---------------------", file = f)
		if(i == 2): 
			List = soup.find_all(text='%s: partition function' % name)
			if(len(List)> 0):
				print("*** Partition function ***", file = f)
				print("--------------------------", file = f)
			else:
				List = soup.find_all(text='%s: partition function (external)' % name)
				if(len(List)> 0):
					print("*** Partition function ***", file = f)
					print("--------------------------", file = f)
		if(i == 3): 
			List = soup.find_all(text='%s: broadening coefficients' % name)
			if(len(List)> 0):
				print("*** Broadening coefficients ***", file = f)
				print("-------------------------------", file = f)
			else:
				List = soup.find_all(text='%s: broadening coefficients (external)' % name)
				if(len(List)> 0):
					print("*** Broadening coefficients ***", file = f)
					print("-------------------------------", file = f)

		#print(List)
		#print(" length of list %d" % len(List))

		if(len(List) > 0):
			#find the References title
			text1 = List[0].findNext(text='References')

			#find the list of entries
			text = text1.findNext('ol')

			#find all entries
			text2 = text.find_all('li')

			for t in text2:
				print(t.text, file = f)

				ref = t.find_all('a', href=True)
				if(len(ref) > 0):
					print(ref[0]['href'], file = f)

				print("", file = f)
		



