# this function scans the isotopologue properties from the hitran webiste
# The website is more complete thant the molparam.txt file and the HAPI interface

#Author: Simon Grimm
#February 2019
import pandas as pd


def main():

	url="http://hitran.org/docs/iso-meta/"


	data = pd.read_html(url, header=0)
	print(len(data))

	abundance= []
	m = []
	Q0 = []
	iso = []
	g = []
	qfile = []
	Formula = []

	for j in range(len(data)):
	#for j in range(3):
		#print(data[j])

		data[j].to_csv('test.csv')

		l= data[j].values.tolist()

		#print(l)
		#print("****")
		#print(l[j])
		#print(l[j][4])


		for i in range(len(l)):
			#Abundance
			a = str(l[i][4])
			#replace E notation
			a1=float(a.replace('\xa0×\xa010', 'E'))
			abundance.append(a1)
	
			#Molar Mass
			m.append(l[i][5])

			#Formula
			Formula.append(l[i][2])


			#Q(296 K)
			q0 = str(l[i][6])
			#replace E notation
			q0 = float(q0.replace('\xa0×\xa010', 'E'))
			
			Q0.append(q0)

			#g
			g.append(l[i][8])

			#iso , molecule + local ID
			index = str(l[i][1])
			if(index == "11"):
				index = "A"
			if(index == "12"):
				index = "B"
			if(j < 9):
				index = " " + str(j+1) + index
			else:
				index = str(j+1) + index
			iso.append(index)
		
			#qfile	
			qfile.append(l[i][7])

	f = open('Hitran_species.dat', 'w')
	for i in range(len(iso)):
		print("%3s %-16.12s %-16.12s %-16s %-16s %-16s %-16s" % (iso[i], abundance[i], m[i], Q0[i], qfile[i], g[i], Formula[i]), file=f)	

	f.close()

if __name__ == '__main__':
	main()

