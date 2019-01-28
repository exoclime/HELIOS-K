# This script downloads and unpacks  the *.states, *.pf, *.trans
# and *.def files from wwww.exomol.com.
# And it generates the < species name >.param files

# July 2018
# Author: Simon Grimm

#run with "python exomol.py id" where id is the number of the molecule

import sys
import os
import subprocess
import numpy as np

def main(m):


	print("Molecule = %s" % m)

	PrintISO=1			#when set to 1, then print the code fot the ISO.h file
	DownloadFiles=1			#0 no download
					#1 download all file
					#2 download only .def files

	M = ""
	P = ""
	s = 0		#file range
	ntcol = 0	#columns in transition files
	npfcol = 0	#columns in partition function file




	if(m == 1):
		#1 H2O
		M="1H2-16O__BT2"
		P="H2O/1H2-16O/BT2"
		s=-1
		ntcol=3
		npfcol=3

	if(m == 1.01):
		#the number of transition files is wrong, correct later
		#1 H2O
		M="1H2-16O__POKAZATEL"
		P="H2O/1H2-16O/POKAZATEL/"
		s=100
		ntcol=3
		npfcol=2

	if(m == 1.02):
		#the number of transition files is wrong, correct later
		#1 H2O 2
		M="1H2-18O__HotWat78"
		P="H2O/1H2-18O/HotWat78"
		s=1000
		ntcol=3
		npfcol=2

	if(m == 1.03):
		#the number of transition files is wrong, correct later
		#1 H2O 3
		M="1H2-17O__HotWat78"
		P="H2O/1H2-17O/HotWat78"
		s=1000
		ntcol=3
		npfcol=2

	if(m == 1.04):
		#the number of transition files is wrong, correct later
		#1 H2O 4
		M="1H-2H-16O__VTT"
		P="H2O/1H-2H-16O/VTT"
		s=1000
		ntcol=3
		npfcol=2

	if(m == 5):
		#5 CO
		M="12C-16O__Li2015"
		P="CO/12C-16O/Li2015"
		s=22000
		ntcol=3
		npfcol=2

	if(m == 6):
		#6 CH4
		M="12C-1H4__YT10to10"
		P="CH4/12C-1H4/YT10to10"
		s=100
		ntcol=3
		npfcol=2

	if(m == 8):
		#8 NO
		M="14N-16O__NOname"
		P="NO/14N-16O/NOname"
		s=40000
		ntcol=4
		npfcol=2

	if(m == 8.02):
		#8 NO 2
		M="15N-16O__NOname"
		P="NO/15N-16O/NOname"
		s=40000
		ntcol=4
		npfcol=2

	if(m == 8.03):
		#8 NO 3
		M="14N-18O__NOname"
		P="NO/14N-18O/NOname"
		s=40000
		ntcol=4
		npfcol=2

	if(m == 8.04):
		#8 NO 4
		M="15N-18O__NOname"
		P="NO/15N-18O/NOname"
		s=40000
		ntcol=4
		npfcol=2

	if(m == 8.05):
		#8 NO 5
		M="14N-17O__NOname"
		P="NO/14N-17O/NOname"
		s=40000
		ntcol=4
		npfcol=2

	if(m == 8.06):
		#8 NO 6
		M="15N-17O__NOname"
		P="NO/15N-17O/NOname"
		s=40000
		ntcol=4
		npfcol=2

	if(m == 9):
		#9 SO2
		M="32S-16O2__ExoAmes"
		P="SO2/32S-16O2/ExoAmes"
		s=100
		ntcol=3
		npfcol=2

	if(m == 11):
		#11 NH3
		M="14N-1H3__BYTe"
		P="NH3/14N-1H3/BYTe"
		s=100
		ntcol=3
		npfcol=2

	if(m == 11.02):
		#11 NH3 2
		M="15N-1H3__BYTe-15"
		P="NH3/15N-1H3/BYTe-15"
		s=100
		ntcol=3
		npfcol=2

	if(m == 12):
		#12 NHO3
		M="1H-14N-16O3__AIJS"
		P="HNO3/1H-14N-16O3/AIJS"
		s=100
		ntcol=3
		npfcol=2

	if(m == 15):
		#15 NCl
		M="1H-35Cl__Yueqi"
		P="HCl/1H-35Cl/Yueqi"
		s=2400
		ntcol=4
		npfcol=2

	if(m == 20):
		#20 H2CO
		M="1H2-12C-16O__AYTY"
		P="H2CO/1H2-12C-16O/AYTY"
		s=100
		ntcol=3
		npfcol=2

	if(m == 23):
		#23 HCN
		M="1H-12C-14N__Harris"
		P="HCN/1H-12C-14N/Harris"
		s=17586			#file range
		ntcol=4
		npfcol=2

	if(m == 23.02):
		#23 HCN 2
		M="1H-13C-14N__Larner"
		P="HCN/1H-13C-14N/Larner"
		s=17596			#file range
		ntcol=4
		npfcol=2

	if(m == 24):
		#24 CH3Cl
		M="12C-1H3-35Cl__OYT"
		P="CH3Cl/12C-1H3-35Cl/OYT"
		s=100
		ntcol=3
		npfcol=2

	if(m == 24.02):
		#24 CH3Cl 2
		M="12C-1H3-37Cl__OYT"
		P="CH3Cl/12C-1H3-37Cl/OYT"
		s=100
		ntcol=3
		npfcol=2

	if(m == 25):
		#25 H2O2
		M="1H2-16O2__APTY"
		P="H2O2/1H2-16O2/APTY"
		s=100
		ntcol=3
		npfcol=2

	if(m == 28):
		#28 PH3
		M="31P-1H3__SAlTY"
		P="PH3/31P-1H3/SAlTY"
		s=100
		ntcol=3
		npfcol=2

	if(m == 31):
		#31 H2S
		M="1H2-32S__AYT2"
		P="H2S/1H2-32S/AYT2"
		s=1000			#file range
		ntcol=3
		npfcol=2

	if(m == 38):
		#38 C2H4
		M="12C2-1H4__MaYTY"
		P="C2H4/12C2-1H4/MaYTY"
		s=100			#file range
		ntcol=3
		npfcol=2

	if(m == 46):
		#46 CS
		M="12C-32S__JnK"
		P="CS/12C-32S/JnK"
		s=11000			#file range
		ntcol=4
		npfcol=2

	if(m == 46.02):
		#46 CS 2
		M="12C-34S__JnK"
		P="CS/12C-34S/JnK"
		s=11000			#file range
		ntcol=4
		npfcol=2

	if(m == 46.03):
		#46 CS 3
		M="13C-32S__JnK"
		P="CS/13C-32S/JnK"
		s=11000			#file range
		ntcol=4
		npfcol=2

	if(m == 46.04):
		#46 CS 4
		M="12C-33S__JnK"
		P="CS/12C-33S/JnK"
		s=11000			#file range
		ntcol=4
		npfcol=2

	if(m == 46.05):
		#46 CS 5
		M="12C-36S__JnK"
		P="CS/12C-36S/JnK"
		s=11000			#file range
		ntcol=4
		npfcol=2

	if(m == 46.06):
		#46 CS 6
		M="13C-33S__JnK"
		P="CS/13C-33S/JnK"
		s=11000			#file range
		ntcol=4
		npfcol=2

	if(m == 46.07):
		#46 CS 7
		M="13C-34S__JnK"
		P="CS/13C-34S/JnK"
		s=11000			#file range
		ntcol=4
		npfcol=2

	if(m == 46.08):
		#46 CS 8
		M="13C-36S__JnK"
		P="CS/13C-36S/JnK"
		s=11000			#file range
		ntcol=4
		npfcol=2

	if(m == 47):
		#47 SO3
		M="32S-16O3__UYT2"
		P="SO3/32S-16O3/UYT2"
		s=10
		ntcol=3
		npfcol=2

	if(m == 78):
		#78 CH3F
		M="12C-1H3-19F__OYKYT"
		P="CH3F/12C-1H3-19F/OYKYT"
		s=100
		ntcol=4
		npfcol=2

	if(m == 79):
		#79 SiH4
		M="28Si-1H4__OY2T"
		P="SiH4/28Si-1H4/OY2T"
		s=100
		ntcol=3
		npfcol=2

	if(m == 80):
		#80 VO
		M="51V-16O__VOMYT"
		P="VO/51V-16O/VOMYT"
		s=5000			#file range
		ntcol=4
		npfcol=2

		#81 TiO

	if(m == 82):
		#82 FeH
		M="56Fe-1H__Yueqi"
		P="FeH/56Fe-1H/Yueqi"
		s=15000			#file range
		ntcol=3
		npfcol=2

	if(m == 83):
		#83 AlO
		M="27Al-16O__ATP"
		P="AlO/27Al-16O/ATP"
		s=35000			#file range
		ntcol=3
		npfcol=2

	if(m == 83.02):
		#83 AlO 2
		M="26Al-16O__ATP"
		P="AlO/26Al-16O/ATP"
		s=35000			#file range
		ntcol=4
		npfcol=2

	if(m == 83.03):
		#83 AlO 3
		M="27Al-17O__ATP"
		P="AlO/27Al-17O/ATP"
		s=35000			#file range
		ntcol=4
		npfcol=2

	if(m == 83.04):
		#83 AlO 4
		M="27Al-18O__ATP"
		P="AlO/27Al-18O/ATP"
		s=35000			#file range
		ntcol=4
		npfcol=2

	if(m == 84):
		#84 SiO
		M="28Si-16O__EBJT"
		P="SiO/28Si-16O/EBJT"
		s=6050			#file range
		ntcol=4
		npfcol=2

	if(m == 84.02):
		#84 SiO 2
		M="28Si-17O__EBJT"
		P="SiO/28Si-17O/EBJT"
		s=5939			#file range
		ntcol=4
		npfcol=2

	if(m == 84.03):
		#84 SiO 3
		M="28Si-18O__EBJT"
		P="SiO/28Si-18O/EBJT"
		s=5838			#file range
		ntcol=4
		npfcol=2

	if(m == 84.04):
		#84 SiO 4
		M="29Si-16O__EBJT"
		P="SiO/29Si-16O/EBJT"
		s=6013			#file range
		ntcol=4
		npfcol=2

	if(m == 84.05):
		#84 SiO 5
		M="30Si-16O__EBJT"
		P="SiO/30Si-16O/EBJT"
		s=5978			#file range
		ntcol=4
		npfcol=2

	if(m == 85):
		#85 CaO
		M="40Ca-16O__VBATHY"
		P="CaO/40Ca-16O/VBATHY"
		s=25000			#file range
		ntcol=4
		npfcol=2

	if(m == 86):
		#101 SiH
		M="28Si-1H__SiGHTLY"
		P="SiH/28Si-1H/SiGHTLY"
		s=32000			#file range
		ntcol=4
		npfcol=2

	if(m == 86.02):
		#101 SiH 2
		M="29Si-1H__SiGHTLY"
		P="SiH/29Si-1H/SiGHTLY"
		s=32000			#file range
		ntcol=4
		npfcol=2

	if(m == 86.03):
		#101 SiH 3
		M="30Si-1H__SiGHTLY"
		P="SiH/30Si-1H/SiGHTLY"
		s=32000			#file range
		ntcol=4
		npfcol=2

	if(m == 86.04):
		#101 SiH 4
		M="28Si-2H__SiGHTLY"
		P="SiH/28Si-2H/SiGHTLY"
		s=32000			#file range
		ntcol=4
		npfcol=2

	if(m == 87):
		#87 caH
		M="40Ca-1H__Yadin"
		P="CaH/40Ca-1H/Yadin"
		s=15278			#file range
		ntcol=4
		npfcol=2

	if(m == 88):
		#88 H3+
		M="1H3_p__MiZATeP"
		P="H3_p/1H3_p/MiZATeP"
		s=25000			#file range
		ntcol=3
		npfcol=3

	if(m == 88.02):
		#88 H3+ 2
		M="1H2-2H_p__ST"
		P="H3_p/1H2-2H_p/ST"
		s=18500			#file range
		ntcol=3
		npfcol=2

	if(m == 89):
		#89 PO
		M="31P-16O__POPS"
		P="PO/31P-16O/POPS"
		s=12000			#file range
		ntcol=4
		npfcol=2

	if(m == 90):
		#90 MgH
		M="24Mg-1H__Yadin"
		P="MgH/24Mg-1H/Yadin"
		s=11000			#file range
		ntcol=4
		npfcol=2

	if(m == 90.02):
		#90 MgH 2
		M="25Mg-1H__Yadin"
		P="MgH/25Mg-1H/Yadin"
		s=11000			#file range
		ntcol=4
		npfcol=2

	if(m == 90.03):
		#90 MgH 3
		M="26Mg-1H__Yadin"
		P="MgH/26Mg-1H/Yadin"
		s=11000			#file range
		ntcol=4
		npfcol=2

	if(m == 90.04):
		#90 MgH 4
		M="24Mg-1H__Yueqi"
		P="MgH/24Mg-1H/Yueqi"
		s=29000			#file range
		ntcol=3
		npfcol=2

	if(m == 91):
		#91 NaH
		M="23Na-1H__Rivlin"
		P="NaH/23Na-1H/Rivlin"
		s=32147			#file range
		ntcol=4
		npfcol=2

	if(m == 91.02):
		#91 NaH
		M="23Na-2H__Rivlin"
		P="NaH/23Na-2H/Rivlin"
		s=32147			#file range
		ntcol=4
		npfcol=2

	if(m == 92):
		#92 AlH
		M="27Al-1H__AlHambra"
		P="AlH/27Al-1H/AlHambra"
		s=28000			#file range
		ntcol=4
		npfcol=2

	if(m == 92.02):
		#92 AlH 2
		M="27Al-2H__AlHambra"
		P="AlH/27Al-2H/AlHambra"
		s=28000			#file range
		ntcol=4
		npfcol=2

	if(m == 92.03):
		#92 AlH 3
		M="26Al-1H__AlHambra"
		P="AlH/26Al-1H/AlHambra"
		s=28000			#file range
		ntcol=4
		npfcol=2

	if(m == 93):
		#93 CrH
		M="52Cr-1H__Yueqi"
		P="CrH/52Cr-1H/Yueqi"
		s=15000			#file range
		ntcol=4
		npfcol=2

	if(m == 94):
		#94 BeH
		M="9Be-1H__Yadin"
		P="BeH/9Be-1H/Yadin"
		s=17000			#file range
		ntcol=4
		npfcol=2

	if(m == 94.01):
		#94 BeH
		M="9Be-1H__Darby-Lewis"
		P="BeH/9Be-1H/Darby-Lewis"
		s=40000			#file range
		ntcol=4
		npfcol=2

	if(m == 94.02):
		#94 BeH 2
		M="9Be-2H__Darby-Lewis"
		P="BeH/9Be-2H/Darby-Lewis"
		s=40000			#file range
		ntcol=4
		npfcol=2

	if(m == 94.03):
		#94 BeH 3
		M="9Be-3H__Darby-Lewis"
		P="BeH/9Be-3H/Darby-Lewis"
		s=40000			#file range
		ntcol=4
		npfcol=2

	if(m == 95):
		#95 TiH
		M="48Ti-1H__Yueqi"
		P="TiH/48Ti-1H/Yueqi"
		s=24000			#file range
		ntcol=4
		npfcol=2

	if(m == 96):
		#96 LiH
		M="7Li-1H__CLT"
		P="LiH/7Li-1H/CLT"
		s=20000			#file range
		ntcol=4
		npfcol=2

	if(m == 97):
		#97 ScH
		M="45Sc-1H__LYT"
		P="ScH/45Sc-1H/LYT"
		s=16000			#file range
		ntcol=4
		npfcol=2

	if(m == 98):
		#98 NiH
		M=""
		P=""
		s=16000			#file range
		ntcol=4
		npfcol=2

	if(m == 99):
		#99 NH
		M="14N-1H__Yueqi"
		P="NH/14N-1H/Yueqi"
		s=17000			#file range
		ntcol=4
		npfcol=2

	if(m == 100):
		#100 CH
		M="12C-1H__Yueqi"
		P="CH/12C-1H/Yueqi"
		s=40000			#file range
		ntcol=4
		npfcol=2

	if(m == 100.02):
		#100 CH 2
		M="13C-1H__Yueqi"
		P="CH/13C-1H/Yueqi"
		s=40000			#file range
		ntcol=4
		npfcol=2

	if(m == 101):
		#101 SH
		M="32S-1H__SNaSH"
		P="SH/32S-1H/SNaSH"
		s=28000			#file range
		ntcol=4
		npfcol=2

	if(m == 101.02):
		#101 SH 2
		M="33S-1H__SNaSH"
		P="SH/33S-1H/SNaSH"
		s=28000			#file range
		ntcol=4
		npfcol=2

	if(m == 101.03):
		#101 SH 3
		M="34S-1H__SNaSH"
		P="SH/34S-1H/SNaSH"
		s=28000			#file range
		ntcol=4
		npfcol=2

	if(m == 101.04):
		#101 SH 4
		M="36S-1H__SNaSH"
		P="SH/36S-1H/SNaSH"
		s=28000			#file range
		ntcol=4
		npfcol=2

	if(m == 101.05):
		#101 SH 5
		M="32S-2H__SNaSH"
		P="SH/32S-2H/SNaSH"
		s=28000			#file range
		ntcol=4
		npfcol=2

	if(m == 102):
		#102 PN
		M="31P-14N__YYLT"
		P="PN/31P-14N/YYLT"
		s=7000			#file range
		ntcol=4
		npfcol=2

	if(m == 102.02):
		#103 PN 2
		M="31P-15N__YYLT"
		P="PN/31P-15N/YYLT"
		s=7000			#file range
		ntcol=4
		npfcol=2

	if(m == 103):
		#103 KCl
		M="39K-35Cl__Barton"
		P="KCl/39K-35Cl/Barton"
		s=3000			#file range
		ntcol=4
		npfcol=2

	if(m == 103.02):
		#103 KCl 2
		M="39K-37Cl__Barton"
		P="KCl/39K-37Cl/Barton"
		s=3000			#file range
		ntcol=4
		npfcol=2

	if(m == 103.03):
		#103 KCl 3
		M="41K-35Cl__Barton"
		P="KCl/41K-35Cl/Barton"
		s=3000			#file range
		ntcol=4
		npfcol=2

	if(m == 103.04):
		#103 KCl 4
		M="41K-37Cl__Barton"
		P="KCl/41K-37Cl/Barton"
		s=3000			#file range
		ntcol=4
		npfcol=2

	if(m == 104):
		#104 NaCl
		M="23Na-35Cl__Barton"
		P="NaCl/23Na-35Cl/Barton"
		s=3000			#file range
		ntcol=4
		npfcol=2

	if(m == 104.02):
		#104 NaCl 2
		M="23Na-37Cl__Barton"
		P="NaCl/23Na-37Cl/Barton"
		s=3000			#file range
		ntcol=4
		npfcol=2
	if(m == 105):
		#105 CN
		M="12C-14N__Yueqi"
		P="CN/12C-14N/Yueqi"
		s=53000			#file range
		ntcol=3
		npfcol=2

	if(m == 106):
		#106 C2
		M="12C2__8states"
		P="C2/12C2/8states"
		s=49000			#file range
		ntcol=4
		npfcol=2

	if(m == 106.02):
		#106 C2 2
		M="12C-13C__8states"
		P="C2/12C-13C/8states"
		s=49000			#file range
		ntcol=4
		npfcol=2

	if(m == 106.03):
		#106 C2 3
		M="13C2__8states"
		P="C2/13C2/8states"
		s=49000			#file range
		ntcol=4
		npfcol=2

	if(m == 108):
		#108 CP
		M="12C-31P__Yueqi"
		P="CP/12C-31P/Yueqi"
		s=16000			#file range
		ntcol=4
		npfcol=2

	if(m == 109):
		#109 PS
		M="31P-32S__POPS"
		P="PS/31P-32S/POPS"
		s=37000			#file range
		ntcol=4
		npfcol=2

	if(m == 110):
		#110 NS
		M="14N-32S__SNaSH"
		P="NS/14N-32S/SNaSH"
		s=39000			#file range
		ntcol=4
		npfcol=2

	if(m == 110.02):
		#110 NS 2
		M="15N-32S__SNaSH"
		P="NS/15N-32S/SNaSH"
		s=39000			#file range
		ntcol=4
		npfcol=2

	if(m == 110.03):
		#110 NS 3
		M="14N-33S__SNaSH"
		P="NS/14N-33S/SNaSH"
		s=39000			#file range
		ntcol=4
		npfcol=2

	if(m == 110.04):
		#110 NS 3
		M="14N-34S__SNaSH"
		P="NS/14N-34S/SNaSH"
		s=39000			#file range
		ntcol=4
		npfcol=2

	if(m == 110.05):
		#110 NS 4
		M="14N-36S__SNaSH"
		P="NS/14N-36S/SNaSH"
		s=39000			#file range
		ntcol=4
		npfcol=2


	if(m == 111):
		#111 SiS
		M="28Si-32S__UCTY"
		P="SiS/28Si-32S/UCTY"
		s=4000			#file range
		ntcol=4
		npfcol=2

	if(m == 111.02):
		#111 SiS 2
		M="28Si-33S__UCTY"
		P="SiS/28Si-33S/UCTY"
		s=4000			#file range
		ntcol=4
		npfcol=2

	if(m == 111.03):
		#111 SiS 3
		M="28Si-34S__UCTY"
		P="SiS/28Si-34S/UCTY"
		s=4000			#file range
		ntcol=4
		npfcol=2

	if(m == 111.04):
		#111 SiS 4
		M="28Si-36S__UCTY"
		P="SiS/28Si-36S/UCTY"
		s=4000			#file range
		ntcol=4
		npfcol=2

	if(m == 111.05):
		#111 SiS 5
		M="29Si-32S__UCTY"
		P="SiS/29Si-32S/UCTY"
		s=4000			#file range
		ntcol=4
		npfcol=2

	if(m == 111.06):
		#111 SiS 6
		M="29Si-34S__UCTY"
		P="SiS/29Si-34S/UCTY"
		s=4000			#file range
		ntcol=4
		npfcol=2

	if(m == 111.07):
		#111 SiS 7
		M="29Si-36S__UCTY"
		P="SiS/29Si-36S/UCTY"
		s=4000			#file range
		ntcol=4
		npfcol=2

	if(m == 111.08):
		#111 SiS 8
		M="30Si-32S__UCTY"
		P="SiS/30Si-32S/UCTY"
		s=4000			#file range
		ntcol=4
		npfcol=2

	if(m == 111.09):
		#111 SiS
		M="30Si-33S__UCTY"
		P="SiS/30Si-33S/UCTY"
		s=4000			#file range
		ntcol=4
		npfcol=2

	if(m == 111.10):
		#111 SiS 10
		M="30Si-34S__UCTY"
		P="SiS/30Si-34S/UCTY"
		s=4000			#file range
		ntcol=4
		npfcol=2

	if(m == 111.11):
		#111 SiS 11
		M="30Si-36S__UCTY"
		P="SiS/30Si-36S/UCTY"
		s=4000			#file range
		ntcol=4
		npfcol=2

	if(m == 111.12):
		#111 SiS 12
		M="29Si-33S__UCTY"
		P="SiS/29Si-33S/UCTY"
		s=4000			#file range
		ntcol=4
		npfcol=2

	if(m == 112):
		#112 HeH+
		M="4He-1H_p__Engel"
		P="HeH_p/4He-1H_p/Engel"
		s=15000			#file range
		ntcol=4
		npfcol=2

	if(m == 112.02):
		#112 HeH+ 2
		M="3He-1H_p__Engel"
		P="HeH_p/3He-1H_p/Engel"
		s=15000			#file range
		ntcol=4
		npfcol=2

	if(m == 112.03):
		#112 HeH+ 3
		M="4He-2H_p__Engel"
		P="HeH_p/4He-2H_p/Engel"
		s=16000			#file range
		ntcol=4
		npfcol=2

	if(m == 112.04):
		#112 HeH+ 4
		M="3He-2H_p__Engel"
		P="HeH_p/3He-2H_p/Engel"
		s=16000			#file range
		ntcol=4
		npfcol=2

	if(m == 113):
		#113 LiCl
		M="7Li-35Cl__Bernath"
		P="LiCl/7Li-35Cl/Bernath"
		s=5000			#file range
		ntcol=4
		npfcol=2

	if(m == 113.02):
		#113 LiCl 2
		M="7Li-37Cl__Bernath"
		P="LiCl/7Li-37Cl/Bernath"
		s=5000			#file range
		ntcol=4
		npfcol=2

	if(m == 113.03):
		#113 LiCl 3
		M="6Li-35Cl__Bernath"
		P="LiCl/6Li-35Cl/Bernath"
		s=6000			#file range
		ntcol=4
		npfcol=2

	if(m == 113.04):
		#113 LiCl 4
		M="6Li-37Cl__Bernath"
		P="LiCl/6Li-37Cl/Bernath"
		s=6000			#file range
		ntcol=4
		npfcol=2

	#if(m == 114):
		#114 NaF 
	#if(m == 115):
		#115 AlCl 
	#if(m == 116):
		#116 AlF 
	#if(m == 117):
		#117 KF 

	if(m == 118):
		#118 LiF 
		M="7Li-19F__Bernath"
		P="LiF/7Li-19F/Bernath"
		s=1900			#file range
		ntcol=4
		npfcol=2

	if(m == 118.02):
		#118 LiF 2
		M="6Li-19F__Bernath"
		P="LiF/6Li-19F/Bernath"
		s=1900			#file range
		ntcol=4
		npfcol=2

	#if(m == 119):
		#119 CaF 
	#if(m == 120):
		#120 MgF 

	print(M)

	if(DownloadFiles == 2):
		if(s!= 0):
			com = "wget http://exomol.com/db/%s/%s.def" % (P, M)
			er=os.system(com)
			com = "wget http://exomol.com/db/%s/%s.pf" % (P, M)
			er=os.system(com)

		return 0

	if(DownloadFiles == 1):
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
			if not "Version number with format" in line:
				continue
			version = line.split()[0]
			print(line) 
			print(version)

	#correct now wrong number of files
	if(m == 38):
		n=71
	if(m == 24.02):
		n=64


	l=np.zeros(n, dtype=int)

	for nu in range(n):
		print("------Download file %d from %d -----" % (nu, n))
		if(m == 1):
			jarray = [0, 250, 500, 750, 1000, 1500, 2000, 2250, 2750, 3500, 4500, 5500, 7000, 9000, 14000, 20000, 30000] 
			if(DownloadFiles == 1):
				com = "wget http://www.exomol.com/db/%s/%s\_\_%05d-%05d.trans.bz2" % (P, M, jarray[nu], jarray[nu+1])
				er=os.system(com)
				if(er != 0):
					print("Error in download .trans file")
					exit()
				else:
					com = "bzip2 -d %s\_\_%05d-%05d.trans.bz2" % (M, jarray[nu], jarray[nu+1])
					er=os.system(com)

			com = "wc -l %s\_\_%d-%d.trans" % (M, jarray[nu], jarray[nu+1])
			l=os.system(com)
			#l[nu]=`wc -l < $M\_\_${jarray[$nu]}-${jarray[$nu + 1]}.trans | awk '{print $1}'`
		elif(m == 104):
			jarray = [0, 250, 500, 750, 1000, 1500, 2000, 2250, 2750, 3500, 4500, 5500, 7000, 9000, 14000, 20000, 26000] 
			if(DownloadFiles == 1):
				com = "wget http://www.exomol.com/db/%s/%s\_\_%05d-%05d.trans.bz2" % (P, M, jarray[nu], jarray[nu+1])
				er=os.system(com)
				if(er != 0):
					print("Error in download .trans file")
					exit()
				else:
					com = "bzip2 -d %s\_\_%05d-%05d.trans.bz2" % (M, jarray[nu], jarray[nu+1])
					er=os.system(com)

			l[nu]=int(subprocess.check_output(['wc', '-l', "%s__%05d-%05d.trans" % (M, nu * s, nu * s + s)]).split()[0])
		else:

			if(n > 1):
				if(DownloadFiles == 1):
					#trans files with inly 4 digits
					if(m == -1):
						com = "wget http://www.exomol.com/db/%s/%s\_\_%04d-%04d.trans.bz2" % (P, M, nu * s, nu * s + s)
						er=os.system(com)
						if(er != 0):
							print("Error in download .trans file")
							exit()
						else:
							com = "bzip2 -d %s\_\_%04d-%04d.trans.bz2" % (M, nu * s, nu * s + s)
							er=os.system(com)
						com = "mv %s__%04d-%04d.trans %s__%05d-%05d.trans" % (M, nu * s, nu * s + s, M, nu * s, nu * s + s)
						er=os.system(com)
					else:
						com = "wget http://www.exomol.com/db/%s/%s\_\_%05d-%05d.trans.bz2" % (P, M, nu * s, nu * s + s)
						er=os.system(com)
						if(er != 0):
							print("Error in download .trans file")
							exit()
						else:
							com = "bzip2 -d %s\_\_%05d-%05d.trans.bz2" % (M, nu * s, nu * s + s)
							er=os.system(com)
				l[nu]=int(subprocess.check_output(['wc', '-l', "%s__%05d-%05d.trans" % (M, nu * s, nu * s + s)]).split()[0])
			else:
				if(DownloadFiles == 1):
					com = "wget http://www.exomol.com/db/%s/%s.trans.bz2" % (P, M)
					er=os.system(com)
					if(er == 0):
						com = "bzip2 -d %s.trans.bz2" % M
						er=os.system(com)
					else:
						com = "wget http://www.exomol.com/db/%s/%s.trans" % (P, M)
						er=os.system(com)
						if(er != 0):
							print("Error in download .trans file")
							exit()
				l[nu]=int(subprocess.check_output(['wc', '-l', "%s.trans" % M]).split()[0])
		print(nu, l[nu])



	print("download finished")

	if(PrintISO == 1):
		f = open(("%s.param" % M),'w')

		lStates=int(subprocess.check_output(['wc', '-l', "%s.states" % M]).split()[0])

		print("Database = 2", file = f)
		print("Molecule number = %d" % m, file = f)
		print("Name = %s" % M, file = f)
		print("Number of Isotopes = 1", file = f)
		print("#Id Abundance      Q(296K)   g     Molar Mass(g)  partition file :", file = f)
		print("0 1.0             0.0       0      %s        %s.pf" % (mass, M), file = f)
		print("Number of columns in partition File = %s" % npfcol, file = f)
		print("Number of line/transition files = %d" % n, file = f)
		print("Number of lines per file :", file = f)
		for nu in range(n):
			print("%d" % l[nu], file = f)
		print("Line file limits :", file = f)
		for nu in range(n + 1):
			print("%d" % (nu * s), file = f)
		print("#ExoMol :", file = f)
		print("Number of states = %d" % lStates, file = f)	
		print("Number of columns in transition files = %s" % ntcol, file = f)
		print("Default value of Lorentzian half-width for all lines = %s" % dL, file = f)
		print("Default value of temperature exponent for all lines = %s" % dn, file = f)
		print("Version = %s" % version, file = f)

if __name__ == "__main__":
	if(len(sys.argv) < 2):
		print("Error, no molecule specified, run python exomol.py <id>")
		exit()

	m = float(sys.argv[1])
	main(m)
