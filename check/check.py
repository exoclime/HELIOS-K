'''
This script runs the test cases and compares the output of a
new code version to an existing output file

Author: Simon Grimm
Date: March 2020
'''

import numpy as np
import os
import filecmp
success = 1

er = 0

checkfile = 'check.dat'

f= open(checkfile,"w")

'''
############################################################
#compare Heliok with older Version
for i in range(1,19):
#for i in range(10,16):
	os.chdir('TestH%03d' % i)
	print('TestH%03d' % i, file = f)
	os.system('pwd')
	os.system('../../heliosk -name 1')
	os.system('../../../HELIOS-K-master/heliosk -name 0')
	er = os.system('python3 ../../tools/CompareFiles.py -f0 0 -f1 1')
	er1 = os.WEXITSTATUS(er)
	print("error", er1)
	print("error", er1, file = f)
	if(er1 == 100):
		success = 0
	if(er1 == 0 and i != 13 and i!= 14 and i != 15):
		os.system('rm Out*')
		os.system('rm Info*')
	os.chdir('../')
#-----------------------------------------------------------
'''

'''
############################################################
#compare bin files
binFiles = ['TestH013', 'TestH015']
for files in binFiles:
	os.chdir(files)
	print(files, 'bin', file = f)
	os.system('pwd')
	er = os.system('python3 ../../tools/CompareFiles.py -f0 0_bin -f1 1_bin')
	er1 = os.WEXITSTATUS(er)
	print("error", er1)
	print("error", er1, file = f)
	if(er1 == 100):
		success = 0
	if(er1 == 0):
		os.system('rm Out*')
		os.system('rm Info*')
	os.chdir('../')

#bins files
binFiles = ['TestH014']
for files in binFiles:
	os.chdir(files)
	print(files, 'bin', file = f)
	os.system('pwd')
	er = 0
	er += os.system('python3 ../../tools/CompareFiles.py -f0 0_bin0000 -f1 1_bin0000')
	er += os.system('python3 ../../tools/CompareFiles.py -f0 0_bin0001 -f1 1_bin0001')
	er += os.system('python3 ../../tools/CompareFiles.py -f0 0_bin0002 -f1 1_bin0002')
	er1 = os.WEXITSTATUS(er)
	print("error", er1)
	print("error", er1, file = f)
	if(er1 == 100):
		success = 0
	if(er1 == 0):
		os.system('rm Out*')
		os.system('rm Info*')
	os.chdir('../')
#-----------------------------------------------------------
'''

'''
############################################################
#exomol scripts
os.chdir('TestExomol001')
print('TestExomol001', file = f)
os.system('pwd')
er = os.system('python3 ../../exomol2.py')
er1 = os.WEXITSTATUS(er)
print("error", er1)
print("error", er1, file = f)
if(er1 == 100):
	success = 0

#compare entries for BT2
search = open("Exomol_species.dat")
for line in search:
	if "BT2" in line:
		l1 = line
search = open("Exomol_speciesOld.dat")

for line in search:
	if "BT2" in line:
		l2 = line
if(l1 != l2):
	print(l1)
	print(l2)
	print(l1, file = f)
	print(l2, file = f)
	success = 0
os.chdir('../')
#-----------------------------------------------------------
'''

'''
#-----------------------------------------------------------
os.chdir('TestExomol002')
print('TestExomol002', file = f)
os.system('pwd')
er = os.system('python3 ../../exomol.py -M 23Na-1H__Rivlin')
er = os.system(' ../../prepareExomol -M 23Na-1H__Rivlin')
er1 = os.WEXITSTATUS(er)
print("error", er1)
print("error", er1, file = f)
if(er1 == 100):
	success = 0
er = filecmp.cmp('23Na-1H__Rivlin.param', '../../data/23Na-1H__Rivlin.param')
print("error", er)
print("error", er1, file = f)
if(er == False):
	success = 0
os.system('../../heliosk -tuning 0 -name 1')
os.system('../../heliosk -tuning 0 -name 0 -path ../../data/')
er = os.system('python3 ../../tools/CompareFiles.py -f0 0 -f1 1')
er1 = os.WEXITSTATUS(er)
print("error", er1)
print("error", er1, file = f)
if(er1 == 100):
	success = 0
if(er1 == 0):
	er = os.system('rm 23Na-1H__Rivlin*')
	os.system('rm Out*')
	os.system('rm Info*')
os.chdir('../')
#-----------------------------------------------------------
'''

'''
############################################################
#NIST scripts
os.chdir('TestNIST001')
print('TestNIST001', file = f)
os.system('pwd')

er = os.system('python3 ../../nist_ELevels2.py -Z 26 -I 0')
er = os.system('python3 ../../nist_partition.py -Z 26 -I 0')

er = filecmp.cmp('NIST2600.pf', '../../data/NIST2600.pf')
print("error", er)
print("error", er, file = f)
if(er == False):
	success = 0

er = os.system('python3 ../../nist_Lines3.py -Z 26 -I 0')
er = os.system('python3 ../../nist_Lines2.py -Z 26 -I 0')
er = filecmp.cmp('NIST2600.param', '../../data/NIST2600.param')
print("error", er)
print("error", er, file = f)
if(er == False):
	success = 0

os.system('../../heliosk -tuning 0 -name 1')
os.system('../../heliosk -tuning 0 -name 0 -path ../../data/')
er = os.system('python3 ../../tools/CompareFiles.py -f0 0 -f1 1')
er1 = os.WEXITSTATUS(er)
print("error", er1)
print("error", er1, file = f)
if(er1 == 100):
	success = 0
if(success == 1):
	os.system('rm NIST*')
	os.system('rm Out*')
	os.system('rm Info*')
	os.system('rm geckodriver.log')


os.chdir('../')
#-----------------------------------------------------------
'''

'''	
############################################################
#Kurucz scripts
os.chdir('TestKurucz001')
print('TestKurucz001', file = f)
os.system('pwd')

er = os.system('python3 ../../Kurucz2.py -D 1 -Z 26 -I 0')

er = filecmp.cmp('gfnew2600.param', '../../data/gfnew2600.param')
print("error", er)
print("error", er, file = f)
if(er == False):
	success = 0

os.system('../../heliosk -tuning 0 -name 1')
os.system('../../heliosk -tuning 0 -name 0 -path ../../data/')
er = os.system('python3 ../../tools/CompareFiles.py -f0 0 -f1 1')
er1 = os.WEXITSTATUS(er)
print("error", er1)
print("error", er1, file = f)
if(er1 == 100):
	success = 0
if(success == 1):
	os.system('rm gf*')
	os.system('rm partfn*')
	os.system('rm Out*')
	os.system('rm Info*')



os.chdir('../')
#-----------------------------------------------------------
'''


'''
############################################################
os.chdir('TestHitran001')
print('TestHitran001', file = f)
os.system('pwd')
os.system('cp ../../data/01_hit16.par .')
os.system('cp ../../Hitran_species.dat .')
os.system('../../hitran -M 01 -in hit16')
er = filecmp.cmp('01_hit16.param', '../../data/01_hit16.param')
print("error", er)
print("error", er, file = f)
if(er == False):
	success = 0
if(success == 1):
	os.system('rm Hitran_species.dat')
	os.system('rm 01_hit*')
os.chdir('../')
#-----------------------------------------------------------
'''
print('success', success)
print('success', success, file = f)
f.close()
