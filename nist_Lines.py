#https://towardsdatascience.com/controlling-the-web-with-python-6fceb22c5f08



'''
Download the proper geckodriver depending upon your operating system and System Architecture from the url here - geckodriver

Now follow the below steps -
Extract the geckodriver file from the downloaded zip file. Now depending upon your operating system do the following.

For Linux system :

Open terminal and login as root user. copy/move the extracted geckodriver to bin direcctory.
In my case I moved the file to /usr/bin directory. Because the driver finds geckodriver binary in '/usr/bin' path and the problem is solved now.

To move the file inside bin directory use command like -

$mv 'geckodriver binary source path' 'destination path'
destination path should be the binary folder path as per user system.

'''

#for pyperclip:
#sudo apt install xsel

import sys
import time
import pyperclip
import subprocess
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import argparse

def Lines(Z, I):

	print(Z, I)

	driver = webdriver.Firefox()
	driver.get('https://physics.nist.gov/PhysRefData/ASD/lines_form.html')

	spectra = driver.find_element_by_name("spectra")
	spectra.send_keys("Z = %d %d" % (Z, I))
	#print("click advanced")
	advanced = driver.find_element_by_name("show_advanced")
	advanced.click()


	fo = driver.find_element_by_id("format")
	fo.send_keys(Keys.DOWN, Keys.DOWN)


	wn = driver.find_element_by_id("show_wn")
	wn.click()
	#print("click wn")

	wl = driver.find_element_by_name("show_av")
	wl.send_keys(Keys.DOWN, Keys.DOWN, Keys.DOWN)


	g = driver.find_element_by_id("g_out")
	g.click()

	submit = driver.find_element_by_name("submit")
	submit.click()

	lenS = 0
	lenSOld = 0
	for t in range(10):
		lenSOld = lenS
		data = driver.find_element_by_css_selector("body")
		data.send_keys(Keys.CONTROL, "a")
		data.send_keys(Keys.CONTROL, "c")


		s = pyperclip.paste()
		lenS = len(s)
		print("lenght of data", len(s))
		if(lenS == lenSOld and lenS > 3000):
			break
		time.sleep(5)
	driver.quit()


	s1 = s.replace('?', '')
	s2 = s1.replace('=', '')
	s3 = s2.replace('[', '')
	s4 = s3.replace(']', '')
	s5 = s4.replace('(', '')
	s6 = s5.replace(')', '')
	with open("test.dat", "w") as f:
	    f.write(s6)



	'''
	echo "T"


	#replace " with ' '
	#sed -i 's/"//g' test.dat


	'''

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


