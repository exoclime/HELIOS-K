#https://towardsdatascience.com/controlling-the-web-with-python-6fceb22c5f08



'''
from:
https://python-forum.io/Thread-Getting-error-geckodriver-executable-needs-to-be-in-PATH

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

def nist(Z, I):

	print(Z, I)

	driver = webdriver.Firefox()
	driver.get('https://physics.nist.gov/PhysRefData/ASD/levels_form.html')



	spectrum = driver.find_element_by_name("spectrum")
	spectrum.send_keys("Z = %d %d" % (Z, I))

	#Format output, set to Tab-delimited
	fo = driver.find_element_by_name("format")
	fo.send_keys(Keys.DOWN, Keys.DOWN, Keys.DOWN, Keys.DOWN)



	#go to g checkbox
	g = driver.find_element_by_name("g_out")
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



	s1 = s.replace('"', '')
	with open("test.dat", "w") as f:
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

