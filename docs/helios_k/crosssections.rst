Generate k-tables from cross section files
==========================================

.. _ktables:

This section describes how HELIOS-K can be used to calculate k-tables from cross section files.
As an example we are using the HITRAN ozone cross section file ``O3_293.0K-0.0Torr_28901.0-40999.0_118.xsc``

Step 1, Species Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~

First, a species param.dat file must be created. The only relevant entry for this process is the molar mass, which is
used to transform the units of the cross sections. 

For our example, the file (named ``03_1_hit20.param``) looks like::

	Database = 0
	Molecule number = 3
	Name = 03_1_hit20
	Number of Isotopologues = 1
	#ID Abundance      Q(296K)   g     Molar Mass(g)  partition file :
	 31            1   3475.0     1   47.984745 q16.txt
	Number of columns in partition File = 2
	Number of line/transition files = 1
	Number of lines per file :
	1
	Line file limits :
	0
	1
	#ExoMol :
	Number of states = 0
	Number of columns in transition files = 0
	Default value of Lorentzian half-width for all lines = 0
	Default value of temperature exponent for all lines = 0
	Version = 0
	 
Step 2, Download the files and set up the param.dat file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the file from the HITRAN website.


Next, the HELIOS-K param.dat file must be set up.

- Set the name, as the downloaded file, without the ending.
- Set the temperature according to the HITRAN file
- Set the specis name, according to step 1
- Set numin according to the HITRAN file 
- Set numax according to the HITRAN file, note that the upper bound numax is not included in HELIOS-K, so set + dnu here 
- Set dnu according to the HITRAN file 
- Set doResampling = 1
- Set doStoreFullK = -1, to load in the converted cross section file (see step 3)
- Set nbins as needed
- Set units as needed. Note that the HITRAN files are in units of cm^2 / molecule
  When unit = 0, is set, then it will be converted to cm^2 / g

Example (param.dat) ::

	name = O3_293.0K-0.0Torr_28901.0-40999.0_118
	T = 293.0
	P = 1.
	PFile = -
	Species Name = 03_1_hit20
	ciaSystem = -
	pathToData =
	numin = 28901.0
	numax = 41001.0
	dnu = 1.0
	Nnu per bin = 0
	cutMode = 0
	cut = 100.0
	doResampling = 1
	nC = 20
	doTransmission = 0
	nTr = 1000
	dTr = 0.05
	doStoreFullK = -1
	pathToK =
	doStoreSK = 1
	nbins = 20
	binsFile = -
	OutputEdgesFile = -
	kmin = 0.0
	qalphaL = 1.0
	doMean = 0
	Units = 0
	ReplaceFiles = 1
	RLOW = 0
	profile = 1


 
Step 3, convert the format
~~~~~~~~~~~~~~~~~~~~~~~~~~

After downloading the .xsc file, it must be converted to a format similar as the HELIOS-K output format.
This can be done with the provided script ``convertXSC.py`` as ::

	python3 convertXSC.py

This script will read the ``param.dat`` file to extract the file names and all relevant information.
It will create a new file Out_<name>.dat containing two columns, the wave number, and the opacity.
This file can now be used with HELIOS-K.

Step 4, run HELIOS-K
~~~~~~~~~~~~~~~~~~~~

Run HELIOS-K to generate the k-tables file::

	./heliosk

A new file Out_<name>_cbin is created.

Step 5 (optional), check the result
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The produced k-table file, can be reconverted to an opacity file to check the results.
This can be done with the provided script  ``recon_bin.py in tools directory`` as::

	python3 recon_bin.py -M Out_O3_293.0K-0.0Torr_28901.0-40999.0_118



