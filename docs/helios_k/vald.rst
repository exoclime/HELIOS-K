VALD3 Database
--------------

Step 1, Request the files
~~~~~~~~~~~~~~~~~~~~~~~~~

The download of the line lists can either be done manually through
the website, or in a more automated way, by using navigation scripts. 

Before the VALD data can be acessed, one has to register on the VALD website.

Attention, there is a limit of 100000 transition lines per query. Some species have more lines
than that. Is this case, the wavenumber ranges have to be split into two parts. This is
especially the case for Fe and Cr.


Method A, manual request
^^^^^^^^^^^^^^^^^^^^^^^^

- Vist ``http://vald.astro.uu.se/``
- Register with email address
- Choose Mirror Server
- Click ``Extract Element``
- Click ``Unit selections`` and select:

  -  Energy level units: cm-1
  -  Give wavelenths in medium: vacuum
  -  Wavelength units: cm-1
  -  Van der Waals syntax: Default
  -  Isotopic scaling: Off
  -  click ``Save settings``

- click again ``Extract Element``
- check the units, they must be ``Energy unit: 1/cm - Medium: vacuum - Wavelength unit: 1/cm - VdW syntax: default``
- Enter Starting wavenumber, e.g. 0.001
- Enter Ending wavenumber : e.g. 1000000
- Enter Element: e.g Fe 1, (ionizaion number = 1 means neutral)
- Extraction format : Long format
- Retrieve data via : FTP
- Select ``Include HFS splitting``
- Enter a comment in ``Optional comment for request``, for example the element and ionization number.
- Click ``Submit request``



Method B, using scripts
^^^^^^^^^^^^^^^^^^^^^^^

Also this mehod requires to first visit the web interface manually and to
choose the correct units, the same way as Method A.
After that, the vald_request.py script can navigate the webiste and fill in the necessary webforms. 
For this, a geckodriver needs to be installed. 

Then type:

::

   python3 vald_request.py -Z <z> -I <i> -u <email address>

where ``<z>``\ is the atomic number and ``<i>``\ the ionization state (i = 0 is neutral)
The email address must be the same as the user account from the VALD website.


Step 2, download the files
~~~~~~~~~~~~~~~~~~~~~~~~~~

After a few minutes you will get an email from VALD with a link to download the data file.


Method A, manual download
^^^^^^^^^^^^^^^^^^^^^^^^^
Download the file and unpack it. Rename the file to ``VALD<Z><I>.dat``.
E.g VALD2600.dat for Z = 26, I = 0 (I = 0 means neutral atoms).


Method B, using scripts
^^^^^^^^^^^^^^^^^^^^^^^

The requested files are stored temporarely on the VALD webiste with a name similar to
``SimonGrimm.4199418``, where the first part is the user name, and the second part is the 
job number. The vald_download.py script can download all requested files, scan which element
they contain and rename them according to the element number. 

::

   python3 vald_download.py -jobstart <start> -jobend <end> -jobname <name>

The name is the user name from the VALD webiste, e.g SimonGrimm.
The jobstart is the first jobnumber (check email from VALD) .e.g 4199418.
The jobend is the last jobnumber (check email from VALD) .e.g 4199420.

The files are renamed as ``VALD<Z><I>.dat``.


Step 3, partition function files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Partition functions are not available in the VALD database.
Therefore we use the partition functions calculated from the NIST database.


Step 4, calculate natural broadening coefficients and create the binary files and ``< species >.param`` files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For many species, the vald database provides broadening coefficients,
but not for all. The missing coefficients are calculated in this step.
That can take a while to complete.

::

   python3 vald.py -Z <z> -I <i>


This step creates the ``< species >.param`` and the `*.bin` files. 


Step 5, data path
~~~~~~~~~~~~~~~~~

Include the path of the directory, which contains the obtained binary
files, the ``*.pf`` partition function files and the ``*.param`` file to
the HELIOS-K ``param.dat`` file under ``pathToData``.



References
~~~~~~~~~~

Vald provides reference files for each species. We rename these ``*.bib`` files
the same way as the data files. When using data from VALD, all these references
must be cited in publications.


We list here also the citation conditions from VALD:

::

	A short note on using VALD data :

	The VALD data and services are provided free of charge on the following conditions:

	    1. If VALD has been used in your research work, you agree to include the acknowledge
	       below as well as citations to the following papers where appropriate:

	    "This work has made use of the VALD database, operated at Uppsala University, the
	     Institute of Astronomy RAS in Moscow, and the University of Vienna."

		Ryabchikova T., Piskunov, N., Kurucz, R.L., et al., Physics Scripta, vol 90, issue 5, article id. 054005 (2015), (VALD-3)
		Kupka F., Ryabchikova T.A., Piskunov N.E., Stempels H.C., Weiss W.W., 2000, Baltic Astronomy, vol. 9, 590-594 (2000), (VALD-2)
		Kupka F., Piskunov N.E., Ryabchikova T.A., Stempels H.C., Weiss W.W., A&AS 138, 119-133 (1999), (VALD-2)
		Ryabchikova T.A. Piskunov N.E., Stempels H.C., Kupka F., Weiss W.W. Proc. of the 6th International Colloquium on Atomic Spectra and Oscillator Strengths, Victoria BC, Canada, 1998, Physica Scripta T83, 162-173 (1999), (VALD-2)
		Piskunov N.E., Kupka F., Ryabchikova T.A., Weiss W.W., Jeffery C.S., A&AS 112, 525 (1995) (VALD-1) 

	    2. You also agree to reference to original source(s) of the line data that are important for your research.
	       A complete list of references is available from the bibliography section of the VALD manual. A short list
	       of references is compiled and sent to you with every reply from VALD. 

