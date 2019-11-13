Kurucz Database
---------------

Step 1 Download the files and create binary files and .param files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The needed files can be downloaded with the following comand:

::

   python3 Kurucz2.py -D 1 -Z <z> -I <i>

where ``<z>``\ is the atomic number and ``<i>``\ the ionization state. z
= -1 means to download all atoms, i = -1 means to download Z = 0, 1 and
2.

The script will downlaod the file gfallwn08oct17.dat and all available
partition function files. If the source file should be changed, the the
script ``Kurucz2.py`` must be modified.

Step 2 data path
~~~~~~~~~~~~~~~~

Include the path of the directory, which contains the obtained binary
files, the ``.pf`` partition function files and the ``.param`` file to
the HELIOS-K ``param.dat`` file under ``pathToData``.
