Kurucz Database
---------------

Step 1, Download the files and create the binary files and ``< species >.param`` files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The needed files can be downloaded with the following command:

::

   python3 Kurucz2.py -D 1 -Z <z> -I <i>

where ``<z>``\ is the atomic number and ``<i>``\ the ionization state. z
= -1 means to download all atoms, i = -1 means to download Z = 0, 1 and
2.

- D 1, means to download the file ``gfallwn08oct17.dat`` and all available partition function files.
- D 2, means to download only the partition function files.
- D 0, means no download from the Kurucz website.

If the source file ``gfallwn08oct17.dat`` should be changed, then the script `Kurucz2.py` must be modified.
This script will also create the ``< species >.param`` and the `*.bin` files.


Step 2, data path
~~~~~~~~~~~~~~~~~

Include the path of the directory, which contains the obtained binary
files, the ``*.pf`` partition function files and the ``*.param`` file to
the HELIOS-K ``param.dat`` file under ``pathToData``.
