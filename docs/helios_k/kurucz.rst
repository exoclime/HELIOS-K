Kurucz Database
---------------

Step 1, Download the files
~~~~~~~~~~~~~~~~~~~~~~~~~~

The needed files can be downloaded with the following command:

::

   python3 Kurucz3.py -D 1 -Z <z> -I <i>

where ``<z>``\ is the atomic number and ``<i>``\ the ionization state

- D 1, means to download the file ``gfallwn08oct17.dat`` and all available partition function files.
- D 2, means to download only the partition function files.
- D 0, means no download from the Kurucz website.

- Z i, means to use species with atomic number Z = i.
- Z -1, means to use all species. 

- I j, means to use ionisation state j, where 0 are neutral atoms.
- I -1, means to use ionisation states 0, 1 and 2.

If the source file ``gfallwn08oct17.dat`` should be changed, then the script ``Kurucz3.py`` and ``Kurucz4.py``
 must be modified.

This step writes a file ``gfnew<Z><I>.dat`` containing the line list data.


Step 2, Correct isotope fractions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The file ``gfallwn08oct17.dat`` contains some incorrect isopte fractions for hyperfine splitted lines.
This can be corrected with:

::

   python3 KuruczHyper.py -Z <z> -I <i>


This wrtie a new file ``gfnewCorr<Z><I>.dat`` with the corrected data.


Step 3, calculate natural broadening and create the binary files and ``< species >.param`` files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   python3 Kurucz4.py -Z <z> -I <i>


This step creates the ``< species >.param`` and the `*.bin` files. It also calculates missing natural
broadening coefficients. This can take a while to complete.


Step 4, data path
~~~~~~~~~~~~~~~~~

Include the path of the directory, which contains the obtained binary
files, the ``*.pf`` partition function files and the ``*.param`` file to
the HELIOS-K ``param.dat`` file under ``pathToData``.
