NIST Database
-------------

Step 1, Download the Energy levels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step can be done either manually or by using a script.

-  For manuall download:

   -  visit
      ``https://physics.nist.gov/PhysRefData/ASD/levels_form.html``
   -  enter the species and ionization state, e.g. ``Z= 3 0``
   -  select the Format output to Tab-delimited
   -  select the ``g`` option
   -  store the data in a file called ``NIST_ELevels<Z><I>.dat``
   -  remove all ``"`` in the file

-  Automatical download:

   -  type ``python3 nist_ELevels2.py -Z <z> -I <i>``, where ``<z>``\ is
      the atomic number and ``<i>``\ the ionization state , e.g. -Z 3 -I 0
   -  The script will download the data and store it in a file called ``NIST_ELevels<Z><I>.dat``.

Step 2, Generate the partition functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run ``python3 nist_partition.py - -Z <z> -I <i>``. This script will read
the previously produces file ``NIST_ELevels<Z><I>.dat``, and writes a ``*.pf`` file.


Step 3, Download atomic masses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  visit
   ``https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses``
-  select ``All Elements``
-  select ``Linearized ASCII Output``
-  select ``All isotopes``
-  click ``Get Data``\ and store the data in a file called
   ``masses.txt``

Step 4, Download the line list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step can be done either manually or by using a script.

-  For manuall download:

   -  visit `https://physics.nist.gov/PhysRefData/ASD/lines_form.html <https://physics.nist.gov/PhysRefData/ASD/lines_form.html>`__
   -  enter the species and ionization state, e.g. ``Z= 3 0``
   -  select ``Show Advanced Settings``
   -  select the Format output to csv
   -  Select ``Wavenumbers (in cm-1)`` (keep ``Observed``, ``Ritz`` and ``Uncertainties`` selected)
   -  Select ``Vacuum (all wavelengths)``
   -  Select ``g``
   -  click ``Retrieve Data``\ and store the data in a file called
      ``NIST_Lines<Z><I>.dat``.
   -  remove all ``?, =, [,],( and )``\ from the file

-  Automatical download:

   -  type ``python3 nist_Lines3.py -Z <z> -I <i>``.
   -  The script will download the data and store it in a file called ``NIST_Lines<Z><I>.dat``.

.. _step-5-create-<-species->.param-file-and-binary-files:

Step 5, Create ``<species_name>.param`` file and the binary files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All necessary files can be created with:

::

   python3 nist_Lines2.py -Z <z> -I <i>

This step includes the calculation of the natural broadening coefficients. This can take o moment to complete.

Step 6, data path
~~~~~~~~~~~~~~~~~

Include the path of the directory, which contains the obtained binary
files, the ``*.pf`` partition function files and the ``*.param`` file to
the HELIOS-K ``param.dat`` file under ``pathToData``.
