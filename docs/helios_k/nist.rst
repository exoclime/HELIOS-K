NIST Database
-------------

Step 1 Download the Energy levels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step can be done either manually or by using a script.

-  For manuall download:

   -  visit
      ``https://physics.nist.gov/PhysRefData/ASD/levels_form.html``
   -  enter the species and ionization state, e.g. Z= 3 0
   -  select the Format output to Tab-delimited
   -  select the ``g`` option
   -  store the data in a file called ``test.dat``
   -  remove all ``"`` in the file

-  Automatical download:

   -  type ``python3 nist_ELevels.py -Z <z> -I <i>``, where ``<z>``\ is
      the atomic number and ``<i>``\ the ionization state
   -  The script will download the data and store it to ``test.dat``.
   -  This script will need a geckodriver to be installed .

Step 2 Generate partition functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run ``python3 nist_partition.py - -Z <z> -I <i>``. This script will read
the previously produces file ``test.dat``.

Step3 Download atomic masses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  visit
   ``https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses``
-  select ``All Elements``
-  select ``Linearized ASCII Output``
-  select ``All isotopes``
-  click ``Get Data``\ and store the data in a file called
   ``masses.txt``

Step 4 Download the line list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step can be done either manually or by using a script.

-  For manuall download:

   -  `https://physics.nist.gov/PhysRefData/ASD/lines_form.html <https://physics.nist.gov/PhysRefData/ASD/lines_form.html>`__
   -  enter the species and ionization state, e.g. Z= 3 0
   -  select ``Show Advanced Settings``
   -  select the Format output to csv
   -  Select ``Wavenumbers (in cm-1)``
   -  Select ``Wavenumbers (all wavelenghts)``
   -  Select ``g``
   -  click ``Retrieve Data``\ and store the data in a file called
      ``test.dat``.
   -  remove all ``?, =, [,],( and )``\ from the file

-  Automatical download:

   -  type ``python3 nist_Lines.py -Z <z> -I <i>``.
   -  This script will need a geckodriver to be installed .

.. _step-5-create-<-species->.param-file-and-binary-files:

Step 5 Create ``< species >.param`` file and binary files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All necessary files can be created with:

::


   python3 nist_Lines2.py -Z <z> -I <i>

Step 6 data path
~~~~~~~~~~~~~~~~

Include the path of the directory, which contains the obtained binary
files, the ``.pf`` partition function files and the ``.param`` file to
the HELIOS-K ``param.dat`` file under ``pathToData``.