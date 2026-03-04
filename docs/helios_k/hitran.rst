.. _hitran:

HITRAN
------

.. _HITRAN_step_1:

Step 1, Species Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~

The HELIOS-K repository provides a file called ``Hitran_species.dat``.
This file contains all available species of the Hitran database.
The file can be recreated or updated with the python code::

	 python3 hitran2.py



The file format is::

	Species ID, Abundance, Molecular Mass in g/mol, Q(296K), partition function file, gi, Isotopologue Formula


The same information can be found at
``https://hitran.org/docs/iso-meta/``.

The species ID consists of a two digits molecule ID and a one digit
local isotopologue ID. Note that the local isotogologue ID can sometimes
consist of non numerical values. e.g A or B,

For identifying a species, the molecule number and the isotopologue
number should be combined, e.g. ``01_1`` for 1H2-16O or ``01`` for all
isotopogolues from H20.


.. _HITRAN_step_2:

Step 2, Download the files
~~~~~~~~~~~~~~~~~~~~~~~~~~

The line list files must be downloaded manually from ``www.hitran.org``,
and note that it is necessary to register on the hitran homepage. To
download the files, select ``DataAcess`` and then ``Line-by-line``. Select the
``molecule id``, select ``all isotopologues`` (single isotopologues can be
filtered later), leave wavenumber range blank, select ``.par``\ file and
store the file on your computer under the name ``Molecule-ID_hit16.par``.
It is important to use this naming convention, otherwise the following steps will not work properly. 

Example for water (molecule id = 01) and a choose line list name = `hit16` (for Hitran 2016)::

	01_hit16.par

Download also all the necessary partition function files from
www.hitran.org under Documentation/Isotogologues.

.. _step-3-create-<-species->.param-file-and-binary-files:

Step 3, create ``<species_name>.param`` file and the binary files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All necessary files can be created with:

::

   ./hitran -M < molecule ID > -ISO < isotopologue ID > -in < line list name >

The ``<molecule ID>`` is the two digit molecule number, e.g. ``01`` for
H2O. The ``<isotopologue ID>`` is the hitran internal isotopologue
identifier, e.g. ``1`` for 1H2-16O.

The ``<line list name>`` is the name that was given in the download section, e.g. ``hit16``.

Example for water (molecule id = 01), isotopolgue = 1  and a choose line list name = `hit16`::

	./hitran -M 1 -ISO 1 -in hit16


This example will produce the species parameter file ``01_1_hit16.param`` and the binary line
list file called ``01_1_hit16.bin``.  

.. _step-4-data-path-1:

Step 4, data path
~~~~~~~~~~~~~~~~~

Include the path of the directory, which contains the obtained binary
files, the ``*.txt`` partition function files and the ``*.param`` file to
the HELIOS-K ``param.dat`` file under ``pathToData``.
