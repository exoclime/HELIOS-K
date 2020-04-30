HITRAN
------

.. _HITRAN_step_1:

Step 1, Species Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~

The HELIOS-K repository provides a file called ``Hitran_species.dat``.
This file contains all available species of the Hitran database. The
file format is:

``Species ID, Abundance, Molecular Mass, Q(296K), partition function file, Isotopologue Formula``
The same information can be found at
``https://hitran.org/docs/iso-meta/``.

The species ID consitst of a two digits molecule ID and a one digit
local isotopologue ID. Note that the local isotogologue ID can sometimes
consist of non numerical values. e.g A or B,

For identifying a species, the molecule number and the isotopologue
number should be combinded, e.g. ``01_1`` for 1H2-16O or ``01`` for all
isotopogolues from H20.

The ``Hitran_species.dat`` file can be recreated or updated with the
python code ``hitran2.py``.

.. _HITRAN_step_2:

Step 2, Download the files
~~~~~~~~~~~~~~~~~~~~~~~~~~

The line list files must be downloaded manually from `www.hitran.org``,
and note that it is necessary to register on the hitran homepage. To
download the files, select DataAcess and then Line-by-line. Select the
molecule id, select all isotopologues (single isotopologues can be
filtered later), leave wavenumber range blank, select ``.par``\ file and
store the file on your computer under the name ``Molecule-ID_hit16.par``,
e.g. 01_hit16.par for H2O.

Download all the necessary partition function files from
Documentation/Isotogologues.

.. _step-3-create-<-species->.param-file-and-binary-files:

Step 3, create ``<species_name>.param`` file and the binary files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All necessary files can be created with:

::

   ./hitran -M < molecule ID > -ISO < isotopologue ID > -in < line list name >

The ``<molecule ID>`` is the two digit molecule number, e.g. ``01`` for
H2O. The ``<isotopologue ID>`` is the hitran internal isotopologue
identifier, e.g. ``1`` for 1H2-16O. The ``<line list name>`` is the name
that was given in the download section, e.g. ``hit16``.

.. _step-4-data-path-1:

Step 4, data path
~~~~~~~~~~~~~~~~~

Include the path of the directory, which contains the obtained binary
files, the ``*.txt`` partition function files and the ``*.param`` file to
the HELIOS-K ``param.dat`` file under ``pathToData``.
