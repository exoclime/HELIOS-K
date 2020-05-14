HITEMP
------

.. _step-1-species-properties-1:

Step 1, Species Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~

For HITEMP, the same ``Hitran_species.dat`` is needed as for HITRAN. See
:ref:`HITRAN_step_1`.

.. _step-2-download-the-files-1:

Step 2, Download the files
~~~~~~~~~~~~~~~~~~~~~~~~~~

The line list files must be downloaded manually from
``https://hitran.org/hitemp/``. The same partition function files are
needed as for HITRAN, see
:ref:`HITRAN_step_2`.

.. _step-3-create-<-species->.param-file-and-binary-files-1:

Step 3, create ``< species >.param`` file and the binary files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All necessary files can be created with:

::

   ./hitran -M < molecule ID > -ISO < isotopologue ID > -in < line list name >

The ``<molecule ID>`` is the two digit molecule number, e.g. ``01`` for
H2O. The ``<isotopologue ID>`` is the hitran internal isotopologue
identifier, e.g. ``1`` for 1H2-16O. The ``<line list name>`` is the name
that was given in the download section, e.g. ``HITEMP2010``.

.. _step-4-data-path-2:

Step 4, data path
~~~~~~~~~~~~~~~~~

Include the path of the directory, which contains the obtained binary
files, the ``*.txt`` partition function files and the ``*.param`` file to
the HELIOS-K ``param.dat`` file under ``pathToData``.

