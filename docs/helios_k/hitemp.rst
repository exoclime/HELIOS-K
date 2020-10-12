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


References
~~~~~~~~~~

HITEMP gives an example how to cite their work:

::

	How to cite HITEMP

	When using the HITEMP database, please cite the appropriate HITEMP article (e.g., [1]) as well as the original sources of data. To assist the user, each line transition has reference codes that are provided for each parameter and are consistent with those in HITRAN for the same molecule. We hope that you will find HITEMP helpful in your research.
	 
	References

	[1] L. S. Rothman, I. E. Gordon, R. J. Barber, H. Dothe, R. R. Gamache, A. Goldman, V. Perevalov, S. A. Tashkun, J. Tennyson, "HITEMP, the high-temperature molecular spectroscopic database", J. Quant. Spectrosc. Radiat. Transfer 111, 2139-2150 (2010). [link to article] [ADS]

	[2] L. S. Rothman, R. B. Wattson, R. R. Gamache, J. Schroeder, A. McCann, "HITRAN, HAWKS and HITEMP: High-Temperature Molecular Database", Proc. SPIE, Atmospheric Propogation and Remote Sensing IV 2471, 105-111 (1995). [link to article] [ADS]

	[3] R. J. Hargreaves, I. E. Gordon, L. S. Rothman, S. A. Tashkun, V. I. Perevalov, S. N. Yurchenko, J. Tennyson, H. S. P. Müller, "Spectroscopic line parameters of NO, NO2, and N2O for the HITEMP database", J. Quant. Spectrosc. Radiat. Transfer 232, 35-53 (2019). [link to article] [ADS]

	[4] R. J. Hargreaves, I. E. Gordon, M. Rey, A. V. Nikitin, V. G. Tyuterev, R. V. Kochanov, L. S. Rothman, "An accurate, extensive, and practical line list of methane for the HITEMP database", Astrophys. J. Supp. Ser. 247, 55 (2020). [link to article] [ADS]

	[5] G. Li, I. E. Gordon, L. S. Rothman, Y. Tan, S.-M. Hu, S. Kassi, A. Campargue, E. S. Medvedev, "Rovibrational line lists for nine isotopologues of the CO molecule in the X1Σ+ ground electronic state", Astrophys. J. Supp. Ser. 216, 15 (2015). [link to article] [ADS]

