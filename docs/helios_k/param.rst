Download and pre-process the line lists
=======================================

If you have access to the hulk cluster at the University of Bern, then most of
the pre-processed files are already available at
``scratch/sigrimm/EXOMOL/``.
This files can then directly be used for the opacity calculation.


Supported line list databases
-----------------------------

HELIOS-K supports line lists from the HITRAN, HITEMP, ExoMol, Kurucz,
NIST or VALD databases. Before the line lists can be used, they have to
be pre-processed into binary files. This saves in generally memory space
and allows HELIOS-K to read the line lists in a more efficient way.

The next section explains the structure of these files and shows how to generate them.



The ``<species_name>.param`` file
---------------------------------

Each molecular or atomic species, which will be processed by HELIOS-K,
needs a ``<species_name>.param`` file, which contains all necessary
information about the line list database and other parameters.
This files can be generated automatically by the tools in HELIOS-K. However, when using
a different database, they must be created manually.


Example
-------

An example for HITRAN 2016 H2O is given in the ``data/01_hit16.param`` file,
which is shown here and available in the repository. An example for
HITEMP and ExoMol are also given in the repository files
``data/01_HITEMP2010.param`` and ``data/1H2-16O__BT2.param``.

::

   01_hit16.param
   -------------------------------------------------------------------
   Database = 0
   Molecule number = 1
   Name = hit16
   Number of Isotopologues = 7
   #ID Abundance      Q(296K)   g     Molar Mass(g)  partition file :
   11 0.997317        174.58    1     18.010565      q1.txt
   12 0.002000        176.05    1     20.014811      q2.txt
   13 3.718840E-04    1052.14   6     19.01478       q3.txt
   14 3.106930E-04    864.74    6     19.01674       q4.txt
   15 6.230030E-07    875.57    6     21.020985      q5.txt
   16 1.158530E-07    5226.79   36    20.020956      q6.txt
   17 2.419700E-07    1027.80   1     20.022915      q129.txt
   Number of columns in partition File = 2
   Number of line/transition files = 1
   Number of lines per file :
   304225
   Line file limits :
   0
   30000
   #ExoMol :
   Number of states = 0
   Number of columns in transition files = 0
   Default value of Lorentzian half-width for all lines = 0
   Default value of temperature exponent for all lines = 0
   Version = 0
   -------------------------------------------------------------------

Arguments
---------

-  Database: Indicate which database format to use.

  -  0 = Hitran or HITEMP
  -  2 = ExoMol
  - 30 = Kurucz atomic database
  - 31 = NIST  atomic database

-  Molecule number: This is optional and follows an old HELIOS-K
   version.
-  Name: name of the line list files.
-  Number of isotopologues in the line list files.
-  For each isotopologue:

   -  ID: for Hitran and HITEMP the same as the first three digits in
      the ``*.par`` files e.g. 11 for 1H2-18O, or 2A for 18O-13C-17O.
      For ExoMol no meaning.
   -  Abundance, only relevant for more than one isotopologue.
   -  Reference partition function value. For ExoMol no meaning.
   -  Molar mass in g.
   -  Name of the partition function file.

-  Number of columns in the partition function file.
-  Number of line/transition files.
-  Line by line, the number of molecular/atomic lines in each
   line/transition file.
-  Line by line, the wavenumber range of each line/transition file,
   must have one entry more than the number of line/transition files to
   cover all lower and upper range limits.
-  Number of states in ExoMol ``*.states``\ files. For other databases no
   meaning.
-  Number of columns in ExoMol ``*.trans``\ files. For other databases no
   meaning.
-  ExoMol default value of Lorentzian half-width. For other databases no
   meaning.
-  ExoMol default value of temperature exponent. For other databases no
   meaning.
-  Version of the line list.
