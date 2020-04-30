
The binsfile options
====================

When a ``binsFile`` name is given in the ``param.dat`` file, then this
file is used to generate the edges of the bins, which can be irregular.
Note that this option does not support the doResampling and
doTransmission options. The binsfile must contain line by line the edges
of the bins in cm^-1.

For example:

::

   0.5
   50.0
   100.0
   200.0

The output edges option
=======================

When a ``outputedgesFile`` name is given in the ``pram.dat`` file, then
this file is used to specify the averaged output positions of the
``Out_<name>_bin.dat`` files. The file must contain line by line the
positions in y.

For example:

::

   0.0
   0.1
   0.45
   0.78
   1.0


The P file option
=================

When a ``PFile`` name is given in the ``param.dat`` file, then this file
is used to read multiple values for P. This option is useful to speed up
the performance, because multiple reads from the data files can be
avoided. Too many entries in the Pfile can lead to a memory shortage.

For example:

::

   1.0
   10.0
   100.0

The Species file option
=======================

This option must be used to calculate opacities for gas mixtures,
containing multiple species. The File contains in the two columns the
species name, and the number fraction.

For example:

::

   01_hit16    0.9
   05_hit16    0.1

This example will produce an opacitiy with 90% H2O and 10% CO.


