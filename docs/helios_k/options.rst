

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



