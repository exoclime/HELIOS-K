Compilation
===========

HELIOS-K can be compiled with the provided Makefile by typing
``make SM=xx`` to the terminal, where xx corresponds to the compute
capability. For example use ``make SM=20`` for compute capability of
2.0, or ``make SM=35`` for 3.5. A table with all compute capabilities
can be found `here <https://developer.nvidia.com/cuda-gpus>`_.

On Windows machines
-------------------

If using Cygwin on Windows, then HELIOS-K can be compiled the same way
with ``make SM=xx``. If using the Windows Command Prompt, type
``nmake -f MakefileW SM=xx``. Note, that the Windows c++ compiler ``cl``
must be installed, and the compiler path must be loaded in the shell. If
this is not the case, it can be loaded similar to:
``call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\Common7\Tools\vsvars32.bat"``
. See also issue #1 for details.

Download and pre-process the line lists
=======================================

If you work an the University of Bern on the Hulk cluster, then most of
the pre-processed files are already available at
``scratch/sigrimm/EXOMOL/``.

Supported line list databases
-----------------------------

HELIOS-K supports line lists from the Hitran, HITEMP, ExoMol, Kurucz,
NIST or VALD databases. Before the line lists can be used, they have to
be pre-processed into binary files. This saves in generally memory space
and allows HELIOS-K to read the line lists in a more efficient way.
