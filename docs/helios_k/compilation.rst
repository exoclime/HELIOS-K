Requirements
============

| The opacity calculation from HELIOS-K is running on GPUs and require a Nvidia GPUi
 with compute capability of 3.0 or higher. It can run on Nvidia Tesla, GeForce and Quadro GPUs.
| The code needs the CUDA toolkit to be installed. This can be downloaded from
 https://developer.nvidia.com/cuda-downloads.
| Helper code for downloading the files and preprocessing are written in C++ and Python3.
 They require the following libraries:

- exomol.py and exomol2.py

  - bs4
  - requests
  - sys
  - os
  - subprocess
  - numpy
  - argparse
  - math

- Kurucz2.py

  - numpy
  - math
  - struct
  - os
  - argparse

- nist_ELevels2.py and nist_Lines3.py

  - sys
  - argparse
  - csv
  - requests

- nist_partition.py

  - numpy
  - sys
  - argparse

- nist_Lines2.py

  - numpy
  - struct
  - math
  - csh
  - pandas
  - argparse
  - re
 

Note, when using a computing cluster, all these libraries can be installed locally in the home directory  with: 

::

  pip3 install --user <package name>



Compilation
===========

HELIOS-K can be compiled with the provided Makefile by typing

::

  make SM=xx

into the terminal, where ``xx`` corresponds to the compute
capability of the installed GPU. For example use ``make SM=20`` for compute capability 2.0, or ``make SM=35`` for 3.5. A table with all compute capabilities
can be found `here <https://developer.nvidia.com/cuda-gpus>`_.
On a computing cluster, eventually the CUDA module must be loaded before compiling the code.

On Windows machines
-------------------

If using Cygwin on Windows, then HELIOS-K can be compiled the same way
with ``make SM=xx``. If using the Windows Command Prompt, type
``nmake -f MakefileW SM=xx``. Note, that the Windows C++ compiler ``cl``
must be installed, and the compiler path must be loaded in the shell. If
this is not the case, it can be loaded similar to this command:
``call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\Common7\Tools\vsvars32.bat"``
, where the exact path and file name must eventually be changed.

