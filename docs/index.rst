HELIOS-K
========

HELIOS-K is an opacity calculator running on GPUs, by Simon Grimm & Kevin Heng.
University of Bern.

Introduction
~~~~~~~~~~~~

HELIOS-K calculates opacity functions for planetary atmopheres by using opacity line lists from different databases. Before the opacity functions can be calculated, the line lists need to be downloaded and preprocessed into files that can be read from HELIOS-K. HELIOS-K provides tools to automatically download and preprocess files from the ExoMol, HITRAN, HITEMP, NIST and Kurucz databases. 
HELIOS-K is running on GPUs and require a Nvidia GPU with compute capability of 2.0 or higher.

Setup
~~~~~

.. toctree::
  :maxdepth: 1

  helios_k/compilation.rst
  helios_k/param.rst

Databases
~~~~~~~~~

.. toctree::
  :maxdepth: 1

  helios_k/exomol.rst
  helios_k/hitran.rst
  helios_k/hitemp.rst
  helios_k/kurucz.rst
  helios_k/nist.rst

Running HELIOS-K
~~~~~~~~~~~~~~~~

.. toctree::
  :maxdepth: 1

  helios_k/heliosk.rst

Citations
~~~~~~~~~

If you make use of HELIOS-K in your work, please cite
`Grimm & Heng (2015) <https://ui.adsabs.harvard.edu/abs/2015ApJ...808..182G/abstract>`_:

.. code-block::

    @ARTICLE{2015ApJ...808..182G,
           author = {{Grimm}, Simon L. and {Heng}, Kevin},
            title = "{HELIOS-K: An Ultrafast, Open-source Opacity Calculator for Radiative Transfer}",
          journal = {\apj},
         keywords = {methods: numerical, planets and satellites: atmospheres, radiative transfer, Astrophysics - Earth and Planetary Astrophysics, Physics - Atmospheric and Oceanic Physics},
             year = "2015",
            month = "Aug",
           volume = {808},
           number = {2},
              eid = {182},
            pages = {182},
              doi = {10.1088/0004-637X/808/2/182},
    archivePrefix = {arXiv},
           eprint = {1503.03806},
     primaryClass = {astro-ph.EP},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2015ApJ...808..182G},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

