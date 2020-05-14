Using HELIOS-K
==============

Before HELIOS-K can be used, the molecular or atomic line-lists must be
downloaded and pre-processed. HELIOS-K provides some scripts for that,
which are described in the previous section.

All necessary parameters for HELIOS-K are set in the file ``param.dat``,
which is described next. Some parameters can also be set by
console arguments.

HELIOS-K can be started with

::

   ./heliosk

followed by optional console arguments, which are listed later.


HELIOS-K Input parameters
-------------------------

The input parameters can be specified in the ``param.dat`` file. The
used parameters are listed here, the order can not be changed.

-  name: This name will appear in the output filenames.
-  T: Temperature in Kelvin
-  P: Pressure in Atmospheres
-  PFile: A '-' ignores this option, otherwise this option specifies a
   filename which contains multiple values for P.
-  Species Name: The name of the molecular or atomic param file. e.g.
   01_hit16 or 1H2-16O__BT2.
-  SpeciesFile: A '-' ignores this option, otherwise this option
   specifies a filename which contains multiple species for gas mixtures.
-  ciaSystem: A '-' ignores this option, otherwise a cia system is read
   here. supported are H2-H2, H2-H2_eq, H2_H2_norm, H2-He, H2-He_eq,
   H2-He_norm, H2-CH4_eq, H2-CH4_norm and H2-H.
-  pathToData: The location where the HITRAN or HITEMP data files are
   located, e.g. pathToData = ../HITEMP/ , pathToData = /data/HITEMP/ or
   empty when the files are in the same directory pathToData =
-  numin: minimum wavenumber in 1/cm
-  numax: maximum wavenumber in 1/cm
-  dnu: spatial resolution in wavenumber in 1/cm
-  Nnu per bin: The number of points in wavenumbers nu per bin.
 
   -  0; Nnu is ignored and dnu is used
   -  otherwise, dnu is overwritten.

-  cutMode: Set the voigt profile cutting mode:

   -  0: cut at absolute wavenumbers
   -  1: cut at multiple values of Lorentz half widths
   -  2: cut at multiple values of Doppler half widths

-  cut: value where to cut the line wings, according to the cutMode,
   if cut = 0, then no cutting is performed.
-  doResampling:
 
   - 0: no resampling is done.
   - 1: all the sorted opacity functions per bin are resampled with a Chebyshev
     polynomial, and the coefficients from all bins are written to the file ``Out_<name>_cbin.dat``.
   - 2: all the sorted opacity functions per bin are resampled with a Chebyshev

 polynomial, and the coefficients from each bin are written to a sepparate file
 ``Out_<name>_cbin<bin_number>.dat``.


-  nC: Number of Chebyshev coefficients use for the resampling.
-  doTransmission: Compute the transmission function tau for a set of column masses m.

   - 0: no transmission function is calculated.
   - 1: the transmission function per bin is calculated,
     and all bins are written to the file ``Out_<name>_tr.dat``.
   - 2: the transmission function per bin is calculated,
     and each bin is written to a sepparate file ``Out_<name>_tr<bin_number>.dat``.

-  nTr: nummber of points used for the transmission function.
-  dTr: spacing of the points used for the transmission function in exp
   space: m_i = exp((i - nTr/2) \* dTr)
-  doStoreFullK: 

   - 0: the opacity functions is not written to a file.
   - 1: the full unsorted opacity function is written to the text file ``Out_<name>.dat``.
   - 2: the full unsorted opacity function is written to the binary file ``Out_<name>.bin``.
   - -1: the opacity functions is not calculated, but read in from a text file.
     The name of the file is specified in the next argument.
   - -2: the opacity function is nor calculated, but read in from a binary file.
     The name of the file is specified in the next argument.

- pathToK: path and filename to read in an opacity function. When left blank, then no file is read in.

-  doStoreSK:

   - 0: the resampled opacity function is not written to a file.
   - 1: all bins from the resampled opacity function are written to the same
     file ``Out_<name>_bin.dat``.
   - 2: all bins from the resampled opacity function are written to a sepparate
     file   with names ``Out_<name>_bin<bin_number>.dat``.

-  nbins: number of bins
-  binsfile: A '-' ignores this option, otherwise this option specifies
   a filename which contains the edges of the bins, which can be
   irregular. This option overrides the numin, numax and nbins
   arguments.
-  OutputEdgesFile: A '-' ignores this option, otherwise this option
   specifies a file name which contains the edges of the output
   locations in y for each bin.
-  kmin: minimal value for the opacity function
-  qalphaL: q value in the Lorentzian half width q = Pself / P
-  gammaF: scaling factor for the Lorentzian half width gamma factor,
   gamma = gamma \* gammaF
-  doMean: 

   - 0: nothing is done here.
   - 1: Calculate the Planck and Rosseland opacity means.

-  Units: The units of the opacities.

   - 0: cm^2 / g,
   - 1: cm^2 / molecule

-  ReplaceFile:

   - 1: all existing output files are overwritten.
   - 0: the data is appended to the existing files.

-  profile:

   - 1: Voigt
   - 2: Lorentzian
   - 3: Gaussia
   - 4: cross section

-  doTuning:

   - 0: nothing is done here
   - 1: use self tuning routines to specify the best kernel parameters.


Console Arguments
-----------------

Instead of using the parameter file, some arguments can also be set
by console arguments. The console arguments have the highest priority
and are overwriting the arguments of the ``param.dat`` file. The options
are:

-  -name ``<c>``: name
-  -T ``<double>`` : T
-  -P ``<double>`` : P
-  -M ``<int>`` : Molecule Name
-  -path ``<c>`` : pathToData
-  -pathK ``<c>`` : pathToK
-  -numin ``<double>`` : numin
-  -numax ``<double>`` : numax
-  -dnu ``<double>`` : dnu
-  -cutM ``<int>`` : cutMode
-  -cut ``<double>`` : cut
-  -dR ``<int>`` : doResampling
-  -nC ``<int>`` : nC
-  -dT ``<int>`` : doTRansmission
-  -nTr ``<int>`` : nTr
-  -dTr ``<double>`` : dTr
-  -dSF ``<int>`` : doStoreFullK
-  -dSS ``<int>`` : doStoreSK
-  -nbins ``<int>`` : nbins
-  -kmin ``<double>`` : kmin
-  -dev ``<int>`` : Device number (For multiple GPU systems)
-  -q ``<double>`` : qalphaL
-  -gammaF ``<double>`` : gammaF
-  -Mean ``<int>`` : doMean
-  -tuning ``<int>`` : doTuning

where ``<c>``\ is a string, ``<double>`` a floating point number, and
``<int>``\ an integer.


Code parameters
---------------

The file define.h contains the physical parameters and some code
parameters. After changing some entry here, the code needs to be recompiled.
The code parameters are:

-  def_TOL: Tolerance parameter in the Voigt function. See Algorithm 916
-  def_nthmax: Maximum number of threads per kernel launch. In 2.0 GPUs
   it can not be larger than 32768.
-  def_nlmax: Maximum number of molecular lines per kernel launch.
   Setting a lower number prevents from a time-out on Desktop machines.
-  def_maxlines: Maximum number of lines stored on the GPU.
-  def_maxfiles: Maximum number of files per molecule.
-  def_NmaxSample: Maximum Number of resample coefficients for K(y)

When using a Desktop GPU running an x session, the runtime of a single
kernel launch can be limited to a few seconds. Choosing smaller values
for nlmax and nthmax splits the kernel into smaller parts. But it makes
the code a bit slower.



