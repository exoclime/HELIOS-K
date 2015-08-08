# HELIOS-K #
#### Authors: Simon Grimm, Kevin Heng ####

# Limitations #
The current version of the code supports only the molecules H2O, CO2, CO and CH4. 
When HITEMP is used, only H20, C02 and C0 are supported.
All the other molecules will be included soon. 


# Requirements #

HELIOS-K runs on Nvidia GPUs with compute capability of 2.0 or higher. To be able
to use the Code one has to install the Cuda Toolkit. This can be downloaded
from https://developer.nvidia.com/cuda-downloads.


# Compilation #
HELIOS-K can be compiled with the Makefile by typing 'make SM=xx' to the
terminal, where xx corresponds to the compute capability. For example use 'make SM=20'
for compute capability of 2.0, or 'make SM=35' for 3.5.

# Starting HELIOS-K #
HELIOS-K can be startet with

```
./heliosk
```
followed by optional arguments listed below.


# Input parameters #
The imput parameters can be specified in the 'param.dat' file. The used
parameters are listed here, the order can not be changed.

 * name: The name is included in th Output filenames
 * T: Temperature in Kelvin
 * P: Pressure in Atmospheres
 * PFile: A '-' ignores this option, otherwise this option specifies a file name which contains multiple values for P
 * useHITEMP: when set to 0, the HITRAN files are read. When set to 1 the HITEMP files are read
 * Molecule: Molecule identity according to HITRAN, 1 = H20, 2 = CO2, ...
 * pathToData: The location where the HITRAN or HITEMP data files are located, e.g. pathToData = ../HITEMP/ , pathToData = /data/HITEMP/ or empty when the files are in the same directory  pathToData = 
 * numin: minimum wavenumber in 1/cm
 * numax: maximum wavenumber in 1/cm
 * dnu: spatial resolution in wavenumber in 1/cm
 * cutMode: Set the voigt profile cutting mode:
   * 0: cut at absolute wavenumbers
   * 1: cut at multiple values of Lorentz widths
   * 2: cut at multiple values of Doppler widths
 * cut: value where to cut, according to the cutMode, if cut = 0, then no cutting is performed
 * doResampling: When set to one, then all the sorted opacity functions per bin are resampled with a Chebyshev polynomial, and the coefficients are written to the file 'Out_< name >_cbin.dat'.
 *nC: Number of Chebyshev coefficients use for the resampling
 * doTransmission: When set to one, then the transmission function per bin is computed for various numbers of m, and written to the file 'Out_< name >_tr.dat'.
 * nTr: nummber of points used for the transmission function
 * dTr: spacing of the points used for the transmission function in exp space: m_i = exp((i - nTr/2) * dTr)
 * doStoreFullK: When set to one, then the full unsorted opacity function is written to the file 'Out_< name >.dat'.
 * doStoreSK: When set to one, then the perbin sorted opacity function is written to the file 'Out_< name >_bin.dat'.
 * nbins: number of bins
 * binsfile: A '-' ignores this option, otherwise this option specifies a file name which contains the edges of the bins, which can be irregular. This option overrides the numin, numax and nbins arguments.
 * OutputEdgesFile: A '-' ignores this option, otherwise this option specifies a file name which contains the edges of the output locations in y for each bin.
 * kmin: minimal value for the opacity function 
 * qalphaL: q value in the Lorentz half width q = Pself / P 
 * doMean: Calculate the Planck and Rosseland opacity means
 * Units: The units of the opacities. 0: cm^2 / g, 1: cm^2 / molecule

# Console Arguments #
Instead of using the parameter file, some arguments can also be passed as console arguments. The console arguments have the highest priority and are overwriting the arguments of the param.dat file. The options are:

 * \-name < c >: name
 * \-T < double > : T
 * \-P < double > : P
 * \-HITEMP < int > : useHITEMP
 * \-M < int > : Molecule
 * \-path < c > : pathToData
 * \-numin < double > : numin
 * \-numax < double > : numax
 * \-dnu < double > : dnu
 * \-cutM < int > : cutMode
 * \-cut < double > : cut
 * \-dR < int > : doResampling
 * \-nC < int > : nC
 * \-dT < int > : doTRansmission
 * \-nTr < int > : nTr
 * \-dTr < double > : dTr
 * \-dSF < int > : doStoreFullK
 * \-dSS < int > : doStoreSK
 * \-nbins < int > : nbins
 * \-kmin < double > : kmin
 * \-dev < int > : Device number (For multiple GPU systems) 	
 * \-q < double > : qalphaL
 * \-Mean < int > : doMean

# Code parameters #
The file define.h contains the physical parameters and some code parameters.
After changing some entry, the Code needs to be recompiled.
The code parameters are:

 * TOL: 	Tolerance parameter in the Voigt function. See Algorithm 916
 * NCheb:       Number of Chebychev coefficients in the q.dat file
 * nthmax:	Maximum number of threads per kernel launch. In 2.0 GPUs it can not be larger than 32768. 
 * nlmax:	Maximum number of molecular lines per kernel launch. Setting a lower number prevents from a time out on Desktop machines
 * maxbins:	Maximum number of bins
 * PROFILE:	1 = Voigt, 2 = Lorentz, 3 = Gauss
 * NmaxSample:	Maximum Number of resample coefficients for K(y)

When using a Desktop GPU running an x session, the runtime of a single kernel
launch can be limited to a few seconds. Chosing smaller values for nlmax and nthmax 
splitts the kernel into smaller parts. But it makes the Code a bit slower.


# Input Files #

The following input files must be provided:

 * Line list files from HITRAN or HITEMP in ASCII format. For example for H20 the "01_hit12.par" file from HITRAN, or all the files "01_00000-00050_HITEMP2010.par" to "01_11000-30000_HITEMP2010.par" from HITEMP. These must be stored at the location, specified with the "pathToData" argument. 
 * "q.dat" file containing line by line for each Isotopologue the Chebychev  coefficinets
  for the partition function. The number of coefficients must correspond to the value "NCheb" in the define.h file
 * The file "ISO.h" must contain the specifications for each Molecule:
   * m.NL: The number of Lines in the line list file
   * m.nISO: The number of Isotopologues per Molecule. It must correspond the HITRNN data files. For Example CO2 is limited to 10, not 11 as the web version is.
   * For each Isotopologue: id, AFGL id, Abundance, Q(296K), gj, Molar Mass(g).
   These values can be found in the "molparam.txt" file from HITRAN. id corresponds to the HITRAN Identity:
the  first 2 characters are the Molecule id, the third character is the order along the abundances.
   The AFGL id is the Molecule code. For example 161: 1H16O1H.
   * The name of the corresponding line list file.


# The binsfile options #
When a 'binsFile' name is given in the param.dat file, then this file is used to generate the edges of the bins, which can be irregular. Note that this option
does no support the doResampling and doTransmission options. 
The binsfile must contain line by line the edges of the bins in cm^-1.   
  For example:  
  0.5  
  50  
  100  
  200  

# The output edges option #
When a 'outputedgesFile' name is given in the pram.dat file, then this file is used to specify the output positions of the Out_<name>_bin.dat files.
The file must contain line by line the positions in y.  
  For example:  
  0.0  
  0.1  
  0.45  
  0.78  
  1.0  

# The P file option #
When a 'PFile' name is given in the pram.dat file, then this file is used to read multiple values for P. This option is useful to speed up the performance, 
because multiple reads from the datafiles can be avoided. To many entries in the Pfile can lead to a memory shortage.
  For example:  
  1.0  
  10.0  
  100.0  

# Output Files #
Different Output files are written, depending to the set values in the 'param.dat' file

# Info_< name >.dat #
Contains the used parameters, and timing information

# Out_< name >.dat #
It contains nu and K(nu)
nu are positions of the wavenumber and K is the unsorted opacity function

# Out_< name >_bin.dat #
It contains the y and K(y) per bin
y goes from 0 to 1. K(y) is the per bin sorted opacity function. The bins are separated by two blank lines, starting with 
the bin with the lowest wavenumbers and ending with the bin with the highest wavenumbers.
When doResampling is set to one, then this file contains the sorted opacity functions computed from the 
Chebyshev coefficients

# Out_< name >_cbin.dat #
It contains the chebyshev coefficients of the per bins resampled sorted logarithm of the opacity functions in the format
kmin_i ystart_i C0_i C1_i ... C(nC - 1)_i
where i refers to the bin index, and C are the Chebyshev coefficients
kmin is the minimal value of K(y), reached when cutting the Voigt profiles
ystart is the position in y when the K(y) starts to be larger than kmin
K(y) can be recomputed as K(y) = sum_(0 <= j < nC) (C[j] * T[j](yy)), where T(y) are the Chebyshev polynomials,
where yy = (2.0 * y - 1.0 - ystart) / (1.0 - ystart), for y in the range [ystart, 1]
The bins are separated with a blank line, starting with 
the bin with the lowest wavenumbers and ending with the bin with the highest wavenumbers.

# Out_< name >_tr.dat #
It contains m and T.
m is the column mass, m_i = exp((i - nTr/2) * dTr)
T is the Transmission function Int_0^1 exp(-K(y)m) dy

# Out_< name >_mean.dat #
When the argument doMean is set to one, this file contains the Planck and Rosseland means.
They are computed over the entire range in wavenumbers from numin to numax with spacing dnu.  
The first line is the Planck mean: kappa_P = Int_0^infty (kappa * B_nu * dnu) / Int_0^infty (B_nu * dnu)  
The second line is the Rosseland mean: kappa_R = (Int_0^infty (kappa^-1 * del(B_nu)/del(T) * dnu) / Int_0^infty ( del(B)/del(T)_nu * dnu))^-1  
The third line is the numerical integral Int_0^infty (B_nu * dnu)  
The fourth line is the analytic integral Int_0^infty (B_nu * dnu) = sigma * T^4 / pi
The fifth line is the numerical integral Int_0^infty ( del(B)/del(T)_nu * dnu)  
The sixth line is the analytic integral Int_0^infty ( del(B)/del(T)_nu * dnu) = 4 * sigma * T^3 / pi


The value of the numerical integrals should converge to the analytic expressions for high resolutions dnu, numin -> 0 and numax -> infinity.


