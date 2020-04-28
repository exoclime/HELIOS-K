# HELIOS-K #
HELIOS-K is an opacity calculator, running on GPUs.
#### Authors: Simon Grimm, Kevin Heng ####

# Updates #
## version 1.67 ##
The `Molecule` and `useHITEMP` arguments in the `param.dat` file are not valid anymore. Species must now be set by `Species Name`. The database is written from the `< species >.param` file.

## version 1.65 ##
All species must have now a `< species >.param` file, which contains all necessary information about the line list.

# Requirements #

HELIOS-K runs on Nvidia GPUs with compute capability of 2.0 or higher.
The code needs the CUDA toolkit to be installed.This can be downloaded
from https://developer.nvidia.com/cuda-downloads.

Besides the main code, HELIOS-K provides various scripts to download and pre-process line list files. These scripts have additional dependencies:

When running on a cluster, these libraries can be installed locally with
```
pip3 install --user <package name>
```

 * exomol.py need the following python3 libraries:
    * bs4
    * requests
    * sys
    * os
    * subprocess
    * numpy
    * argparse
 
 * Kurucz2.py
     * numpy
     * math
     * struct
     * os
     * argparse
     
  * nist_ELevels.py
    * sys
    * time
    * pyperclip
    * subprocess
    * selenium
    * argparse


# Compilation #
HELIOS-K can be compiled with the provided Makefile by typing `make SM=xx` to the
terminal, where xx corresponds to the compute capability. For example use
`make SM=20`
for compute capability of 2.0, or `make SM=35` for 3.5. A table with all compute capabilities can be found here: https://developer.nvidia.com/cuda-gpus

## On Windows machines ##
If using Cygwin on Windows, then HELIOS-K can be compiled the same way with `make SM=xx`.
If using the Windows Command Prompt, type `nmake -f MakefileW SM=xx`. Note, that the Windows c++ compiler `cl` must be installed, and the compiler path must be loaded in the shell. If this is not the case, it can be loaded similar to:
`call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\Common7\Tools\vsvars32.bat"` .
See also issue #1 for details.



# Download and pre-process the line lists #
If you work an the University of Bern on the Hulk cluster, then most of the pre-processed files are already available at `scratch/sigrimm/EXOMOL/`.

## Supported line list databases ##
HELIOS-K supports line lists from the Hitran, HITEMP, ExoMol, Kurucz, NIST  or VALD
databases. Before the line lists can be used, they have to be pre-processed into binary
files. This saves in generally memory space and allows HELIOS-K to read the line lists in a more efficient way.

## The `< species >.param` file
Each molecular or atomic species, which will be processed by HELIOS-K, needs a
`< species >.param` file, which contains all necessary information about the line list database and other parameters.

Scripts are provided to download and pre-process some of the databases.

### Example 
An example for Hitran 2016 H2O is given in the `01_hit16.param` file, which is
shown here and available in the repository. An example for HITEMP and ExoMol are given in the repository files `01_HITEMP2010.param` and `1H2-16O__BT2.param`.

```
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
```
### Description
 * Database: It must be indicated which database format to use. 0 = Hitran or HITEMP, 2 = ExoMol, 30 = Kurucz atomic database, 31 = NIST atomic database.
 * Molecule number: This is optional and follows an old HELIOS-K version.
 * Name: name of the line list files.
 * Number of isotopologues in the line list files.
 * For each isotopologue:
    * ID: for Hitran and HITEMP the same as the first three digits in the `.par`files e.g. 11 for 1H2-18O, or 2A for 18O-13C-17O. For ExoMol no meaning.
    * Abundance, only relevant for more than one isotopologue
    * Reference partition function value. For Exomol no meaning.
    * Molar mass in g
    * Name of partition function file
 * Number of columns in partition file.
 * Number of line/transition files.
 * Line by line, the number of molecular/atomic lines in each line/transition files
 * Line by line, the wavenumber range of each line/transition files, must have one more entry than the number of line/transition files to cover all lower and upper range limits.
 * Number of states in ExoMol `.states`files. For other databases no meaning.
 * Number of columns in ExoMol `.trans`files. For other databases no meaning.
 * ExoMol default value of Lorentzian half-width. For other databases no meaning.
 * ExoMol default value of temperature exponent. For other databases no meaning.
 * Version of the line list.



## ExoMol ##

### Step 1, Species Properties ###
The HELIOS-K repository provides a file called `Exomol_species.dat`. This file contains all
available species of the ExoMol database. The file format is:

`Molecule name , Isotopologue name, Full name,  path on exomol.com, range of .trans files, number of .trans files, number of digits in .trans file ranges.`

The full name of the species contains the isotopologue and the line list name. This full name should be used, when specifying a species for the opacity calculation, e.g. `1H2-16O__BT2`.


The `Exomol_species.dat` file can be recreated or updated with the python code `exomol2.py`.

### Step 2 Download the files and create `< species >.param` file###
The ExoMol files can be downloaded with the python script `exomol.py` as:

```
python3 exomol.py -M <id>
``` 

where `<id>` is the full  species name e.g. `1H2-16O__BT2`. The script needs the file `Exomol_species.dat` to be available.
If the `Exomol_species.dat` file needs to be updated, it can be recreated by running `python3 exomol2.py`.

The `exomol.py` script automatically writes the `< species >.param` files for each molecule.

### Step 3 create binary files ###
The downloaded line list files must be pre-processed into binary files with the following code:

```
./prepareExomol -M < id >
```

where < id > is the full species name, e.g. `1H2-16O__BT2`.
After this step, the `.trans` and `.states` files from ExoMol can be deleted.

### Step 4 data path ###
Include the path of the directory, which contains the obtained binary files, the `.pf`file and the  `.param` file to the HELIOS-K `param.dat` file under `pathToData`.

## HITRAN ##

### Step 1  Species Properties ###
The HELIOS-K repository provides a file called `Hitran_species.dat`. This file contains all
available species of the Hitran database. The file format is:

`Species ID, Abundance, Molecular Mass, Q(296K), partition function file, Isotopologue Formula`
The same information can be found at `https://hitran.org/docs/iso-meta/`.

The species ID consitst of a two digits molecule ID and a one digit local isotopologue ID. Note that the local isotogologue ID can sometimes consist of non numerical values. e.g A or B,

For identifying a species, the molecule number and the isotopologue number should be combinded, e.g. `01_1` for 1H2-16O or `01` for all isotopogolues from H20.


The `Hitran_species.dat` file can be recreated or updated with the python code `hitran2.py`.

### Step 2 Download the files ###
The line list files must be downloaded manually from `www.hitran.org`, and note that it is necessary to register on the hitran homepage. To download the files select DataAcess and then Line-by-line. Select the molecule id, select all isotopologues (single isotopologues can be filtered later), leave wavenumber range blank, select `.par`file and store the file on your computer with the name `Molecule-ID_hit16.par`, e.g. 01_hit16.par for H2O.  


Download all the necessary partition function files from Documentation/Isotogologues.


### Step 3 create `< species >.param` file and binary files ###
All necessary files can be created with:
```

./hitran -M < molecule ID > -ISO < isotopologue ID > -in < line list name > 
```
The `<molecule ID>` is the two digit molecule number, e.g. `01` for H2O.
The `<isotopologue ID>` is the hitran internal isotopologue identifier, e.g. `1` for 1H2-16O.
The `<line list name>` is the name that was given in the download section, e.g. `hit16`.


### Step 4 data path ###
Include the path of the directory, which contains the obtained binary files, the `.txt` partition function files and the  `.param` file to the HELIOS-K `param.dat` file under `pathToData`.


## HITEMP ##

### Step 1  Species Properties ###
For HITEMP, the same `Hitran_species.dat` is needed as for Hitran. See section Hitran step 1


### Step 2 Download the files ###
The line list files must be downloaded manually from `https://hitran.org/hitemp/`.  The same partition function files are needed as for Hitran, see section Hitran step 2.


### Step 3 create `< species >.param` file and binary files ###
All necessary files can be created with:
```

./hitran -M < molecule ID > -ISO < isotopologue ID > -in < line list name > 
```
The `<molecule ID>` is the two digit molecule number, e.g. `01` for H2O.
The `<isotopologue ID>` is the hitran internal isotopologue identifier, e.g. `1` for 1H2-16O.
The `<line list name>` is the name that was given in the download section, e.g. `HITEMP2010`.


### Step 4 data path ###
Include the path of the directory, which contains the obtained binary files, the `.txt` partition function files and the  `.param` file to the HELIOS-K `param.dat` file under `pathToData`.


## Kurucz Database ##

###  Step 1 Download the files and create binary files and <species >.param files ###
The needed files can be downloaded with the following comand:
```
python3 Kurucz2.py -D 1 -Z <z> -I <i>
```
where `<z>`is the atomic number and `<i>`the ionization state. z = -1 means to download all atoms, i = -1 means to download Z = 0, 1 and 2.
-D 1, means to download the file gfallwn08oct17.dat and all available partition function files.
-D 2, means to download only the partition function files.
-D 0, means no download from the Kurucz website. 

If the source file `gfallwn08oct17.dat` should be changed, then the script `Kurucz2.py` must be modified. 
This script will also create the `.param`and the `.bin` file.

### Step 2 data path ###
Include the path of the directory, which contains the obtained binary files, the `.pf` partition function files and the  `.param` file to the HELIOS-K `param.dat` file under `pathToData`.


## NIST Database ##

### Step 1 Download the Energy levels ###
This step can be done either manually or by using a script. 

 * For manuall download:
    * visit `https://physics.nist.gov/PhysRefData/ASD/levels_form.html`
    * enter the species and ionization state, e.g. `Z= 3 0`
    * select the Format output to Tab-delimited
    * select the `g` option
    * click on `Retrieve Data`
    * store the data in a file called `NIST_ELevels<Z><I>.dat`
    * remove all `"` in the file
 * Automatical download:
    * type `python3 nist_ELevels.py -Z <z> -I <i> `, 
       where `<z>`is the atomic number and `<i>`the ionization state, e.g. -Z 3 -I 0
    * The script will download the data and store it in a file `NIST_ELevels<Z><I>.dat`.
    * This script will need geckodriver to be installed.
 
### Step 2 Generate partition functions ###

Run  `python3 nist_partition.py -Z <z> -I <i> `. This script will
read the previously produces file `NIST_ELevels<Z><I>.dat`, and writes a '*.pf' file.

### Step3 Download atomic masses ###
 * visit `https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses`
 * select `All Elements`
 * select `Linearized ASCII Output`
 * select `All isotopes`
 * click `Get Data`and store the data in a file called `masses.txt`
 
### Step 4 Download the line list ### 
This step can be done either manually or by using a script. 

 * For manuall download:
    * https://physics.nist.gov/PhysRefData/ASD/lines_form.html
    * enter the species and ionization state in the `Spectrum field`, e.g. `Z= 3 0`
    * select `Show Advanced Settings`
    * select the Format output to csv
    * Select `Wavenumbers (in cm-1)` (keep `Observed`, `Ritz` and `Uncertainties`selected)
    * Select `Wavenumbers (all wavelenghts)`
    * Select `g`
    * click `Retrieve Data`and store the data in a file called `NIST_Lines<Z><I>.dat`.
    * remove all  `?, =, [,],( and )`from the file
 * Automatical download:
    * type `python3 nist_Lines.py -Z <z> -I <i> `.
    * The script will download the data and store it in a file `NIST_Lines<Z><I>.dat`.
    * This script will need a geckodriver to be installed .

### Step 5 Create `< species >.param` file and binary files ###
All necessary files can be created with:
```

python3 nist_Lines2.py -Z <z> -I <i> 
```

### Step 6 data path ###
Include the path of the directory, which contains the obtained binary files, the `.pf` partition function files and the  `.param` file to the HELIOS-K `param.dat` file under `pathToData`.

# Using HELIOS-K #
Before HELIOS-K can be used, the molecular or atomic line-lists must be downloaded and pre-processed. HELIOS-K provides some scripts for that, which are described later in in this manual.

All necessary parameters for HELIOS-K are set in the file `param.dat`, which is also described later. Some parameters can be set also by console arguments.


HELIOS-K can be started with

```
./heliosk
```
followed by optional console arguments, which are listed below.



# HELIOS-K Input parameters #
The input parameters can be specified in the `param.dat` file. The used
parameters are listed here, the order can not be changed.

 * name: This name will appear in the output filenames.
 * T: Temperature in Kelvin
 * P: Pressure in Atmospheres
 * PFile: A '-' ignores this option, otherwise this option specifies a filename which contains multiple values for P
 * Species Name: The name of the molecule or atomic param file. e.g. 01_hit16 or 1H2-16O__BT2
 * SpeciesFile: A '-' ignores this option, otherwise this option specifies a filename which contains multiple species for gas mixtures
 * ciaSystem: A '-' ignores this option, otherwise a cia system is read here. supported are H2-H2, H2-H2_eq, H2_H2_norm, H2-He, H2-He_eq,
  H2-He_norm, H2-CH4_eq, H2-CH4_norm and H2-H.
 * pathToData: The location where the HITRAN or HITEMP data files are located, e.g. pathToData = ../HITEMP/ , pathToData = /data/HITEMP/ or empty when the files are in the same directory  pathToData = 
 * numin: minimum wavenumber in 1/cm
 * numax: maximum wavenumber in 1/cm
 * dnu: spatial resolution in wavenumber in 1/cm
 * Nnu per bin: The number of points in nu per bin. When this is set to zero, then dnu is used, otherwise the value of dnu is overwritten.
 * cutMode: Set the voigt profile cutting mode:
   * 0: cut at absolute wavenumbers
   * 1: cut at multiple values of Lorentz widths
   * 2: cut at multiple values of Doppler widths
 * cut: value where to cut, according to the cutMode, if cut = 0, then no cutting is performed
 * doResampling: When set to one, then all the sorted opacity functions per bin are resampled with a Chebyshev polynomial, and the coefficients are written to the file `Out_<name>_cbin.dat`.
 *nC: Number of Chebyshev coefficients use for the resampling
 * doTransmission: When set to one, then the transmission function per bin is computed for various numbers of m, and written to the file `Out_<name>_tr.dat`.
 * nTr: nummber of points used for the transmission function
 * dTr: spacing of the points used for the transmission function in exp space: m_i = exp((i - nTr/2) * dTr)
 * doStoreFullK: When set to one, then the full unsorted opacity function is written to the file `Out_<name>.dat`.
 * doStoreSK: When set to 1, then the per bin sorted opacity function is written to the file `Out_<name>_bin.dat`. When set to 2, then the per bin sorted opacity function is written to a different bin files `Out_<name>_bin< bin number >.dat`.
 * nbins: number of bins
 * binsfile: A '-' ignores this option, otherwise this option specifies a filename which contains the edges of the bins, which can be irregular. This option overrides the numin, numax and nbins arguments.
 * OutputEdgesFile: A '-' ignores this option, otherwise this option specifies a file name which contains the edges of the output locations in y for each bin.
 * kmin: minimal value for the opacity function
 * qalphaL: q value in the Lorentzian half width q = Pself / P
 * gammaF: scaling factor for the Lorentzian half width gamma factor, gamma = gamma * gammaF
 * doMean: Calculate the Planck and Rosseland opacity means
 * Units: The units of the opacities. 0: cm^2 / g, 1: cm^2 / molecule
 * ReplaceFile: When set to 1, then all output files are overwritten, when set to 0 then the data is appended to the existing files.
 * profile: 1 = Voigt, 2 = Lorentz, 3 = Gauss, 4 = cross section
 * doTuning: 1 = use self tuning routines to specify the kernels with best parameters.

# Console Arguments #
Instead of using the parameter file, some arguments can also be passed as console arguments. The console arguments have the highest priority and are overwriting the arguments of the `param.dat` file. The options are:

 * \-name `<c>`: name
 * \-T `<double>` : T
 * \-P `<double>` : P
 * \-M `<int>` : Molecule Name
 * \-path `<c>` : pathToData
 * \-pathK `<c>` : pathToK
 * \-numin `<double>` : numin
 * \-numax `<double>` : numax
 * \-dnu `<double>` : dnu
 * \-cutM `<int>` : cutMode
 * \-cut `<double>` : cut
 * \-dR `<int>` : doResampling
 * \-nC `<int>` : nC
 * \-dT `<int>` : doTRansmission
 * \-nTr `<int>` : nTr
 * \-dTr `<double>` : dTr
 * \-dSF `<int>` : doStoreFullK
 * \-dSS `<int>` : doStoreSK
 * \-nbins `<int>` : nbins
 * \-kmin `<double>` : kmin
 * \-dev `<int>` : Device number (For multiple GPU systems)
 * \-q `<double>` : qalphaL
 * \-gammaF `<double>` : gammaF
 * \-Mean `<int>` : doMean
 * \-tuning `<int>` : doTuning
where `<c>`is a string, `<double>` a floating point number, and `<int>`an integer.

# Code parameters #
The file define.h contains the physical parameters and some code parameters.
After changing some entry, the Code needs to be recompiled.
The code parameters are:

 * def_TOL:         Tolerance parameter in the Voigt function. See Algorithm 916
 * def_nthmax:      Maximum number of threads per kernel launch. In 2.0 GPUs it can not be larger than 32768.
 * def_nlmax:       Maximum number of molecular lines per kernel launch. Setting a lower number prevents from a time-out on Desktop machines.
 * def_maxlines:    Maximum number of lines stored on the GPU.
 * def_maxfiles:    Maximum number of files per molecule.
 * def_NmaxSample:	Maximum Number of resample coefficients for K(y)

When using a Desktop GPU running an x session, the runtime of a single kernel
launch can be limited to a few seconds. Choosing smaller values for nlmax and nthmax
splits the kernel into smaller parts. But it makes the code a bit slower.


# The binsfile options #
When a `binsFile` name is given in the `param.dat` file, then this file is used to generate the edges of the bins, which can be irregular. Note that this option
does no support the doResampling and doTransmission options.
The binsfile must contain line by line the edges of the bins in cm^-1.

For example:
```
0.5
50.0
100.0
200.0
```

# The output edges option #
When a `outputedgesFile` name is given in the `pram.dat` file, then this file is used to specify the averaged output positions of the `Out_<name>_bin.dat` files. The file must contain line by line the positions in y.

For example:
```
0.0
0.1
0.45
0.78
1.0
```

# The P file option #
When a `PFile` name is given in the `param.dat` file, then this file is used to read multiple values for P. This option is useful to speed up the performance, 
because multiple reads from the data files can be avoided. Too many entries in the Pfile can lead to a memory shortage.

For example:
```
1.0
10.0
100.0
```

# The Species file option #
This option must be used to calculate opacities for gas mixtures, containing multiple species.
The File contains in the two columns the species name, and the number fraction.

For example:
```
01_hit16	0.9
05_hit16	0.1
```
This example will produce an opacitiy with 90% H2O and 10% CO.

# Output Files #
Different Output files are written, depending to the set values in the `param.dat` file

# `Info_<name>.dat` #
Contains the used parameters, and timing information

# `Out_<name>.dat` #
It contains nu and K(nu), where
nu are the wavenumbers  and K(nu) is the full opacity function.
When the `PFile` option is used, then the files contain also the values of T and P.

# `Out_<name>_bin.dat` #
It contains the values of y and K(y) per bin. y goes from 0 to 1. K(y) is the per bin sorted opacity function. The bins are separated by two blank lines, starting with the bin with the lowest wavenumber and ending with the bin with the highest wavenumber.

When `doResampling` is set to one, then this file contains the sorted opacity functions, recomputed from the Chebyshev coefficients.
When the `PFile` option is used, then the files contains also the values of T, P and point index.

When the `OutputEdgesFile` option is used, then the file contains not all points in y, but the averaged values between the edges, given in the `OutputEdgesFile`.

When `doStoreSK` is set to 2, then the bins are stored in different files with names `Out_<name>_bin< bin index>.dat` .

# `Out_<name>_cbin.dat` #
It contains the Chebyshev coefficients of the per bins sorted
natural logarithm of the opacity functions in the format

 `kmin_i ystart_i C0_i C1_i ... C(nC - 1)_i`,
 where i refers to the bin index, and `C` are the Chebyshev coefficients.

`kmin` is the minimal value of `K(y)`, used e.g. in holes in the opacity funtion.

`ystart` is the position in y when the value of`K(y)` starts to be larger than `kmin`.

`K(y)` can be recomputed as

`K(y) = sum_(0 <= j < nC) (C[j] * T[j](yy))`,

 where `T(y)` are the Chebyshev polynomials and `yy = (2.0 * y - 1.0 - ystart) / (1.0 - ystart)`, for y in the range `[ystart, 1]`.
The bins are separated with a blank line, starting with the bin with the lowest wavenumber and ending with the bin with the highest wavenumber.
When the `PFile` option is used, then the files contains also the values of T and P.
When `doResampling` is set to 2, then the bins are stored in different files with names `Out_<name>_cbin< bin index>.dat` .

The following python script can be used to reconstruct the per bin sorted opacity function from the Chepyshev coefficients:

```
import numpy as np
from numpy.polynomial.chebyshev import chebval

#change here the name of the file
data_c = np.loadtxt('Out_name_cbin.dat')

#change here the bin index and the bin size:
binIndex = 0
binSize = 300

#extract Chebyshev coefficients
c = data_c[binIndex,2:]
#extract starting point in x of opacity function
xs = data_c[binIndex,1]

#rescale x to the standard Chebychev polynomial range [-1:1]
x1 = x * 2.0 - 1.0
k_res = chebval(x1,c,tensor=False)
x2 = x * (1.0 - xs) + xs

#result is in k_res for x values in x2
k_res = np.exp(k_res)

```
 

# `Out_<name>_tr.dat` #
It contains m and T.
m is the column mass, m_i = exp((i - nTr/2) * dTr)
T is the Transmission function Int_0^1 exp(-K(y)m) dy
When the PFile is used then the files contains also the values of T, P and point index.
When doTransmission is set to 2, then the bins are stored in different files with names `Out_<name>_tr< bin index>.dat`

# `Out_<name>_mean.dat` #
When the argument doMean is set to one, this file contains the Planck and Rosseland means.
They are computed over the entire range in wavenumbers from numin to numax with spacing dnu.
The first line is the Planck mean: kappa_P = Int_0^infty (kappa * B_nu * dnu) / Int_0^infty (B_nu * dnu).
The second line is the Rosseland mean: kappa_R = (Int_0^infty (kappa^-1 * del(B_nu)/del(T) * dnu) / Int_0^infty ( del(B)/del(T)_nu * dnu))^-1  
The third line is the numerical integral Int_0^infty (B_nu * dnu)  
The fourth line is the analytic integral Int_0^infty (B_nu * dnu) = sigma * T^4 / pi
The fifth line is the numerical integral Int_0^infty ( del(B)/del(T)_nu * dnu)  
The sixth line is the analytic integral Int_0^infty ( del(B)/del(T)_nu * dnu) = 4 * sigma * T^3 / pi


The value of the numerical integrals should converge to the analytic expressions for high resolutions dnu, numin -> 0 and numax -> infinity.


