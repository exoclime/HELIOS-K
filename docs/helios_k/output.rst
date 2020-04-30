Output Files
============

Different Output files are written, depending to the values in the
``param.dat`` file

.. _info_<name>.dat:

``Info_<name>.dat``
===================

Contains all used parameters, and timing information

.. _out_<name>.dat:

``Out_<name>.dat``
==================

It contains nu and K(nu), where nu are the wavenumbers and K(nu) is the
full opacity function. When the ``PFile`` option is used, then the files
contain also the values of T and P.

.. _out_<name>_bin.dat:

``Out_<name>_bin.dat``
======================

It contains the values of y and K(y) per bin. y goes from 0 to 1. K(y)
is the per bin sorted opacity function. The bins are separated by two
blank lines, starting with the bin with the lowest wavenumber and ending
with the bin with the highest wavenumber.

When ``doResampling`` is set to one, then this file contains the sorted
opacity functions, recomputed from the Chebyshev coefficients. When the
``PFile`` option is used, then the files contain also the values of T,
P and point index.

When the ``OutputEdgesFile`` option is used, then the file contains not
all points in y, but the averaged values between the edges, given in the
``OutputEdgesFile``.

When ``doStoreSK`` is set to 2, then the bins are stored in different
files with names ``Out_<name>_bin< bin index>.dat`` .

.. _out_<name>_cbin.dat:

``Out_<name>_cbin.dat``
=======================

It contains the Chebyshev coefficients of the per bins sorted natural
logarithm of the opacity functions in the format

``kmin_i ystart_i C0_i C1_i ... C(nC - 1)_i``, where i refers to the bin
index, and ``C`` are the Chebyshev coefficients.

``kmin`` is the minimal value of ``K(y)``, used e.g. in holes in the
opacity funtion.

``ystart`` is the position in y when the value of\ ``K(y)`` starts to be
larger than ``kmin``.

``K(y)`` can be recomputed as

``K(y) = sum_(0 <= j < nC) (C[j] * T[j](yy))``,

where ``T(y)`` are the Chebyshev polynomials and
``yy = (2.0 * y - 1.0 - ystart) / (1.0 - ystart)``, for y in the range
``[ystart, 1]``. The bins are separated with a blank line, starting with
the bin with the lowest wavenumber and ending with the bin with the
highest wavenumber. When the ``PFile`` option is used, then the files
contains also the values of T and P. When ``doResampling`` is set to 2,
then the bins are stored in different files with names
``Out_<name>_cbin< bin index>.dat`` .

The following python script can be used to reconstruct the per bin
sorted opacity function from the Chebyshev coefficients:

::

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

.. _out_<name>_tr.dat:

``Out_<name>_tr.dat``
=====================

It contains m and T. m is the column mass, m_i = exp((i - nTr/2) \* dTr)
T is the Transmission function Int_0^1 exp(-K(y)m) dy When the PFile is
used then the files contains also the values of T, P and point index.
When doTransmission is set to 2, then the bins are stored in different
files with names ``Out_<name>_tr< bin index>.dat``

.. _out_<name>_mean.dat:

``Out_<name>_mean.dat``
=======================

| When the argument doMean is set to one, this file contains the Planck
  and Rosseland means. They are computed over the entire range in
  wavenumbers from numin to numax with spacing dnu. 
| The first line is the Planck mean: kappa_P = Int_0^infty (kappa \* B_nu \* dnu) /
  Int_0^infty (B_nu \* dnu). 
| The second line is the Rosseland mean:
  kappa_R = (Int_0^infty (kappa^-1 \* del(B_nu)/del(T) \* dnu) /
  Int_0^infty ( del(B)/del(T)_nu \* dnu))^-1
| The third line is the numerical integral Int_0^infty (B_nu \* dnu)
| The fourth line is the analytic integral Int_0^infty (B_nu \* dnu) =
  sigma \* T^4 / pi 
| The fifth line is the numerical integral Int_0^infty
  ( del(B)/del(T)_nu \* dnu)
| The sixth line is the analytic integral Int_0^infty ( del(B)/del(T)_nu
  \* dnu) = 4 \* sigma \* T^3 / pi

The value of the numerical integrals should converge to the analytic
expressions for high resolutions dnu, numin -> 0 and numax -> infinity.

