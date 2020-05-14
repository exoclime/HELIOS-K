Planck and Rosseland Means
==========================

.. _mean:

With the ``doMean`` option in the ``param.dat`` file, the Plank and Rosseland means, 
:math:`\kappa_P` and :math:`\kappa_R` respectively: 

.. math::
   :label: eq_a

   \kappa_P = \frac{\int_0^\infty \kappa  B_{\nu}(T) d\nu}{\int_0^\infty B_{\nu}(T)  d\nu},

and
 
.. math::
   :label: eq_b

   \kappa_R = \left(  \frac{\int_0^\infty \kappa^{-1} \frac{\partial B_{\nu}(T)}{\partial T} d\nu}  {\int_0^\infty  \frac{\partial B_{\nu}(T)}{\partial T} d\nu} \right)^{-1} 

can be calculated, where the infinity integral is truncated to the specified wavenumber limits, and :math:`d\nu` is 
equal to ``dnu`` set in the ``param.dat`` file.

Note that the denominators in :math:`\kappa_P` and :math:`\kappa_R` can be computed analytically as

.. math::
   :label: eq_0

   \int_0^\infty B_{\nu}(T) d\nu = \frac{\sigma  T^4} {\pi}

and

.. math::
   :label: eq_1

   \int_0^\infty  \frac{\partial(B_{\nu}(T))}{ \partial (T)}  d\nu = \frac{4  \sigma  T^3} {\pi}.

Therefor, it is usefull to compare those analytic results to the numerical integration

.. math::
   :label: eq_2

   \int_0^\infty B_{\nu}(T) d\nu = \sum_i  B_{\nu}(T) dnu

and 

.. math::
   :label: eq_3

   \int_0^\infty  \frac{\partial(B_{\nu}(T))}{ \partial (T)}  d\nu =  \sum_i  \frac{\partial(B_{\nu}(T))}{ \partial (T)}  dnu



The results of the Planck and Resseland means are stored in the file ``Out<name>_mean.dat``, together with 
the analytic and numerical expressions :eq:`eq_0` to :eq:`eq_3`. If the numerical expressions deviate strongly
from the analytical expression, then it is a hint that the wavenumber resolution is not set fine enough.


An example of the Planck and Rosseland means is shown in :numref:`figmean`.


| Relevant parameters for this example:

 - doStoreFullK = 1
 - doMean = 1


.. figure:: ../plots/p007/plot001.png  
   :name: figmean

   Planck and Rosseland means


