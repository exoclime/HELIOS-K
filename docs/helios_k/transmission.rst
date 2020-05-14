The transmission function
=========================

.. _transmission:

With the ``doTransmission`` option in the ``param.dat`` file, the transmissions function
:math:`\tau` 

.. math::

   \tau = \int_0^\infty e^{-\kappa \tilde{m}} d\nu = \int_0^1 e^{-\kappa \tilde{m}} dy

can be calculated for a set of column masses :math:`\tilde{m}`, where :math:`\tilde{m}`
is set according to:

.. math::

    \tilde{m}_i = e^{((i - nTr/2) \cdot dTr)}.

The parameters ``nTr`` and ``dTr`` can be set in the ``param.dat`` file.
The transmission function is stored in the file ``Out<name>_tr.dat``,


An example of the transmission function for fife bins is shown in :numref:`figtransmission`.


| Relevant parameters for this example:

 - doStoreFullK = 1
 - doTransmission = 2
 - nbins = 5
 - nTr = 1000
 - dTr = 0.05


.. figure:: ../plots/p006/plot001.png  
   :name: figtransmission

   Transmission function 


