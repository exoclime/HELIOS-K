Cut the line wings
==================


The line wings can be cut by the ``cutMode`` and ``cut`` parameters. Reasons
for cutting the line wings can either be that the line wings would be 
overestimated, or just to save computing time. Line wings or cut at the 
specified distance from their central peak location. Choosing a broader cutting
lenght can increase the computational time.
In :numref:`figcut` is shown an example with three different cutting lenghts. 


| Relevant parameters for this example:

 - doStoreFullK = 1
 - cutMode = 0
 - cut = 25, 100 or 0

 

.. figure:: ../../plots/cut/plot001.png  
   :name: figcut

   Three cutting lengths of the line wings. An infinity cuttoff lengths can be set by `cut = 0`.
    

