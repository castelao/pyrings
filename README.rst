pyrings
=======

`PyRings <http://pyrings.castelao.net>`_ is a Python package to handle coherent vortex (rings and eddies) in the ocean.

A set of functions and classes to deal with rings/eddies identified by the velocity field. It does not require any specific data distribution, i.e. works fine with irregularly distributed data, like ship tracks, or even with hybrid mixed dataset composed from different instruments.

The center is estimated and the data is re-estructured to a natural coordinate system.

This concept was developed during my PhD, and I'm planning to release it in three major blocks. This first one is just to identify the center and transform the usual cartesian data distribution into a cylindrical system for better suit the Cyclo-Geostrophic approximation.

Acknowledgment
--------------

Thanks to Bill Johns who initially introduce me to this problem, and Roberto de Almeida.
