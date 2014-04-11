pyrings
=======

`PyRings <http://pyrings.castelao.net>`_ is a Python package to handle coherent vortex (rings and eddies) in the ocean.

A set of functions and classes to deal with rings/eddies identified by the velocity field. It does not require any specific data distribution, i.e. works fine with irregularly distributed data, like ship tracks, or even with hybrid mixed dataset composed from different instruments.

The center is estimated and the data is re-estructured to a natural coordinate system.

This concept was developed during my PhD, and I'm planning to release it in three major blocks. This first one is just to identify the center and transform the usual cartesian data distribution into a cylindrical system for better suit the Cyclo-Geostrophic approximation.

Reference
---------

If you use it, please cite us: http://dx.doi.org/10.1016/j.cageo.2013.07.004

| @article{castelao2013,
|  title = "An objective reference system for studying rings in the ocean",
|  author = "Guilherme P. Castel{\~a}o and Luiz C. Irber Jr and Ana Villas Boas",
|  journal = "Computers \& Geosciences",
|  volume = "61",
|  pages = "43--49",
|  year = "2013",
|  publisher = "Pergamon"
|  }

Acknowledgment
--------------

Thanks to Bill Johns who initially introduce me to this problem, and Roberto de Almeida.
