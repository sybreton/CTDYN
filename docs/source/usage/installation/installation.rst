Installation
#############

Fortran CTDYN executable
-------------------------

The first thing is to compile the source code.
At the root of the repository, simply do

``cd src/CTDYN``

``make``

An executable should be created at the root of the repository under the name
``ctdyn``.

Python module
--------------

The *py_ctdyn* module provides a library of Python tools designed to analyse
and display the outputs of CTDYN computations. The module can be installed by
running the following command at the root of the repository.

``pip install .``

