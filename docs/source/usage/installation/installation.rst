Installation
#############

Fortran CTDYN executable
-------------------------

The first thing is to compile the Fortran source code.

``make``

An executable should be created in the ``bin`` repository under the name
``ctdyn``. It is recommended to set a ``CTDYN_DIR`` environment variable
in your local profile

``export CTDYN_DIR="CTDYN_repository"``  

Python module
--------------

The *py_ctdyn* module provides a library of Python tools designed to analyse
and display the outputs of CTDYN computations. The module can be installed by
running the following command at the root of the repository.

``pip install .``

