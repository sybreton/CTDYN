Installation
#############

Fortran CTDYN executable
-------------------------

Before compiling the Fortran executable, you should make sure that the
necessary pre-requisit, in particular the `gfortran <https://gcc.gnu.org/fortran/>`_
compiler and the `LAPACK <https://www.netlib.org/lapack/>`_ linear algebra library
are installed on your system. A convenient way to get these pre-requisite is to
install the `MESA Software Development Kit <http://user.astro.wisc.edu/~townsend/static.php?ref=mesasdk>`_
which provides them as a bundle. The MESA SDK is the recommended way to get the
installation pre-requisites for several commonly-used codes in stellar astrophysics,
such as `MESA <https://docs.mesastar.org/>`_ or `GYRE <https://gyre.readthedocs.io/>`_.

``$ make``

An executable should be created in the ``bin`` repository under the name
``ctdyn``. It is recommended to set a ``CTDYN_DIR`` environment variable
in your local profile

``$ export CTDYN_DIR="CTDYN_repository"``  

Python module
--------------

The *py_ctdyn* module provides a library of Python tools designed to analyse
and display the outputs of CTDYN computations. The module can be installed by
running the following command at the root of the repository.

``$ pip install .``

