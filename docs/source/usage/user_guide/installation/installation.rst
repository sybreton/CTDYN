Installation
#############

Requirements 
------------

In order to compile the Fortran executable, you should make sure that the
necessary requirements are available on your system. You will need:
* a recent version `gfortran <https://gcc.gnu.org/fortran/>`_ compiler ;
* the `Fortran Standard Library <https://github.com/fortran-lang/stdlib>`_ ;
* an installation of the `LAPACK <https://www.netlib.org/lapack/>`_` linear algebra library. 
CTDYN will perform better if use an optimised library  (e.g. the 
`openBLAS <http://www.openmathlib.org/OpenBLAS/>`_ library). 

Fortran CTDYN executable
-------------------------

Before installing CTDYN, you will need to set the `CTDYN_DIR` environment variable
in your local profile:

.. code-block:: console

  $ export CTDYN_DIR="[CTDYN_repository]" 

The next step is to run the installation script to compile the Fortran source code:

.. code-block:: console

  $ ./install

An executable should be created in the `bin` repository under the name `ctdyn`.
When re-installing CTDYN, you might want to delete the `bin` repository beforehand
to avoid interferences from previous `cmake` configurations.

The install scripts accepts two verbosity levels:

.. code-block:: console

  $ ./install -v
  $ ./install -V

The default build is set to be the release build, you might want to build in debug mode:

.. code-block:: console

  $ ./install -b debug


Python module
--------------

The *py_ctdyn* module provides a library of Python tools designed to analyse
and display the outputs of CTDYN computations. The module can be installed by
running the following command at the root of the repository.

.. code-block:: console

  $ pip install .

