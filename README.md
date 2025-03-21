# CTDYN

This repository is a fork of the Catania dynamo code project originally
developed by Alfio Bonanno. The code has been refactored in depth in
order to upgrade its structure and syntax towards modern Fortran standards.
A Python module dedicated to analysis and display of code outputs is also
included with this distribution.  

## Installation

### Requirements 

In order to compile the Fortran executable, you should make sure that the
necessary requirements are available on your system. You will need:
- a recent version [gfortran](https://gcc.gnu.org/fortran/) compiler ;
- the [Fortran Standard Library](https://github.com/fortran-lang/stdlib) ;
- an installation of the [LAPACK](https://www.netlib.org/lapack/) linear algebra library. 
CTDYN will perform better if use an optimised library  (e.g. the 
[openBLAS](http://www.openmathlib.org/OpenBLAS/) library). 

### Fortran CTDYN executable

Before installing CTDYN, you will need to set the `CTDYN_DIR` environment variable
in your local profile:

`$ export CTDYN_DIR="[CTDYN_repository]"`  

The next step is to run the installation script to compile the Fortran source code:

`$ ./install`

An executable should be created in the `bin` repository under the name `ctdyn`.
When re-installing CTDYN, you might want to delete the `bin` repository beforehand
to avoid interferences from previous `cmake` configurations.

The install scripts accepts two verbosity levels:

`$ ./install -v`
`$ ./install -V`

The default build is set to be the release build, you might want to build in debug mode:

`$ ./install -b debug`

### Python module

The *py_ctdyn* module provides a library of Python tools designed to analyse
and display the outputs of CTDYN computations. The module can be installed by
running the following command at the root of the repository. 

`$ pip install .`
 
## Get started

The executable can be run by selecting any correctly formatted input file, e.g.
by doing

`$ $CTDYN_DIR/bin/ctdyn test/inputs/inlist_default`

CTDYN can also be run directly within a Python session using the API provided
by *py_ctdyn*, an example is provided in the `quickstart.ipynb` notebook.

## Documentation

A documentation for the code is available [there](https://ctdyn.netlify.app/) 
(work in progress!).
 
## Contributors

- Sylvain N. Breton
- Alfio Bonanno 
- Giovanni Licciardello 

## Useful references

- [An accurate numerical approach for the kinematic dynamo problem](https://ui.adsabs.harvard.edu/abs/2004MSAIS...4...17B/abstract),
Bonanno 2004.
- [A solar mean field dynamo benchmark](https://ui.adsabs.harvard.edu/abs/2008A%26A...483..949J/abstract), Jouve et al. 2008.

