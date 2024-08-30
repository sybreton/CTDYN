# CTDYN

This repository is a fork of the Catania dynamo code project originally
developed by Alfio Bonanno. The code has been refactored in depth in
order to upgrade its structure and syntax towards modern Fortran standards.
A Python module dedicated to analysis and display of code outputs is also
included with this distribution.  

## Installation

### Fortran CTDYN executable

Before compiling the Fortran executable, you should make sure that the
necessary pre-requisit, in particular the [gfortran](https://gcc.gnu.org/fortran/) 
compiler and the [LAPACK](https://www.netlib.org/lapack/) linear algebra library 
are installed on your system. A convenient way to get these pre-requisite is to 
install the [MESA Software Development Kit](http://user.astro.wisc.edu/~townsend/static.php?ref=mesasdk)
which provides them as a bundle. The MESA SDK is the recommended way to get the
installation pre-requisites for several commonly-used codes in stellar astrophysics,
such as [MESA](https://docs.mesastar.org/) or [GYRE](https://gyre.readthedocs.io/).  

The next step is to compile the Fortran source code

`make`

An executable should be created in the `bin` repository under the name
`ctdyn`. It is recommended to set a `CTDYN_DIR` environment variable
in your local profile

`export CTDYN_DIR="CTDYN_repository"`  

### Python module

The *py_ctdyn* module provides a library of Python tools designed to analyse
and display the outputs of CTDYN computations. The module can be installed by
running the following command at the root of the repository. 

`pip install .`
 
## Get started

The executable can be run by selecting any correctly formatted input file, e.g.
by doing

`$CTDYN_DIR/bin/ctdyn test/inputs/inlist_default`

CTDYN can also be run directly within a Python session using the API provided
by *py_ctdyn*, an example is provided in the `quickstart.ipynb` notebook.
 
## Contributors

- Sylvain N. Breton
- Alfio Bonanno 
- Giovanni Licciardello 

## Useful references

- [An accurate numerical approach for the kinematic dynamo problem](https://ui.adsabs.harvard.edu/abs/2004MSAIS...4...17B/abstract),
Bonanno 2004.
- [A solar mean field dynamo benchmark](https://ui.adsabs.harvard.edu/abs/2008A%26A...483..949J/abstract), Jouve et al. 2008.

