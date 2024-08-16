# CTDYN

This repository is a fork of the Catania dynamo code project originally
developed by Alfio Bonanno. The code has been refactored in depth in
order to upgrade its structure and syntax towards modern Fortran standards.
A Python module dedicated to analysis and display of code outputs is also
included with this distribution.  

## Installation

### Fortran CTDYN executable

The first thing is to compile the source code

`cd src`

`make`

An executable should be created at the root of the repository under the name
`ctdyn`.  

### Python module

The *py_ctdyn* module provides a library of Python tools designed to analyse
and display the outputs of CTDYN computations. The module can be installed by
running the following command at the root of the repository. 

`pip install .`
 

## Get started

The executable can be run  by selecting any correctly formatted input file, e.g.
by doing

`./ctdyn test/inputs/test_ref`
 
## Contributors

- Sylvain N. Breton
- Alfio Bonanno 
- Giovanni Licciardello 
