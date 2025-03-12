.. CTDYN documentation master file, created by
   sphinx-quickstart on Fri Jul  5 11:58:58 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CTDYN: the Catania dynamo code
=================================

CTDYN provides a Fortran implementation of the stellar mean field dynamo problem
(see `Bonanno 2004 <https://ui.adsabs.harvard.edu/abs/2004MSAIS...4...17B/abstract>`_
and `Jouve et al. 2008 <https://ui.adsabs.harvard.edu/abs/2008A%26A...483..949J/abstract>`_).
It is released together with the auxiliary Python module `py_ctdyn` which offers
a dedicated interface for running CTDYN and analysing its output.

For now the source code is hosted on a private `GitHub repository <https://github.com/sybreton/CTDYN>`_. 

.. toctree::
   :maxdepth: 1
   :caption: User guide:

   usage/installation/installation
   usage/quickstart/quickstart
   usage/namelist/namelist_description
   usage/parallel_run/parallel_run

.. toctree::
   :maxdepth: 1
   :caption: Detailed API:

   usage/api/interface.rst
   usage/api/visualisation.rst


