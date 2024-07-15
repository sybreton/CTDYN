# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
from pkg_resources import DistributionNotFound, get_distribution

try:
    __version__ = get_distribution("py_ctdyn").version
except DistributionNotFound:
    __version__ = "unknown version"
sys.path.insert(0, os.path.abspath('.'))

autodoc_mock_imports = ['numpy', 'scipy', 'matplotlib']


project = 'CTDYN'
copyright = '2024, A. Bonanno, S.N. Breton'
author = 'A. Bonanno, S.N. Breton, G. Licciardello'
release = '0.1'
master_doc = 'index'
source_suffix = ".rst"

release = __version__
version = __version__


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.mathjax',
              'sphinx_book_theme',
              'IPython.sphinxext.ipython_console_highlighting']

napoleon_numpy_docstring = True

exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_title = "star-privateer"
html_logo = "logo.png"
html_theme_options = {
     "path_to_docs": "docs",
     "repository_url": "https://gitlab.com/sybreton/star_privateer/",
     "use_repository_button": True,
     "use_download_button": True,
     }
html_sidebars = {
    "**": [
     "navbar-logo.html",
     "search-field.html",
     "sbt-sidebar-nav.html",
          ]
     }
numpydoc_show_class_members = False

html_static_path = []
