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

autodoc_mock_imports = ['numpy', 'scipy', 'matplotlib',
                        'pandas']


project = 'CTDYN'
copyright = '2024, A. Bonanno, S.N. Breton'
author = 'A. Bonanno, S.N. Breton, G. Licciardello'
master_doc = 'index'
source_suffix = ".rst"

release = __version__
version = __version__


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.mathjax',
              'IPython.sphinxext.ipython_console_highlighting']

napoleon_numpy_docstring = True

exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_title = "CTDYN"
html_logo = "logo.png"

json_url = "https://ctdyn.netlify.app/latest/_static/switcher.json"

html_context = {"github_repo":"https://github.com/sybreton/CTDYN"}
html_theme_options = dict ()
html_theme_options["navbar_start"] = ["navbar-logo", 
                                      "version-switcher"]
html_theme_options["navbar_end"] = ["navbar-icon-links"] 

current_version = __version__
# Dev branch set to be the latest 
if current_version=="dev" :
  current_version = "latest"
html_theme_options["switcher"] = {"json_url":json_url, 
                                  "version_match":current_version}

html_theme_options["check_switcher"] = False

html_sidebars = {
    "**": ["sidebar-nav-bs", "sidebar-ethical-ads"]
}

numpydoc_show_class_members = False
html_static_path = []
