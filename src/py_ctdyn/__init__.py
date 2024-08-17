from importlib.metadata import version

__version__ = version ('py_ctdyn')

from .outputs import (read_field, 
                      read_field_text_file,
                      read_butterfly_diagram,
                      read_butterfly_diagram_text_file,
                      plot_meridional_mesh,
                      plot_butterfly_diagram)

