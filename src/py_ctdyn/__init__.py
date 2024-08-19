from importlib.metadata import version

__version__ = version ('py_ctdyn')

import py_ctdyn.templates 

from .ctdyn import (get_ctdyn_dir,
                    run_ctdyn,
                    set_default_inlist)

from .outputs import (read_field, 
                      read_field_text_file,
                      read_butterfly_diagram,
                      read_butterfly_diagram_text_file,
                      plot_meridional_mesh,
                      plot_butterfly_diagram)

