from importlib.metadata import version

__version__ = version ('py_ctdyn')

import py_ctdyn.templates 

from .ctdyn import (get_ctdyn_dir,
                    run_ctdyn,
                    set_default_inlist,
                    make_inlist,
                    load_inlist_template,
                    parse_input_file)

from .outputs import (read_summary_file,
                      read_field_map, 
                      read_field_map_text_file,
                      read_butterfly_diagram,
                      read_butterfly_diagram_text_file,
                      read_radial_profiles,
                      plot_meridional_map,
                      plot_butterfly_diagram,
                      plot_alpha, plot_eta)

