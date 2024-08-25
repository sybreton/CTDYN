import os, pytest
import py_ctdyn as dyn

@pytest.fixture(scope="session")
def tmp_dir (tmp_path_factory) :
  tmp = tmp_path_factory.mktemp ("test_outputs")
  return tmp

class TestUnitary :

  def testGetDir (self) :
    ctdyn_dir = dyn.get_ctdyn_dir ()
    assert type(ctdyn_dir)==str

  def testInlistDefault (self) :
    parameters = dyn.set_default_inlist ()
    assert type(parameters)==dict

  def testLoadInlist (self) :
    template = dyn.make_inlist ()

  def testMakeInlist (self) :
    template = dyn.make_inlist ()

  def testInlistTuned (self) :
    dir_out = "ctdyn_output"
    parameters = {"outputs" : {"dir":"'{}'".format (dir_out)},
                   }
    parameters = dyn.set_default_inlist (parameters=parameters)
    assert parameters["outputs"]["dir"] == "'{}'".format (dir_out)

class TestExecutionVisualisation :

  def testExecution (self, tmp_dir) :
    ctdyn_param = {"outputs" : {"dir":"'{}'".format (tmp_dir)},
               }
    dyn.run_ctdyn (ctdyn_param=ctdyn_param, verbose=False)

  def testPlotMeridionalMap (self, tmp_dir) : 
    ii, time = 1, 1
    filename = "{}/pfld.{}.t{}.A00".format (tmp_dir, str (ii).zfill (6), 
                                            str (time).zfill (2))
    r, theta, mesh = dyn.read_field_map (filename)
    fig = dyn.plot_meridional_map (r, theta, mesh, label=r"$B_r$")

  def testPlotButterflyDiagram (self, tmp_dir) :
    filename = "{}/butf.000001.a00".format (tmp_dir)
    t, theta, mesh = dyn.read_butterfly_diagram_text_file (filename)
    fig = dyn.plot_butterfly_diagram (t, theta, mesh)
