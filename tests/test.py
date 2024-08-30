import os, pytest
import py_ctdyn as dyn

@pytest.fixture(scope="session")
def tmp_dir (tmp_path_factory) :
  tmp = tmp_path_factory.mktemp ("test_outputs")
  return tmp

class TestUnitary :

  @pytest.fixture(scope="class")
  def namelists (self) :
    return ["&global", "&profiles", "&brent",
            "&boundaries", "&fields", "&physics",
            "&outputs", "&controls"]

  def testGetDir (self) :
    ctdyn_dir = dyn.get_ctdyn_dir ()
    assert type(ctdyn_dir)==str

  def testInlistDefault (self) :
    parameters = dyn.set_default_inlist ()
    assert type(parameters)==dict

  def check_template (self, template,
                      namelists) :
    for namelist in namelists :
      missing = True
      ii = 0
      while missing and ii<len(template) :
        if namelist in template[ii] :
          missing = False
        else :
          ii += 1
      assert not missing 

  def testLoadInlist (self, namelists) :
    template = dyn.make_inlist ()
    self.check_template (template, namelists)

  def testMakeInlistFilenameNone (self, namelists) :
    template = dyn.make_inlist (filename=None)
    self.check_template (template, namelists)
  
  def testMakeInlistWithFilename (self, tmp_dir,
                                  namelists) :
    filename = os.path.join (tmp_dir, "inlist")
    template = dyn.make_inlist (filename=filename)
    self.check_template (template, namelists)
    assert os.path.exists (filename)

  def testInlistTuned (self, tmp_dir) :
    parameters = {"outputs" : {"dir":"'{}'".format (tmp_dir)},
                   }
    parameters = dyn.set_default_inlist (parameters=parameters)
    assert parameters["outputs"]["dir"] == "'{}'".format (tmp_dir)

class TestExecutionVisualisation :

  def testExecution (self, tmp_dir) :
    ctdyn_param = {"outputs" : {"dir":"'{}'".format (tmp_dir)},
               }
    dyn.run_ctdyn (ctdyn_param=ctdyn_param, verbose=False)

  def testRerun (self, tmp_dir) :
    ctdyn_param = {"outputs" : {"dir":"'{}'".format (tmp_dir)},
               }
    dyn.run_ctdyn (ctdyn_param=ctdyn_param, verbose=False,
                   rerun=False)

  def testRadialProfiles (self, tmp_dir) :
    filename = os.path.join (tmp_dir, "alpha.dat")
    df = dyn.read_radial_profiles (filename) 
    fig = dyn.plot_alpha (df)
    fig = dyn.plot_eta (df)

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
