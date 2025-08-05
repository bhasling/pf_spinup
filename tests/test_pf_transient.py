"""
Unit tests for pf_spinup_util module.
"""

import sys
import os
import datetime

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))
import pf_transient_util




def test_create_model():
  """
  Test that the initialize_model function initializes the model correctly.
  """

  try:
      runname = "pf_transient"
      directory_path = "./pf_transient"


      start_time = "2001-01-01"
      end_time = "2001-01-02"
      options = {"grid_bounds": [4020, 1964, 4022, 1967], "grid": "conus2", "start_time":start_time, "end_time": end_time}

      model = pf_transient_util.create_model(runname, options, directory_path)


      # Check if the model has been initialized with default settings
      assert hasattr(model, "FileVersion")
      assert model.FileVersion == 4
      assert model.ComputationalGrid.DX == 1000.0

      os.chdir(directory_path)
      model.run()
      return

    # Check if the run was successful
      if not model.is_successful():
          raise RuntimeError("ParFlow run failed.")
  except Exception as e:
    raise e
        
