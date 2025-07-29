"""
Unit tests for pf_spinup_util module.
"""

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))
import pf_spinup_util


def test_create_model():
    """
    Test that the happy_path function returns the expected path.
    """

    runname = "pf_spinup"
    model = pf_spinup_util.get_pf_model(".")

    assert os.path.exists(os.path.join(".", runname, f"{runname}.yaml"))

    pf_spinup_util.save_model(model)
    model = pf_spinup_util.get_pf_model(".")


def test_initialize_model():
    """
    Test that the initialize_model function initializes the model correctly.
    """

    runname = "pf_spinup"
    directory_path = "."
    model_path = f"{directory_path}/{runname}"
    model = pf_spinup_util.get_pf_model(directory_path)

    options = {"grid_bounds": [4020, 1964, 4022, 1967], "grid": "conus2"}
    pf_spinup_util.initialize_model(model, options, model_path)


    # Check if the model has been initialized with default settings
    assert hasattr(model, "FileVersion")
    assert model.FileVersion == 4
    assert model.ComputationalGrid.DX == 1000.0
    pf_spinup_util.save_model(model)

    pf_spinup_util.run_model(model)

