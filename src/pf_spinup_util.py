import os
import numpy as np
import parflow
import shutil
import parflow.tools.fs
import parflow.tools.settings
import hf_hydrodata as hf


def get_pf_model(directory_path=".", runname="pf_spinup"):
    """
    Get or create a parflow model to use to execute a parflow run.

    Parameters:
    runname (str): The name of the run, used to create the directory.

    directory_path (str): The base path where the directory will be stored.

    Returns:
        The parflow model object
    """

    # Create the directory path
    directory_path = os.path.join(directory_path, runname)
    model_path = os.path.join(directory_path, runname)
    if os.path.exists(f"{model_path}.yaml"):
        model = parflow.Run.from_definition(f"{model_path}.yaml")
    else:
        # Create the directory if it does not exist
        os.makedirs(directory_path, exist_ok=True)

        # Create Parfow run object
        model = parflow.Run(runname, directory_path)
        model.FileVersion = 4
        model.write(file_format="yaml")

    return model


def save_model(model):
    """
    Save the Parflow model to a file.

    Parameters:
    model (parflow.Run): The Parflow model to save.

    Returns:
    None
    """

    # Save the model to a YAML file
    model.write(file_format="yaml")


def initialize_model(model, options: dict, model_path:str):
    intialize_slope(model, options, model_path)
    initialize_geometries(model, options, model_path)

    initialize_topology(model, options)
    initialize_porosity(model, options)
    initialize_permeability(model, options)
    initialize_saturation(model, options)
    initialize_mannings(model, options)
    initialize_wells(model, options)
    initialize_timing(model, options)
    initialize_evapotranspiration(model, options)
    initialize_pressure(model, options)
    initialize_solver(model, options)

def initialize_topology(model, options: dict):
    """
    Initialize the Parflow model with default settings.

    Parameters:
    model (parflow.Run): The Parflow model to initialize.

    Returns:
    None
    """
    if options.get("grid_bounds") is None:
        raise ValueError("Grid bounds must be provided in options dictionary.")
    if options.get("grid") is None:
        raise ValueError("Grid must be provided in options dictionary.")
    grid = options["grid"]
    if grid not in ["conus1", "conus2"]:
        raise ValueError(
            f"Grid '{grid}' is not supported. Supported grids are 'conus1' and 'conus2'."
        )

    # Processor topology: This is the way that the problem will be split across processors if you want to run in parallel
    # The domain is divided in x,y and z dimensions by P, Q and R. The total number of processors is P*Q*R.
    model.Process.Topology.P = 1
    model.Process.Topology.Q = 1
    model.Process.Topology.R = 1



    # Locate the origin in the domain.
    grid_bounds = options["grid_bounds"]
    model.ComputationalGrid.Lower.X = grid_bounds[0]
    model.ComputationalGrid.Lower.Y = grid_bounds[1]
    model.ComputationalGrid.Lower.Z = 0.0

    # Define the size of each grid cell. The length units are the same as those on hydraulic conductivity, here that is meters.
    model.ComputationalGrid.DX = 1000.0
    model.ComputationalGrid.DY = 1000.0
    model.ComputationalGrid.DZ = 200.0

    # Define the number of grid blocks in the domain.
    model.ComputationalGrid.NX = grid_bounds[2] - grid_bounds[0]
    model.ComputationalGrid.NY = grid_bounds[3] - grid_bounds[1]
    if grid == "conus1":
        model.ComputationalGrid.NZ = 5
    elif grid == "conus2": 
        model.ComputationalGrid.NZ = 10

def initialize_geometries(model, options:dict, model_path:str):

    #Declare the geometries that you will use for the problem
    model.GeomInput.Names = "indi_input boxinput"

    #Define the solid_input geometry.  
    #Note the naming convention here GeomInput.{GeomName}.key
    model.GeomInput.boxinput.InputType = "Box"
    model.GeomInput.boxinput.GeomName = "domain"
    model.Geom.domain.Lower.X = -1.0
    model.Geom.domain.Lower.Y = -1.0
    model.Geom.domain.Lower.Z = -1.0
    model.Geom.domain.Upper.X = options["grid_bounds"][2] - options["grid_bounds"][0] + 1.0
    model.Geom.domain.Upper.Y = options["grid_bounds"][3] - options["grid_bounds"][1]
    model.Geom.domain.Upper.Z = 11.0


    #First set the name for your `Domain` and setup the patches for this domain
    #model.Domain.GeomName = "domain"
    #model.Geom.domain.Patches = "top bottom side"

    # Next setup the indicator file geometry
    model.GeomInput.indi_input.InputType =   "IndicatorField"
    if options.get("grid") == "conus2":

        #filter_options = {"variable": "slope_x", "grid":options["grid"], "grid_bounds":options["grid_bounds"]}
        #data = hf.get_gridded_data(filter_options)
        grid_bounds = options["grid_bounds"]
        data = np.zeros((grid_bounds[3] - grid_bounds[1], grid_bounds[2] - grid_bounds[0], 10))
        data[:, :, :]= 1000.0
        print("**** SHAPE", data.shape)
        parflow.write_pfb(f"{model_path}/bedrock.pfb", data)

        if True:
            model.GeomInput.indi_input.GeomNames = "s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 g1 g2 g3 g4 g5 g6 g7 g8 g9 g10"
            model.Geom.indi_input.FileName = "bedrock.pfb"

            model.GeomInput.s1.Value =    1
            model.GeomInput.s2.Value =    2
            model.GeomInput.s3.Value =    3
            model.GeomInput.s4.Value =    4
            model.GeomInput.s5.Value =    5
            model.GeomInput.s6.Value =    6
            model.GeomInput.s7.Value =    7
            model.GeomInput.s8.Value =    8
            model.GeomInput.s9.Value =    9
            model.GeomInput.s10.Value =   10
            model.GeomInput.s11.Value =   11
            model.GeomInput.s12.Value =   12
            model.GeomInput.s13.Value =   13

            model.GeomInput.g1.Value =    19
            model.GeomInput.g2.Value =    20
            model.GeomInput.g3.Value =    21
            model.GeomInput.g4.Value =    22
            model.GeomInput.g5.Value =    23
            model.GeomInput.g6.Value =    24
            model.GeomInput.g7.Value =    25
            model.GeomInput.g8.Value =    26
            model.GeomInput.g9.Value =    27
            model.GeomInput.g10.Value =    28

        model.Solver.Nonlinear.VariableDz = True
        model.dzScale.GeomNames = "domain"
        model.dzScale.Type = "nzList"
        model.dzScale.nzListNumber = 10

        model.Cell._0.dzScale.Value = 5
        model.Cell._1.dzScale.Value = 0.5
        model.Cell._2.dzScale.Value = 0.25
        model.Cell._3.dzScale.Value = 0.125
        model.Cell._4.dzScale.Value = 0.05
        model.Cell._5.dzScale.Value = 0.025
        model.Cell._6.dzScale.Value = 0.005
        model.Cell._7.dzScale.Value = 0.003
        model.Cell._8.dzScale.Value = 0.0015
        model.Cell._9.dzScale.Value = 0.0005

def intialize_slope(model, options:dict, model_path:str):
    filter_options = {"variable": "slope_x", "grid":options["grid"], "grid_bounds":options["grid_bounds"]}
    data = hf.get_gridded_data(filter_options)
    parflow.write_pfb(f"{model_path}/slope_x.pfb", data)
    filter_options["variable"] = "slope_y"
    data = hf.get_gridded_data(filter_options)
    parflow.write_pfb(f"{model_path}/slope_y.pfb", data)

    model.TopoSlopesX.Type = "PFBFile"
    model.TopoSlopesX.GeomNames = "domain"
    model.TopoSlopesX.FileName = "slope_x.pfb"

    model.TopoSlopesY.Type = "PFBFile"
    model.TopoSlopesY.GeomNames = "domain"
    model.TopoSlopesY.FileName = "slope_y.pfb"

def initialize_permeability(model, options:dict):
    model.Geom.Perm.Names = "domain s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 g1 g2 g3 g4 g5 g6 g7 g8 g9 g10"

    if options["grid"] == "conus2":
        model.Geom.domain.Perm.Type = "Constant"
        model.Geom.domain.Perm.Value = 0.02

        model.Geom.s1.Perm.Type = "Constant"
        model.Geom.s1.Perm.Value = 0.269022595

        model.Geom.s2.Perm.Type = "Constant"
        model.Geom.s2.Perm.Value = 0.043630356

        model.Geom.s3.Perm.Type = "Constant"
        model.Geom.s3.Perm.Value = 0.015841225

        model.Geom.s4.Perm.Type = "Constant"
        model.Geom.s4.Perm.Value = 0.007582087

        model.Geom.s5.Perm.Type = "Constant"
        model.Geom.s5.Perm.Value = 0.01818816

        model.Geom.s6.Perm.Type = "Constant"
        model.Geom.s6.Perm.Value = 0.005009435

        model.Geom.s7.Perm.Type = "Constant"
        model.Geom.s7.Perm.Value = 0.005492736

        model.Geom.s8.Perm.Type = "Constant"
        model.Geom.s8.Perm.Value = 0.004675077

        model.Geom.s9.Perm.Type = "Constant"
        model.Geom.s9.Perm.Value = 0.003386794

        model.Geom.s10.Perm.Type = "Constant"
        model.Geom.s10.Perm.Value = 0.004783973

        model.Geom.s11.Perm.Type = "Constant"
        model.Geom.s11.Perm.Value = 0.003979136

        model.Geom.s12.Perm.Type = "Constant"
        model.Geom.s12.Perm.Value = 0.006162952

        model.Geom.s13.Perm.Type = "Constant"
        model.Geom.s13.Perm.Value = 0.005009435

        model.Geom.g1.Perm.Type = "Constant"
        model.Geom.g1.Perm.Value = 5e-3

        model.Geom.g2.Perm.Type = "Constant"
        model.Geom.g2.Perm.Value = 1e-2

        model.Geom.g3.Perm.Type = "Constant"
        model.Geom.g3.Perm.Value = 2e-2

        model.Geom.g4.Perm.Type = "Constant"
        model.Geom.g4.Perm.Value = 3e-2

        model.Geom.g5.Perm.Type = "Constant"
        model.Geom.g5.Perm.Value = 4e-2

        model.Geom.g6.Perm.Type = "Constant"
        model.Geom.g6.Perm.Value = 5e-2

        model.Geom.g7.Perm.Type = "Constant"
        model.Geom.g7.Perm.Value = 6e-2

        model.Geom.g8.Perm.Type = "Constant"
        model.Geom.g8.Perm.Value = 8e-2

        model.Geom.g9.Perm.Type = "Constant"
        model.Geom.g9.Perm.Value = 0.1

        model.Geom.g10.Perm.Type = "Constant"
        model.Geom.g10.Perm.Value = 0.2

    model.Perm.TensorType = "TensorByGeom"
    model.Geom.Perm.TensorByGeom.Names = "domain"
    model.Geom.domain.Perm.TensorValX = 1.0
    model.Geom.domain.Perm.TensorValY = 1.0
    model.Geom.domain.Perm.TensorValZ = 1.0

    model.SpecificStorage.Type = "Constant"
    model.SpecificStorage.GeomNames = "domain"
    model.Geom.domain.SpecificStorage.Value = 0.0001

def initialize_porosity(model, options:dict):
    if options["grid"] == "conus2":
        model.Geom.Porosity.GeomNames = "domain s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13"

        model.Geom.domain.Porosity.Type = "Constant"
        model.Geom.domain.Porosity.Value = 0.33

        model.Geom.s1.Porosity.Type = "Constant"
        model.Geom.s1.Porosity.Value = 0.375

        model.Geom.s2.Porosity.Type = "Constant"
        model.Geom.s2.Porosity.Value = 0.39

        model.Geom.s3.Porosity.Type = "Constant"
        model.Geom.s3.Porosity.Value = 0.387

        model.Geom.s4.Porosity.Type = "Constant"
        model.Geom.s4.Porosity.Value = 0.439

        model.Geom.s5.Porosity.Type = "Constant"
        model.Geom.s5.Porosity.Value = 0.489

        model.Geom.s6.Porosity.Type = "Constant"
        model.Geom.s6.Porosity.Value = 0.399

        model.Geom.s7.Porosity.Type = "Constant"
        model.Geom.s7.Porosity.Value = 0.384

        model.Geom.s8.Porosity.Type = "Constant"
        model.Geom.s8.Porosity.Value = 0.482

        model.Geom.s9.Porosity.Type = "Constant"
        model.Geom.s9.Porosity.Value = 0.442

        model.Geom.s10.Porosity.Type = "Constant"
        model.Geom.s10.Porosity.Value = 0.385

        model.Geom.s11.Porosity.Type = "Constant"
        model.Geom.s11.Porosity.Value = 0.481

        model.Geom.s12.Porosity.Type = "Constant"
        model.Geom.s12.Porosity.Value = 0.459

        model.Geom.s13.Porosity.Type = "Constant"
        model.Geom.s13.Porosity.Value = 0.399

        model.Phase.RelPerm.Type =              "VanGenuchten"
        model.Phase.RelPerm.GeomNames =     "domain s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13"

        model.Geom.domain.RelPerm.Alpha =    1.0
        model.Geom.domain.RelPerm.N =        3.0

        model.Geom.s1.RelPerm.Alpha =        3.548
        model.Geom.s1.RelPerm.N =            4.162

        model.Geom.s2.RelPerm.Alpha =        3.467
        model.Geom.s2.RelPerm.N =            2.738

        model.Geom.s3.RelPerm.Alpha =        2.692
        model.Geom.s3.RelPerm.N =            2.445

        model.Geom.s4.RelPerm.Alpha =        0.501
        model.Geom.s4.RelPerm.N =            2.659

        model.Geom.s5.RelPerm.Alpha =        0.661
        model.Geom.s5.RelPerm.N =            2.659

        model.Geom.s6.RelPerm.Alpha =        1.122
        model.Geom.s6.RelPerm.N =            2.479

        model.Geom.s7.RelPerm.Alpha =        2.089
        model.Geom.s7.RelPerm.N =            2.318

        model.Geom.s8.RelPerm.Alpha =        0.832
        model.Geom.s8.RelPerm.N =            2.514

        model.Geom.s9.RelPerm.Alpha =        1.585
        model.Geom.s9.RelPerm.N =            2.413

        model.Geom.s10.RelPerm.Alpha =        3.311
        model.Geom.s10.RelPerm.N =            2.202

        model.Geom.s11.RelPerm.Alpha =        1.622
        model.Geom.s11.RelPerm.N =            2.318

        model.Geom.s12.RelPerm.Alpha =        1.514
        model.Geom.s12.RelPerm.N =            2.259

        model.Geom.s13.RelPerm.Alpha =        1.122
        model.Geom.s13.RelPerm.N =            2.479        

def initialize_saturation(model, options:dict):
    model.Phase.Saturation.Type =             "VanGenuchten"
    model.Phase.Saturation.GeomNames =         "domain s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13"

    model.Geom.domain.Saturation.Alpha =        1.0
    model.Geom.domain.Saturation.N =            3.0
    model.Geom.domain.Saturation.SRes =         0.001
    model.Geom.domain.Saturation.SSat =         1.0

    model.Geom.s1.Saturation.Alpha =        3.548
    model.Geom.s1.Saturation.N =            4.162
    model.Geom.s1.Saturation.SRes =         0.0001
    model.Geom.s1.Saturation.SSat =         1.0

    model.Geom.s2.Saturation.Alpha =        3.467
    model.Geom.s2.Saturation.N =            2.738
    model.Geom.s2.Saturation.SRes =         0.0001
    model.Geom.s2.Saturation.SSat =         1.0

    model.Geom.s3.Saturation.Alpha =        2.692
    model.Geom.s3.Saturation.N =            2.445
    model.Geom.s3.Saturation.SRes =         0.0001
    model.Geom.s3.Saturation.SSat =         1.0

    model.Geom.s4.Saturation.Alpha =        0.501
    model.Geom.s4.Saturation.N =            2.659
    model.Geom.s4.Saturation.SRes =         0.1
    model.Geom.s4.Saturation.SSat =         1.0

    model.Geom.s5.Saturation.Alpha =        0.661
    model.Geom.s5.Saturation.N =            2.659
    model.Geom.s5.Saturation.SRes =         0.0001
    model.Geom.s5.Saturation.SSat =         1.0

    model.Geom.s6.Saturation.Alpha =        1.122
    model.Geom.s6.Saturation.N =            2.479
    model.Geom.s6.Saturation.SRes =         0.0001
    model.Geom.s6.Saturation.SSat =         1.0

    model.Geom.s7.Saturation.Alpha =        2.089
    model.Geom.s7.Saturation.N =            2.318
    model.Geom.s7.Saturation.SRes =         0.0001
    model.Geom.s7.Saturation.SSat =         1.0

    model.Geom.s8.Saturation.Alpha =        0.832
    model.Geom.s8.Saturation.N =            2.514
    model.Geom.s8.Saturation.SRes =         0.0001
    model.Geom.s8.Saturation.SSat =         1.0

    model.Geom.s9.Saturation.Alpha =        1.585
    model.Geom.s9.Saturation.N =            2.413
    model.Geom.s9.Saturation.SRes =         0.0001
    model.Geom.s9.Saturation.SSat =         1.0

    model.Geom.s10.Saturation.Alpha =        3.311
    model.Geom.s10.Saturation.N =            2.202
    model.Geom.s10.Saturation.SRes =         0.0001
    model.Geom.s10.Saturation.SSat =         1.0

    model.Geom.s11.Saturation.Alpha =        1.622
    model.Geom.s11.Saturation.N =            2.318
    model.Geom.s11.Saturation.SRes =         0.0001
    model.Geom.s11.Saturation.SSat =         1.0

    model.Geom.s12.Saturation.Alpha =        1.514
    model.Geom.s12.Saturation.N =            2.259
    model.Geom.s12.Saturation.SRes =         0.0001
    model.Geom.s12.Saturation.SSat =         1.0

    model.Geom.s13.Saturation.Alpha =        1.122
    model.Geom.s13.Saturation.N =            2.479
    model.Geom.s13.Saturation.SRes =         0.0001
    model.Geom.s13.Saturation.SSat =         1.0

def initialize_mannings(model, options:dict):
    model.Mannings.Type = "Constant"
    model.Mannings.GeomNames = "domain"
    model.Mannings.Geom.domain.Value = 0.0000044

def initialize_wells(model, options:dict):
    # Phases
    model.Phase.Names = "water"
    model.Phase.water.Density.Type = "Constant"
    model.Phase.water.Density.Value = 1.0
    model.Phase.water.Viscosity.Type = "Constant"
    model.Phase.water.Viscosity.Value = 1.0
    model.Phase.water.Mobility.Type = "Constant"
    model.Phase.water.Mobility.Value = 1.0

    # Contaminants
    model.Contaminants.Names = ""

    # Gravity
    model.Gravity = 1.0

    #Wells
    model.Wells.Names = ""

    # Phase Sources
    model.PhaseSources.water.Type = "Constant"
    model.PhaseSources.water.GeomNames = "domain"
    model.PhaseSources.water.Geom.domain.Value = 0.0

def initialize_timing(model, options:dict):
    # Time units ans start and stop times
    model.TimingInfo.BaseUnit = 1.0
    model.TimingInfo.StartCount = 0
    model.TimingInfo.StartTime = 0.0
    model.TimingInfo.StopTime = 10000000.0
    model.TimingInfo.DumpInterval = 100.0

    # Growth timestep properties
    model.TimeStep.Type = "Growth"
    model.TimeStep.InitialStep = 1.0
    model.TimeStep.GrowthFactor = 1.1
    model.TimeStep.MaxStep = 100
    model.TimeStep.MinStep = 1

    #Time cycles
    model.Cycle.Names ="constant"
    model.Cycle.constant.Names = "alltime"
    model.Cycle.constant.alltime.Length = 1
    model.Cycle.constant.Repeat = -1

def initialize_evapotranspiration(model, option:dict):
    model.BCPressure.PatchNames = "top bottom side"

    model.Patch.bottom.BCPressure.Type = "FluxConst"
    model.Patch.bottom.BCPressure.Cycle = "constant"
    model.Patch.bottom.BCPressure.alltime.Value = 0.0

    model.Patch.side.BCPressure.Type = "FluxConst"
    model.Patch.side.BCPressure.Cycle = "constant"
    model.Patch.side.BCPressure.alltime.Value = 0.0

    model.Patch.top.BCPressure.Type = "OverlandKinematic"
    model.Patch.top.BCPressure.Cycle = "constant"
    model.Patch.top.BCPressure.alltime.Value = -2.1e-5  

def initialize_pressure(model, options:dict):
    # Starting from a constant head values
    model.ICPressure.Type = "HydroStaticPatch"
    model.ICPressure.GeomNames = "domain"
    model.Geom.domain.ICPressure.Value = 0.0
    model.Geom.domain.ICPressure.RefGeom = "domain"
    model.Geom.domain.ICPressure.RefPatch = "bottom"

def initialize_solver(model, options:dict):
    return
    model.Solver.PrintSubsurfData = True
    model.Solver.PrintPressure = True
    model.Solver.PrintSaturation = True
    model.Solver.PrintMask = True
    model.Solver.PrintVelocities = False

    # Solver types
    model.Solver = "Richards"
    model.Solver.TerrainFollowingGrid = True
    model.Solver.Linear.Preconditioner = "PFMG"
    model.Solver.Linear.Preconditioner.PCMatrixType = "FullJacobian"

    # Exact solution
    model.KnownSolution = "NoKnownSolution"

    # Solver settings
    model.Solver.MaxIter = 25000
    model.Solver.Drop = 1e-20
    model.Solver.AbsTol = 1e-8

    model.Solver.MaxConvergenceFailures = 8
    model.Solver.Nonlinear.MaxIter = 1000
    model.Solver.Nonlinear.ResidualTol = 1e-6
    model.Solver.Nonlinear.EtaChoice =  "EtaConstant"
    model.Solver.Nonlinear.EtaValue = 0.001
    model.Solver.Nonlinear.UseJacobian = True
    model.Solver.Nonlinear.DerivativeEpsilon = 1e-16
    model.Solver.Nonlinear.StepTol = 1e-15
    model.Solver.Nonlinear.Globalization = "LineSearch"


    model.Solver.Linear.KrylovDimension = 70
    model.Solver.Linear.MaxRestarts = 2

    #model.Solver.OverlandFlowSpinUp = 1

def run_model(model):
    """
    Run the Parflow model.

    Parameters:
    model (parflow.Run): The Parflow model to run.

    Returns:
    None
    """
    # Run the model
    model.run()

    # Check if the run was successful
    if not model.is_successful():
        raise RuntimeError("ParFlow run failed.")
    