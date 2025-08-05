import os
import numpy as np
import parflow
import shutil
import parflow.tools.fs
import parflow.tools.settings
import hf_hydrodata as hf
import subsettools as st


def create_model(
    runname: str,
    options: dict,
    directory_path: str,
    template_path="../src/conus2_transient_solid.yaml",
):
    directory_path = os.path.abspath(directory_path)

    os.makedirs(directory_path, exist_ok=True)
    template_path = os.path.abspath(template_path)
    parflow.tools.settings.set_working_directory(directory_path)
    model_path_ext = os.path.abspath(f"{directory_path}/{runname}.yaml")

    # Create Parfow run object
    shutil.copy(template_path, model_path_ext)
    model = parflow.Run.from_definition(model_path_ext)
    model.write(file_format="yaml")

    if options.get("grid_bounds") is None:
        raise ValueError("Grid bounds must be provided in options dictionary.")
    if options.get("grid") is None:
        raise ValueError("Grid must be provided in options dictionary.")
    grid = options["grid"]
    start_time = options.get("start_time", "2001-01-01")
    end_time = options.get("end_time", "2001-01-02")

    if grid not in ["conus1", "conus2"]:
        raise ValueError(
            f"Grid '{grid}' is not supported. Supported grids are 'conus1' and 'conus2'."
        )

    P = 1
    Q = 1
    R = 1
    model.Process.Topology.P = P
    model.Process.Topology.Q = Q
    model.Process.Topology.R = R
    model.FileVersion = 4

    # Locate the origin in the domain.
    grid_bounds = options.get("grid_bounds", None)
    latlon_bounds = options.get("latlon_bounds", None)
    if grid_bounds:
        lat_min, lon_min = hf.to_latlon(grid, grid_bounds[0], grid_bounds[1])
        lat_max, lon_max = hf.to_latlon(grid, grid_bounds[2] - 1, grid_bounds[3] - 1)
        latlon_bounds = [[lat_min, lon_min], [lat_max, lon_max]]
    ij_bounds, mask = st.define_latlon_domain(latlon_bounds, grid)

    model.ComputationalGrid.Lower.X = ij_bounds[0]
    model.ComputationalGrid.Lower.Y = ij_bounds[1]
    model.ComputationalGrid.Lower.Z = 0.0

    # Define the size of each grid cell. The length units are the same as those on hydraulic conductivity, here that is meters.
    model.ComputationalGrid.DX = 1000.0
    model.ComputationalGrid.DY = 1000.0
    model.ComputationalGrid.DZ = 200.0

    # Define the number of grid blocks in the domain.
    model.ComputationalGrid.NX = ij_bounds[2] - ij_bounds[0]
    model.ComputationalGrid.NY = ij_bounds[3] - ij_bounds[1]
    if grid == "conus1":
        model.ComputationalGrid.NZ = 5
    elif grid == "conus2":
        model.ComputationalGrid.NZ = 10
    mask_solid_paths = st.write_mask_solid(
        mask=mask, grid=grid, write_dir=directory_path
    )

    model.write(file_format="yaml")
    model = parflow.Run.from_definition(model_path_ext)

    var_ds = "conus2_domain"
    static_paths = st.subset_static(ij_bounds, dataset=var_ds, write_dir=directory_path)
    clm_paths = st.config_clm(
        ij_bounds,
        start=start_time,
        end=end_time,
        dataset=var_ds,
        write_dir=directory_path,
    )

    forcing_ds = "CW3E"
    st.subset_forcing(
        ij_bounds,
        grid=grid,
        start=start_time,
        end=end_time,
        dataset=forcing_ds,
        write_dir=directory_path,
    )

    st.edit_runscript_for_subset(
        ij_bounds,
        runscript_path=model_path_ext,
        runname=runname,
        forcing_dir=directory_path,
    )
    model = parflow.Run.from_definition(model_path_ext)

    init_press_path = os.path.basename(static_paths["ss_pressure_head"])
    depth_to_bedrock_path = os.path.basename(static_paths["pf_flowbarrier"])

    st.change_filename_values(
        runscript_path=model_path_ext,
        init_press=init_press_path,
        depth_to_bedrock=depth_to_bedrock_path,
    )

    model = parflow.Run.from_definition(model_path_ext)

    for fname in os.listdir(directory_path):
        if fname.endswith(".pfb"):
            data = parflow.read_pfb(f"{directory_path}/{fname}")
            parflow.write_pfb(f"{directory_path}/{fname}", data, dist=True)

    st.dist_run(
        topo_p=P,
        topo_q=Q,
        runscript_path=model_path_ext,
        dist_clim_forcing=True,
    )
    model = parflow.Run.from_definition(model_path_ext)
    model.TimingInfo.StopTime = 24
    model.Solver.CLM.MetFileName = "CW3E"
    model.write(file_format="yaml")

    return model
