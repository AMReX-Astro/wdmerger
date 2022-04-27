#!/bin/bash

source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

    if [ -z $ncell ]; then
	echoerr "ncell not set; exiting."
	exit
    fi

    to_run=1

    # Set the processor count and walltime.

    if [[ "$MACHINE" == "SUMMIT" ]]; then

        if [ $stellar_refinement -eq 1 ]; then
            nprocs=16
        elif [ $stellar_refinement -eq 2 ]; then
            nprocs=16
        elif [ $stellar_refinement -eq 4 ]; then
            nprocs=16
        fi

        walltime="2:00:00"

        amr__blocking_factor="16"
        amr__max_grid_size="128"

    else

        echoerr "This machine is not set up for this job."

    fi

    # Set up the geometry appropriately.

    amr__n_cell="$ncell $ncell $ncell"
    geometry__prob_lo="$prob_lo $prob_lo $prob_lo"
    geometry__prob_hi="$prob_hi $prob_hi $prob_hi"

    # Set up refinement.

    amr__refinement_indicators="density temperature"

    amr__refine__density__field_name="density"
    amr__refine__density__value_greater="1.e3"

    amr__refine__temperature__field_name="Temp"
    amr__refine__temperature__value_greater="1.e9"

    if [ $stellar_refinement -eq 1 ]; then

        amr__max_level="0"

    elif [ $stellar_refinement -eq 2 ]; then

        amr__max_level="1"
        amr__ref_ratio="2"

    elif [ $stellar_refinement -eq 4 ]; then

        amr__max_level="1"
        amr__ref_ratio="4"

    elif [ $stellar_refinement -eq 8 ]; then

        amr__max_level="2"
        amr__ref_ratio="4 2"

    elif [ $stellar_refinement -eq 16 ]; then

        amr__max_level="2"
        amr__ref_ratio="4 4"

    elif [ $stellar_refinement -eq 32 ]; then

        amr__max_level="3"
        amr__ref_ratio="4 4 2"

    elif [ $stellar_refinement -eq 64 ]; then

        amr__max_level="3"
        amr__ref_ratio="4 4 4"

    elif [ $stellar_refinement -eq 128 ]; then

        amr__max_level="4"
        amr__ref_ratio="4 4 4 2"

    elif [ $stellar_refinement -eq 256 ]; then

        amr__max_level="4"
        amr__ref_ratio="4 4 4 4"

    elif [ $stellar_refinement -eq 512 ]; then

        amr__max_level="5"
        amr__ref_ratio="4 4 4 4 2"

    elif [ $stellar_refinement -eq 1024 ]; then

        amr__max_level="5"
        amr__ref_ratio="4 4 4 4 4"

    elif [ $stellar_refinement -eq 2048 ]; then

        amr__max_level="6"
        amr__ref_ratio="4 4 4 4 4 2"

    elif [ $stellar_refinement -eq 4096 ]; then

        amr__max_level="6"
        amr__ref_ratio="4 4 4 4 4 4"

    fi
}

# Specify the problem directory.

exec_dir=$CASTRO_HOME/Exec/science/wdmerger

use_first_castro_ex="1"

# Needed for the makefile: we want to compile in 3D for these tests.

DIM="3"

# Use the aprox19 network for all following tests.

NETWORK_DIR="aprox19"

# Get the right inputs file.

inputs="inputs"

# Abort if we run out of GPU memory.

amrex__abort_on_out_of_gpu_memory="1"

# Variables we need to set up the merger.

problem__problem="1"
problem__roche_radius_factor="1.0"

# Allow first-order interpolations to fine levels.

castro__state_interp_order="1"

# Full plotfiles.

amr__plot_per="-1.0"
amr__plot_vars="ALL"
amr__derive_plot_vars="ALL"

# Small plotfiles.

amr__small_plot_per="0.01"
amr__small_plot_vars="density Temp rho_e rho_c12 rho_o16 rho_si28 rho_ni56 enuc"
amr__derive_small_plot_vars="pressure soundspeed x_velocity y_velocity t_sound_t_enuc"

castro__plot_per_is_exact="0"
castro__small_plot_per_is_exact="0"

# Checkpoints.

amr__check_per="0.5"

# Initial timestep shortening factor.

castro__init_shrink="0.01"

# CFL number.

castro__cfl="0.8"

# Maximum number of subcycles.

castro__max_subcycles="128"

# Enable efficient regridding (don't actually regrid if the grids haven't changed.)

amr__use_efficient_regrid="1"

# Set regridding buffer.

amr__regrid_int="4"
amr__error_buf="4"

# Enable reactions.

castro__do_react="1"

# Disable rotation.

castro__do_rotation="0"

# Set ambient density and temperature.

castro__ambient_density="1.e-4"
castro__ambient_temp="1.e6"

# Clamp ambient temperature for all densities below 1 g/cc.

castro__ambient_safety_factor="1.e4"

# We do not need to calculate the gravitational wave signature for this science.

castro__gw_dist="-1.e0"

# Set the C/O ratio.

problem__co_wd_c_frac="0.4e0"
problem__co_wd_o_frac="0.6e0"

# Add a helium shell.

problem__co_wd_he_shell_mass="0.01e0"

# Allow the timestep to change by up to 25% per advance.

castro__change_max="1.25"

# Some defaults.

castro__dtnuc_e="1.e200"
castro__dtnuc_X="1.e200"
castro__react_T_min="1.0e8"
castro__react_rho_min="1.0e6"

# The flag we will use to determine whether to run the job.

to_run=1

# Simulate 1.0 + 0.6 M_solar WDs.

problem__mass_P="1.00"
problem__mass_S="0.60"

# Stop final stopping time.

stop_time="100.0"

# Set geometry.

prob_lo="-7.0e9"
prob_hi="7.0e9"

# Start with a base resolution of 546.9 km, and refine from there.

ncell="256"
stellar_refinement_list="4"

for stellar_refinement in $stellar_refinement_list
do

    dir=$results_dir/r$stellar_refinement
    set_run_opts

    if [ $to_run -eq 1 ]; then
        run
    fi

done
