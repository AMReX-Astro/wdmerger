#!/bin/bash

source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

    if [ -z $ncell ]; then
	echoerr "ncell not set; exiting."
	exit
    fi

    to_run=1

    # Set the processor count and walltime.

    if [[ "$MACHINE" == "PERLMUTTER" ]]; then

        if [ $stellar_refinement -eq 1 ]; then
            nprocs=4
        elif [ $stellar_refinement -eq 2 ]; then
            nprocs=4
        elif [ $stellar_refinement -eq 4 ]; then
            nprocs=4
        elif [ $stellar_refinement -eq 8 ]; then
            nprocs=16
        elif [ $stellar_refinement -eq 16 ]; then
            nprocs=32
        elif [ $stellar_refinement -eq 32 ]; then
            nprocs=64
        elif [ $stellar_refinement -eq 64 ]; then
            nprocs=128
        fi

        walltime="2:00:00"

        amr__blocking_factor="32"
        amr__max_grid_size="1024"

    elif [[ "$MACHINE" == "FRONTIER" ]]; then

        if [ $stellar_refinement -eq 1 ]; then
            nprocs=8
        elif [ $stellar_refinement -eq 2 ]; then
            nprocs=8
        elif [ $stellar_refinement -eq 4 ]; then
            nprocs=8
        elif [ $stellar_refinement -eq 8 ]; then
            nprocs=16
        elif [ $stellar_refinement -eq 16 ]; then
            nprocs=32
        elif [ $stellar_refinement -eq 32 ]; then
            nprocs=64
        elif [ $stellar_refinement -eq 64 ]; then
            nprocs=128
        fi

        walltime="2:00:00"

        amr__blocking_factor="32"
        amr__max_grid_size="1024"

    else

        echoerr "This machine is not set up for this job."

    fi

    # Set up the geometry appropriately.

    amr__n_cell="$ncell $ncell"
    geometry__prob_lo="$prob_lo $prob_lo"
    geometry__prob_hi="$prob_hi $prob_hi"

    castro__lo_bc="3 3"

    # Set up refinement.

    amr__refinement_indicators="density"
    amr__refine__density__field_name="density"
    amr__refine__density__value_greater="1.e3"

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

# Needed for the makefile: we want to compile in 2D for these tests.

DIM="2"

# Use the aprox13 network for all following tests.

NETWORK_DIR="aprox13"

# Get the right inputs file.

inputs="inputs"

# Abort if we run out of GPU memory.

amrex__abort_on_out_of_gpu_memory="1"

# Limit GPU memory footprint.

castro__hydro_memory_footprint_ratio="3.0"

# Make it easier to see the real memory footprint.

amrex__the_arena_init_size="0"

# Disable flux limiting.

castro__limit_fluxes_on_small_dens="0"

# Variables we need to set up the collision.

problem__problem="0"
problem__collision_separation="2.0"
problem__collision_impact_parameter="0.0"

# Disable AMR subcycling.

amr__subcycling_mode="None"

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

castro__init_shrink="0.1"

# CFL number.

castro__cfl="0.8"

# Set a maximum dt to ensure that we don't jump too far ahead of the plotfile interval.

castro__max_dt="0.0025"

# Maximum number of subcycles.

castro__max_subcycles="128"

# Enable efficient regridding (don't actually regrid if the grids haven't changed.)

amr__use_efficient_regrid="1"

# Don't regrid on the fine levels.

amr__regrid_int="4 -1 -1 -1 -1 -1 -1 -1 -1 -1"

amr__error_buf="4 0 0 0 0 0 0 0 0 0"

# Enable reactions.

castro__do_react="1"

# Disable rotation.

castro__do_rotation="0"

# Ease up on the gravity tolerance since we're in axisymmetric and at high resolution.

gravity__abs_tol="5.e-8"

# We do not need to calculate the gravitational wave signature for this science.

castro__gw_dist="-1.e0"

# Many of the collision papers in the literature use an equal
# C/O ratio  by mass in the initial white dwarfs. We will do
# this too for comparison purposes.

problem__co_wd_c_frac="0.5e0"
problem__co_wd_o_frac="0.5e0"

# Allow the timestep to change by up to 25% per advance.

castro__change_max="1.25"

# Some defaults.

castro__dtnuc_e="1.e200"
castro__dtnuc_X="1.e200"
castro__react_T_min="1.0e8"
castro__react_rho_min="1.0e6"

# The flag we will use to determine whether to run the job.

to_run=1

# Simulate 0.64 + 0.64 M_solar WDs.

problem__mass_P="0.64"
problem__mass_S="0.64"

# Run four helium shell mass cases. These are the same
# values from Holcomb and Kushnir (2016), so that we
# can make a direct comparison (see Figure 2).

helium_shell_mass_list="0.00 0.04 0.08 0.16"
helium_shell_mass_list="0.00 0.16"

stop_time="4.0"
prob_lo="0.00e9"
prob_hi="5.12e9"

# Start with a base resolution of 400 km, and refine from there.

ncell="2048"
stellar_refinement_list="1 2 4 8 16 32 64"

for helium_shell_mass in $helium_shell_mass_list
do

    co_wd_he_shell_mass=$helium_shell_mass

    for stellar_refinement in $stellar_refinement_list
    do

        dir=$results_dir/he$helium_shell_mass/r$stellar_refinement
        set_run_opts

        if [ $to_run -eq 1 ]; then
            run
        fi

    done

done
