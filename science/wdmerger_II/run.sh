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

        if [ $mass_refinement == "1.0e200" ]; then
            nprocs=18
        elif [ $mass_refinement == "1.0e28" ]; then
            nprocs=36
        elif [ $mass_refinement == "1.0e27" ]; then
            nprocs=72
        elif [ $mass_refinement == "1.0e26" ]; then
            nprocs=108
        fi

        walltime="2:00:00"

        amr__blocking_factor="32"
        amr__max_grid_size="64"

    else

        echoerr "This machine is not set up for this job."

    fi

    # Set up the geometry appropriately.

    amr__n_cell="$ncell $ncell $ncell"
    geometry__prob_lo="$prob_lo $prob_lo $prob_lo"
    geometry__prob_hi="$prob_hi $prob_hi $prob_hi"

    # Set up refinement.

    amr__refinement_indicators="mass tste"

    amr__refine__mass__field_name="density"
    amr__refine__mass__value_greater=$mass_refinement
    amr__refine__mass__volume_weighting="1"

    amr__refine__tste__field_name="t_sound_t_enuc"
    amr__refine__tste__value_greater=$tste_refinement

    amr__max_level="12"
    amr__ref_ratio="2"

}

# Specify the problem directory.

exec_dir=$CASTRO_HOME/Exec/science/wdmerger

use_first_castro_ex="1"

# Needed for the makefile: we want to compile in 3D for these tests.

DIM="3"

# Use the aprox13 network for all following tests.

NETWORK_DIR="aprox13"

# Get the right inputs file.

inputs="inputs"

# Variables we need to set up the merger.

problem__problem="1"
problem__roche_radius_factor="1.00"

# Limit GPU memory footprint.

castro__hydro_memory_footprint_ratio="3.0"

# Make it easier to see the real memory footprint.

amrex__the_arena_init_size="0"

# Full plotfiles.

amr__plot_per="-1.0"
amr__plot_vars="ALL"
amr__derive_plot_vars="ALL"

# Small plotfiles.

amr__small_plot_per="1.0"
amr__small_plot_vars="density Temp rho_e rho_c12 rho_o16 rho_si28 rho_ni56 enuc"
amr__derive_small_plot_vars="pressure soundspeed x_velocity y_velocity t_sound_t_enuc"

castro__plot_per_is_exact="0"
castro__small_plot_per_is_exact="0"

# Only plot level 0 in the plotfiles to save space.

amr__plot_max_level="0"
amr__small_plot_max_level="0"

# Checkpoints.

amr__check_per="10.0"

# Initial timestep shortening factor.

castro__init_shrink="0.01"

# CFL number.

castro__cfl="0.8"

# Maximum number of subcycles for retries.

castro__max_subcycles="128"

# Don't do subcycling at the AMR level.

amr__subcycling_mode="None"

# Use simplified SDC timestepper.

castro__time_integration_method="3"

# Since there will be explosive burns, add another iteration to SDC.

castro__sdc_iters="3"

# Burning timestep limiter is not as important when using SDC.

castro__dtnuc_e="1.0e200"
castro__dtnuc_X="1.0e200"

# Enable efficient regridding (don't actually regrid if the grids haven't changed.)

amr__use_efficient_regrid="1"

# Enable reactions.

castro__do_react="1"

# Limit maximum number of reaction integration steps.

integrator__ode_max_steps="15000"

# Set ambient density and temperature.

castro__ambient_density="1.e-4"
castro__ambient_temp="1.e6"

# Clamp ambient temperature for all densities below 1 g/cc.

castro__ambient_safety_factor="1.e4"

# We do not need to calculate the gravitational wave signature for this science.

castro__gw_dist="-1.e0"

# Add a helium shell.

problem__co_wd_he_shell_mass="0.01e0"

# Disable artificial viscosity since it causes problems with sharp composition gradients.

castro__difmag="0.0"

# Allow the timestep to change by up to 25% per advance.

castro__change_max="1.25"

# Set when burning turns on.

castro__react_T_min="1.0e8"
castro__react_rho_min="1.0e6"

# Simulate 1.1 + 0.9 M_solar WDs.

problem__mass_P="1.10"
problem__mass_S="0.90"

# Stop final stopping time.

stop_time="200.0"

# Set geometry.

prob_lo="-5.12e9"
prob_hi="5.12e9"

# Start with a base resolution of 400 km, and refine from there.

ncell="256"
mass_refinement_list="1.0e31 1.0e30 1.0e29 1.0e28 1.0e27 1.0e26"

# The flag we will use to determine whether to run the job.

to_run=1

for mass_refinement in $mass_refinement_list
do
    base_dir=$results_dir/m_$mass_refinement
    start_dir=$base_dir/start

    start_done="0"

    if [ -d $start_dir ]; then
        dir=$start_dir
        start_done=$(is_dir_done)
    fi

    if [ $start_done -ne 1 ]; then

        # First, run up until significant burning starts.

        castro__stopping_criterion_field="t_sound_t_enuc"
        castro__stopping_criterion_value="0.01"

        tste_refinement="0.1"

        dir=$start_dir
        set_run_opts

        if [ $to_run -eq 1 ]; then
            run
        fi

        unset castro__stopping_criterion_field
        unset castro__stopping_criterion_value
        unset tste_refinement

    else

        # Now, run with varying resolutions on t_sound / t_enuc.

        tste_refinement_list="1.0e200 1.0e-1"

        for tste_refinement in $tste_refinement_list
        do
            dir=$base_dir/t_$tste_refinement
            set_run_opts

            copy_checkpoint

            if [ $to_run -eq 1 ]; then
                run
            fi
        done

        unset tste_refinement

    fi

done
