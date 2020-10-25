#!/bin/bash

source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

    if [ -z $ncell ]; then
	echoerr "ncell not set; exiting."
	exit
    fi

    to_run=1

    # Set the processor count and walltime.

    if [ "$MACHINE" == "SUMMIT" ]; then

        queue="batch"

        if [ $ncell -le 1024 ]; then
            nprocs=48
        elif [ $ncell -le 2048 ]; then
            nprocs=96
        elif [ $ncell -le 4096 ]; then
            nprocs=192
        elif [ $ncell -le 8192 ]; then
            nprocs=384
        elif [ $ncell -le 16384 ]; then
            nprocs=768
        fi

        walltime="2:00:00"

        amr__blocking_factor="128"
        amr__max_grid_size="1024"

    else

        echoerr "This machine is not set up for this job."

    fi

    # Set up the geometry appropriately.

    amr__n_cell="$ncell $ncell"
    geometry__prob_lo="$prob_lo $prob_lo"
    geometry__prob_hi="$prob_hi $prob_hi"

    castro__lo_bc="3 3"

    if [ $hydro_refinement -eq 1 ]; then
        amr__max_level=0
    elif [ $hydro_refinement -eq 2 ]; then
        amr__max_level=1
        amr__ref_ratio="2"
    elif [ $hydro_refinement -eq 4 ]; then
        amr__max_level=1
        amr__ref_ratio="4"
    elif [ $hydro_refinement -eq 8 ]; then
        amr__max_level=2
        amr__ref_ratio="4 2"
    elif [ $hydro_refinement -eq 16 ]; then
        amr__max_level=2
        amr__ref_ratio="4 4"
    elif [ $hydro_refinement -eq 32 ]; then
        amr__max_level=3
        amr__ref_ratio="4 4 2"
    elif [ $hydro_refinement -eq 64 ]; then
        amr__max_level=3
        amr__ref_ratio="4 4 4"
    elif [ $hydro_refinement -eq 128 ]; then
        amr__max_level=4
        amr__ref_ratio="4 4 4 2"
    elif [ $hydro_refinement -eq 256 ]; then
        amr__max_level=4
        amr__ref_ratio="4 4 4 4"
    elif [ $hydro_refinement -eq 512 ]; then
        amr__max_level=5
        amr__ref_ratio="4 4 4 4 2"
    elif [ $hydro_refinement -eq 1024 ]; then
        amr__max_level=5
        amr__ref_ratio="4 4 4 4 4"
    elif [ $hydro_refinement -eq 2048 ]; then
        amr__max_level=6
        amr__ref_ratio="4 4 4 4 4 2"
    elif [ $hydro_refinement -eq 4096 ]; then
        amr__max_level=6
        amr__ref_ratio="4 4 4 4 4 4"
    elif [ $hydro_refinement -eq 8192 ]; then
        amr__max_level=7
        amr__ref_ratio="4 4 4 4 4 4 2"
    elif [ $hydro_refinement -eq 16384 ]; then
        amr__max_level=7
        amr__ref_ratio="4 4 4 4 4 4 4"
    elif [ $hydro_refinement -eq 32768 ]; then
        amr__max_level=8
        amr__ref_ratio="4 4 4 4 4 4 4 2"
    elif [ $hydro_refinement -eq 65536 ]; then
        amr__max_level=8
        amr__ref_ratio="4 4 4 4 4 4 4 4"
    elif [ $hydro_refinement -eq 131072 ]; then
        amr__max_level=9
        amr__ref_ratio="4 4 4 4 4 4 4 4 2"
    elif [ $hydro_refinement -eq 262144 ]; then
        amr__max_level=9
        amr__ref_ratio="4 4 4 4 4 4 4 4 4"
    elif [ $hydro_refinement -eq 524288 ]; then
        amr__max_level=10
        amr__ref_ratio="4 4 4 4 4 4 4 4 4 2"
    elif [ $hydro_refinement -eq 1048576 ]; then
        amr__max_level=10
        amr__ref_ratio="4 4 4 4 4 4 4 4 4 4"
    elif [ $hydro_refinement -eq 2097152 ]; then
        amr__max_level=11
        amr__ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
    else
        echo "Unknown refinement factor: "$hydro_refinement"; exiting."
        exit
    fi

    # Only solve gravity on the coarse level. The additional accuracy
    # from solving on the fine levels would provide no benefit on this problem.

    gravity__max_solve_level=0

    # Only solve reactions on the coarse level. On levels above this,
    # where we're effectively subcycling the hydro, we just interpolate
    # the reactions source rather than calculate it.

    castro__reactions_max_solve_level=0

    # The refinement strategy is to maximally refine above a given temperature.

    temp_min="1.0e9"

    amr__refinement_indicators="burning"
    amr__refine__burning__field_name="Temp"
    amr__refine__burning__max_level="$amr__max_level"
    amr__refine__burning__value_greater="$temp_min"

}

function copy_checkpoint() {

    # Copy the initial checkpoint to the directory if it's not there yet.

    if [ -e "$dir/jobIsDone" ]; then
        return
    fi

    job_flag=$(is_job_running $dir)

    if [ $job_flag -ne 1 ]; then

        done_flag=$(is_dir_done)

        if [ $done_flag -ne 1 ]; then

            no_submit=1

            if [ $to_run -eq 1 ]; then
                run
            fi

            unset no_submit

            # Skip the saved checkpoint copy if there's already a checkpoint
            # in the target directory (either it's the saved checkpoint, or it
            # is a later one because we've already started the run).

            checkpoint=$(get_last_checkpoint $dir)

            if [ "$checkpoint" == "" ]; then

                # By default the starting checkpoint is the last checkpoint
                # from the start directory, but this can be overridden if
                # the user prefers a different starting checkpoint.

                if [ -z $start_checkpoint ]; then
                    start_checkpoint=$(get_last_checkpoint $start_dir)
                fi

                if [ -z $checkpoint_dir ]; then
                    checkpoint_dir=$start_dir
                fi

                echo "Copying initial checkpoint to" $dir"."

                if [ -d "$checkpoint_dir/$start_checkpoint" ]; then

                    cp -r "$checkpoint_dir/$start_checkpoint" $dir
                    rm -f "$dir/$start_checkpoint/jobIsDone"

                    # Also, copy the diagnostic output files
                    # so there is an apparently continuous record.

                    cp $checkpoint_dir/*diag.out $dir

                    # Strip out any data in the diagnostics files
                    # that is after the copied checkpoint, to
                    # ensure apparent continuity.

                    chk_step=$(echo $start_checkpoint | cut -d"k" -f2)
                    chk_step=$(echo $chk_step | bc) # Strip leading zeros

                    for diag in $(find $dir -name "*diag.out");
                    do
                        sed -i "/\ $chk_step \ /q" $diag
                    done

                else
                    echoerr "No initial checkpoint available, exiting."
                    exit
                fi

            fi

        fi

    fi

    if [ ! -z $start_checkpoint ]; then
        unset start_checkpoint
    fi

    if [ ! -z $checkpoint_dir ]; then
        unset checkpoint_dir
    fi

}

# Specify the problem directory.

problem_dir=$CASTRO_HOME/Exec/science/wdmerger

use_first_castro_ex="1"

# Needed for the makefile: we want to compile in 2D for these tests.

DIM="2"

# Use the aprox13 network for all following tests.

NETWORK_DIR="aprox13"

# Get the right inputs and probin files.

inputs="inputs_2d"
probin="probin"

# Abort if we run out of GPU memory.

amrex__abort_on_out_of_gpu_memory="1"

# Disable flux limiting.

castro__limit_fluxes_on_small_dens="0"

# Variables we need to set up the collision.

problem="0"
collision_separation="2.0"
collision_impact_parameter="0.0"

# Disable refinement based on non-burning criteria.
# These can be re-enabled later for specific tests.

max_stellar_tagging_level="0"
max_temperature_tagging_level="0"
max_center_tagging_level="0"

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

# Maximum number of subcycles.

castro__max_subcycles="128"

# Enable efficient regridding (don't actually regrid if the grids haven't changed.)

amr__use_efficient_regrid="1"

# Enable reactions.

castro__do_react="1"

# Disable rotation.

castro__do_rotation="0"

# Ease up on the gravity tolerance since we're in axisymmetric and at high resolution.

gravity__abs_tol="5.e-8"

# Many of the collision papers in the literature use an equal
# C/O ratio  by mass in the initial white dwarfs. We will do
# this too for comparison purposes.

co_wd_c_frac="0.5d0"
co_wd_o_frac="0.5d0"

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

mass_P="0.64"
mass_S="0.64"

# Run four helium shell mass cases. These are the same
# values from Holcomb and Kushnir (2016), so that we
# can make a direct comparison (see Figure 2).

helium_shell_mass_list="0.00 0.04 0.08 0.16"
helium_shell_mass_list="0.16"

stop_time="4.0"
prob_lo="0.00e9"
prob_hi="5.12e9"

ncell_list="256 512 1024 2048 4096"

for helium_shell_mass in $helium_shell_mass_list
do

    co_wd_he_shell_mass=$helium_shell_mass

    for ncell in $ncell_list
    do

        hydro_refinement_list="1 2 4 8 16 32 64 128 256 512"

        for hydro_refinement in $hydro_refinement_list
        do

            dir=$results_dir/he$helium_shell_mass/n$ncell/r$hydro_refinement
            set_run_opts

            if [ $to_run -eq 1 ]; then
                run
            fi

        done

    done

done
