#!/bin/bash

source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

    if [ -z $ncell ]; then
	echoerr "ncell not set; exiting."
	exit
    fi

    to_run=1

    # The following assumes we are using Titan.

    if [ $MACHINE == "TITAN" ]; then

        queue="batch"

        threads_per_task=4

        if   [ $ncell -eq 128 ]; then
            nprocs="64"
            walltime="2:00:00"
        elif [ $ncell -eq 256 ]; then
            nprocs="256"
            walltime="2:00:00"

            if [ ! -z $stellar_refinement ]; then

                if [ $stellar_refinement -eq 256 ]; then
                    nprocs="2048"
                    walltime="2:00:00"
                fi

            fi
        elif [ $ncell -eq 512 ]; then
            nprocs="1024"
            walltime="2:00:00"
        elif [ $ncell -eq 1024 ]; then
            nprocs="2048"
            walltime="6:00:00"
        else
            echoerr "Unknown number of cells per dimension."
        fi

    else

        echoerr "This machine is not set up for this job."

    fi

    # Set up the point mass.

    point_mass_for_bc=$(printf "%f" $point_mass) # Convert to floating point for bc
    castro_point_mass=$(echo "$point_mass_for_bc * 1.9884*10^33" | bc -l)
    castro_point_mass=$(printf "%e" $castro_point_mass) # Convert back to scientific notation

    # Set up the geometry.

    amr_n_cell="$ncell $ncell $ncell"
    geometry_prob_lo="$prob_lo $prob_lo $prob_lo"
    geometry_prob_hi="$prob_hi $prob_hi $prob_hi"

    if   [ $ncell -eq 128 ]; then

        amr_blocking_factor="16"
        amr_max_grid_size="32 32 64 64 128 128 256 256 512 512"

    elif [ $ncell -eq 256 ]; then

        amr_blocking_factor="16"
	amr_max_grid_size="32 32 64 64 128 128 256 256 512 512"

    elif [ $ncell -eq 512 ]; then

	amr_blocking_factor="32"
        amr_max_grid_size="32 32 64 64 128 128 256 256 512 512"

    elif [ $ncell -eq 1024 ]; then

        amr_blocking_factor="64"
        amr_max_grid_size="64 64 128 128 256 256 512 512 1024 1024"

    elif [ $ncell -eq 2048 ]; then

        amr_blocking_factor="64"
	amr_max_grid_size="64 64 128 128 256 256 512 512 1024 1024"

    elif [ $ncell -eq 4096 ]; then

        amr_blocking_factor="64"
	amr_max_grid_size="128 128 256 256 512 512 1024 1024 2048 2048"

    elif [ $ncell -eq 8192 ]; then

        amr_blocking_factor="128"
	amr_max_grid_size="128 128 256 256 512 512 1024 1024 2048 2048"

    elif [ $ncell -eq 16384 ]; then

        amr_blocking_factor="128"
	amr_max_grid_size="128 128 256 256 512 512 1024 1024 2048 2048"

    fi

    if [ ! -z $stellar_refinement ]; then

        if [ $stellar_refinement -eq 1 ]; then

            max_stellar_tagging_level="0"

            if [ ! -z $refinement ]; then

                if [ $refinement -eq 1 ]; then
                    amr_max_level=0
                elif [ $refinement -eq 2 ]; then
                    amr_max_level=1
                    amr_ref_ratio="2"
                elif [ $refinement -eq 4 ]; then
                    amr_max_level=1
                    amr_ref_ratio="4"
                elif [ $refinement -eq 8 ]; then
                    amr_max_level=2
                    amr_ref_ratio="4 2"
                elif [ $refinement -eq 16 ]; then
                    amr_max_level=2
                    amr_ref_ratio="4 4"
                elif [ $refinement -eq 32 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 2"
                elif [ $refinement -eq 64 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 4"
                elif [ $refinement -eq 128 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 2"
                elif [ $refinement -eq 256 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 4"
                elif [ $refinement -eq 512 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 2"
                elif [ $refinement -eq 1024 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 4"
                fi

            fi

        elif [ $stellar_refinement -eq 2 ]; then

            max_stellar_tagging_level="1"

            if [ ! -z $refinement ]; then

                if [ $refinement -eq 1 ]; then
                    amr_max_level=1
                    amr_ref_ratio="2"
                elif [ $refinement -eq 2 ]; then
                    amr_max_level=2
                    amr_ref_ratio="2 2"
                elif [ $refinement -eq 4 ]; then
                    amr_max_level=2
                    amr_ref_ratio="2 4"
                elif [ $refinement -eq 8 ]; then
                    amr_max_level=3
                    amr_ref_ratio="2 4 2"
                elif [ $refinement -eq 16 ]; then
                    amr_max_level=3
                    amr_ref_ratio="2 4 4"
                elif [ $refinement -eq 32 ]; then
                    amr_max_level=4
                    amr_ref_ratio="2 4 4 2"
                elif [ $refinement -eq 64 ]; then
                    amr_max_level=4
                    amr_ref_ratio="2 4 4 4"
                elif [ $refinement -eq 128 ]; then
                    amr_max_level=5
                    amr_ref_ratio="2 4 4 4 2"
                elif [ $refinement -eq 256 ]; then
                    amr_max_level=5
                    amr_ref_ratio="2 4 4 4 4"
                elif [ $refinement -eq 512 ]; then
                    amr_max_level=6
                    amr_ref_ratio="2 4 4 4 4 2"
                fi

            fi

        elif [ $stellar_refinement -eq 4 ]; then

            max_stellar_tagging_level="1"

            if [ ! -z $refinement ]; then

                if [ $refinement -eq 1 ]; then
                    amr_max_level=1
                    amr_ref_ratio="4"
                elif [ $refinement -eq 2 ]; then
                    amr_max_level=2
                    amr_ref_ratio="4 2"
                elif [ $refinement -eq 4 ]; then
                    amr_max_level=2
                    amr_ref_ratio="4 4"
                elif [ $refinement -eq 8 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 2"
                elif [ $refinement -eq 16 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 4"
                elif [ $refinement -eq 32 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 2"
                elif [ $refinement -eq 64 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 4"
                elif [ $refinement -eq 128 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 2"
                elif [ $refinement -eq 256 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 4"
                fi

            fi

        elif [ $stellar_refinement -eq 8 ]; then

            max_stellar_tagging_level="2"

            if [ ! -z $refinement ]; then

                if [ $refinement -eq 1 ]; then
                    amr_max_level=2
                    amr_ref_ratio="4 2"
                elif [ $refinement -eq 2 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 2 2"
                elif [ $refinement -eq 4 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 2 4"
                elif [ $refinement -eq 8 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 2 4 2"
                elif [ $refinement -eq 16 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 2 4 4"
                elif [ $refinement -eq 32 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 2 4 4 2"
                elif [ $refinement -eq 64 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 2 4 4 4"
                elif [ $refinement -eq 128 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 2 4 4 4 2"
                fi

            fi

        elif [ $stellar_refinement -eq 16 ]; then

            max_stellar_tagging_level="2"

            if [ ! -z $refinement ]; then

                if [ $refinement -eq 1 ]; then
                    amr_max_level=2
                    amr_ref_ratio="4 4"
                elif [ $refinement -eq 2 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 2"
                elif [ $refinement -eq 4 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 4"
                elif [ $refinement -eq 8 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 2"
                elif [ $refinement -eq 16 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 4"
                elif [ $refinement -eq 32 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 2"
                elif [ $refinement -eq 64 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 4"
                fi

            fi

        elif [ $stellar_refinement -eq 32 ]; then

            max_stellar_tagging_level="3"

            if [ ! -z $refinement ]; then

                if [ $refinement -eq 1 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 2"
                elif [ $refinement -eq 2 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 2 2"
                elif [ $refinement -eq 4 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 2 4"
                elif [ $refinement -eq 8 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 2 4 2"
                elif [ $refinement -eq 16 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 2 4 4"
                elif [ $refinement -eq 32 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 2 4 4 2"
                elif [ $refinement -eq 64 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 2 4 4 4"
                fi

            fi

        elif [ $stellar_refinement -eq 64 ]; then

            max_stellar_tagging_level="3"

            if [ ! -z $refinement ]; then

                if [ $refinement -eq 1 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 4"
                fi

            fi

        elif [ $stellar_refinement -eq 128 ]; then

            max_stellar_tagging_level="4"

            if [ ! -z $refinement ]; then

                if [ $refinement -eq 1 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 2"
                fi

            fi

        elif [ $stellar_refinement -eq 256 ]; then

            max_stellar_tagging_level="4"

            if [ ! -z $refinement ]; then

                if [ $refinement -eq 1 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 4"
                fi

            fi

        fi

        # Allow stellar refinement all the way up to the max level.

        if [ ! -z $full_stellar_refinement ]; then
            if [ $full_stellar_refinement -eq 1 ]; then
                max_stellar_tagging_level=$amr_max_level
            fi
        fi

    fi

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

                start_checkpoint=$(get_last_checkpoint $start_dir)

                echo "Copying initial checkpoint to" $dir"."

                if [ -d "$start_dir/$start_checkpoint" ]; then
                    cp -r "$start_dir/$start_checkpoint" $dir
                    rm -f "$dir/$start_checkpoint/jobIsDone"

                    # Also, copy the diagnostic output files
                    # so there is an apparently continuous record.

                    cp $start_dir/*diag.out $dir
                else
                    echoerr "No initial checkpoint available, exiting."
                    exit
                fi

            fi

        fi

    fi

}

# Specify the problem directory.

problem_dir=$CASTRO_HOME/Exec/science/wdmerger

DIM="3"

NETWORK_DIR="aprox13"

# Get the right inputs and probin files.

inputs="inputs_3d"
probin="probin"

# Use a large domain to capture the initial distance between the WD and BH.

prob_lo="-1.28e11"
prob_hi=" 1.28e11"

# Variables we need to set up the TDE.

problem="5"
tde_separation="8.0"
tde_beta="6.0"

# Set up refinement initially only for the stellar material.

max_stellar_tagging_level="20"
max_temperature_tagging_level="0"
max_center_tagging_level="0"

# Make a full plotfile every second.

amr_plot_per="1.0"
amr_plot_vars="ALL"
amr_derive_plot_vars="ALL"

# Make small plotfiles more rapidly.

amr_small_plot_per="0.1"
amr_small_plot_vars="density Temp rho_e rho_c12 rho_o16 rho_si28 rho_ni56 enuc"
amr_derive_small_plot_vars="pressure soundspeed x_velocity y_velocity t_sound_t_enuc"

# Ensure that plotting intervals are hit exactly.
# This helps in resolving a detonation.

castro_plot_per_is_exact="1"
castro_small_plot_per_is_exact="1"

# Save checkpoints every second.

amr_check_per="1.0"

# Initial timestep shortening factor.

castro_init_shrink="0.1"

# Set the interval for doing diagnostic global sums.

castro_sum_interval="1"

# The interesting part of the evolution in this problem finishes by t = 100s.

stop_time="100.0"

# Enable reactions.

castro_do_react="1"

# Disable rotation.

castro_do_rotation="0"

# Equal C/O ratio by mass.

co_wd_c_frac="0.5d0"
co_wd_o_frac="0.5d0"

# Use the sponge to damp out noise in the ambient medium.

castro_do_sponge="1"
sponge_lower_density="1.0d-1"
sponge_upper_density="1.0d0"

# The timesteps can get quite small if you're fully
# resolving burning, so allow for this.

castro_dt_cutoff="1.0e-12"

# Allow the timestep to change by up to 5% per advance.

castro_change_max="1.05"

# Some defaults.

ncell_default="256"
dtnuc_e_default="1.e200"
dtnuc_X_default="1.e200"
dxnuc_default="1.0e200"
dxnuc_max_default="1.0e200"
mass_P_default="0.6"
mass_S_default="-1.0"

ncell=$ncell_default
mass_P=$mass_P_default
mass_S=$mass_S_default
castro_dtnuc_e=$dtnuc_e_default
castro_dtnuc_X=$dtnuc_X_default
castro_dxnuc=$dxnuc_default
castro_dxnuc_max=$dxnuc_max_default
castro_react_T_min="1.0e8"
castro_react_rho_min="1.0e6"
castro_small_temp="1.0e7"
point_mass="1.0e3"
castro_use_point_mass="1"

spec_tol_default="1.d-6"
temp_tol_default="1.d-6"
enuc_tol_default="1.d-6"

rtol_spec=$spec_tol_default
atol_spec=$spec_tol_default
rtol_temp=$temp_tol_default
atol_temp=$temp_tol_default
rtol_enuc=$enuc_tol_default
atol_enuc=$enuc_tol_default



# The flag we will use to determine whether to run the job.

to_run=1



ncell_list="256"

for ncell in $ncell_list
do

    stellar_refinement_list=""

    if   [ $ncell -eq 128 ]; then
        stellar_refinement_list="8 16 32 64 128"
    elif [ $ncell -eq 256 ]; then
        stellar_refinement_list="128 256"
    elif [ $ncell -eq 512 ]; then
        stellar_refinement_list="2 4"
    elif [ $ncell -eq 1024 ]; then
        stellar_refinement_list="64"
    fi

    for stellar_refinement in $stellar_refinement_list
    do

        # First do a run for a fixed amount of time with the point mass disabled,
        # to test how well the stars stay in HSE as a function of resolution.
        # Do this both for stars at rest, and at the same initial velocity as the orbit.

        stop_time="5.0"

        for tde_initial_velocity in 0 1
        do

            dir=$results_dir/WD_mass_$mass_P/BH_mass_$point_mass/HSE/motion$tde_initial_velocity/n$ncell/r$stellar_refinement

            refinement="1"

            castro_use_point_mass="0"
            castro_do_react="0"

            # Do not limit the timestep for this test.

            castro_plot_per_is_exact="0"
            castro_small_plot_per_is_exact="0"

            set_run_opts

            if [ $to_run -eq 1 ]; then
                run
            fi

            castro_use_point_mass="1"
            castro_do_react="1"

            castro_plot_per_is_exact="1"
            castro_small_plot_per_is_exact="1"

        done

        unset tde_initial_velocity

    done

done
