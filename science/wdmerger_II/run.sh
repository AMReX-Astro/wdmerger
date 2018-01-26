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

        if   [ $ncell -eq 16384 ]; then
            nprocs="2048"
            walltime="6:00:00"
        elif [ $ncell -eq 8192 ]; then
            nprocs="2048"
            walltime="6:00:00"
        elif [ $ncell -eq 4096 ]; then
            nprocs="1024"
            walltime="2:00:00"

            if [ $stellar_refinement -eq 8 ]; then
                nprocs="2048"
                walltime="6:00:00"
            elif [ $stellar_refinement -eq 16 ]; then
                nprocs="2048"
                walltime="6:00:00"
            fi
        elif [ $ncell -eq 2048 ]; then
            nprocs="512"
            walltime="2:00:00"

            if [ ! -z $stellar_refinement ]; then
                if [ $stellar_refinement -eq 16 ]; then
                    nprocs="1024"
                    walltime="2:00:00"
                fi
            fi
        elif [ $ncell -eq 1024 ]; then
            nprocs="256"
            walltime="2:00:00"

            if [ ! -z $stellar_refinement ]; then
                if [ $stellar_refinement -eq 16 ]; then
                    nprocs="512"
                    walltime="2:00:00"
                fi
            fi
        elif [ $ncell -eq 512 ]; then
            nprocs="128"
            walltime="2:00:00"

            if [ ! -z $stellar_refinement ]; then
                if [ $stellar_refinement -eq 8 ]; then
                    nprocs="512"
                    walltime="2:00:00"
                fi
            fi
        elif [ $ncell -eq 256 ]; then
            nprocs="32"
            walltime="2:00:00"

            if [ ! -z $stellar_refinement ]; then
                if [ $stellar_refinement -eq 16 ]; then
                    nprocs="128"
                    walltime="2:00:00"
                elif [ $stellar_refinement -eq 32 ]; then
                    nprocs="256"
                    walltime="2:00:00"
                fi
            fi
        else
            echoerr "Unknown number of cells per dimension."
        fi

    else

        echoerr "This machine is not set up for this job."

    fi

    # Set up the 2D geometry appropriately.

    amr_n_cell="$(echo "$ncell / 2" | bc) $ncell"
    geometry_prob_lo="0.0 $prob_lo"
    geometry_prob_hi="$prob_hi $prob_hi"

    if [ $ncell -eq 256 ]; then

        amr_blocking_factor="32"
	amr_max_grid_size="64 64 64 64 128 128 256 256 512 512"

    elif [ $ncell -eq 512 ]; then

	amr_blocking_factor="32"
        amr_max_grid_size="64 64 64 64 128 128 256 256 512 512"

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

# Needed for the makefile: we want to compile in 2D for these tests.

DIM="2"

# Use the aprox13 network for all following tests.

NETWORK_DIR="aprox13"

# Get the right inputs and probin files.

inputs="inputs_2d"
probin="probin"

# Use a narrower domain than the usual default.
# This problem is insensitive to accuracy in the boundary conditions.

prob_lo="-4e9"
prob_hi="4e9"

# Variables we need to set up the collision.

problem="0"
collision_separation="4.0"
collision_impact_parameter="0.0"

# Disable refinement based on non-burning criteria.
# These can be re-enabled later for specific tests.

max_stellar_tagging_level="0"
max_temperature_tagging_level="0"
max_center_tagging_level="0"

# Make a full plotfile every second.

amr_plot_per="1.0"
amr_plot_vars="ALL"
amr_derive_plot_vars="ALL"

# Make small plotfiles rapidly.

amr_small_plot_per="0.01"
amr_small_plot_vars="density Temp rho_e rho_c12 rho_o16 rho_si28 rho_ni56 enuc"
amr_derive_small_plot_vars="pressure soundspeed x_velocity y_velocity t_sound_t_enuc"

# Ensure that plotting intervals are hit exactly.
# This helps in resolving the detonation point.

castro_plot_per_is_exact="1"
castro_small_plot_per_is_exact="1"

# Save checkpoints every second.

amr_check_per="1.0"

# Initial timestep shortening factor.

castro_init_shrink="0.1"

# Set the interval for doing diagnostic global sums.

castro_sum_interval="1"

# The interesting part of the evolution in this problem finishes by t = 9s.

stop_time="9.0"

# Enable reactions.

castro_do_react="1"

# Disable rotation.

castro_do_rotation="0"

# The collision papers in the literature all use an equal C/O ratio 
# by mass in the initial white dwarfs. We will do this too for 
# comparison purposes.

co_wd_c_frac="0.5d0"
co_wd_o_frac="0.5d0"

# Use the sponge to damp out noise in the ambient medium,
# which can be rather violent due to how much churn
# the colliding white dwarfs create.

castro_do_sponge="1"
sponge_lower_density="1.0d-1"
sponge_upper_density="1.0d0"

# The timesteps can get quite small if you're fully 
# resolving the burning, so allow for this.

castro_dt_cutoff="1.0e-12"

# Allow the timestep to change by up to 5% per advance.
# We don't need to be super strict on the timestep control
# for this problem, we have empirically found.

castro_change_max="1.05"

# Some defaults.

ncell_default="256"
dtnuc_e_default="1.e200"
dtnuc_X_default="1.e200"
dxnuc_default="1.0e200"
dxnuc_max_default="1.0e200"
mass_P_default="0.64"
mass_S_default="0.64"

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



ncell_list="256 512 1024 2048 4096 8192"

for ncell in $ncell_list
do

    # The tolerance on the gravity solve needs to
    # be loosened with resolution, we have found
    # experimentally. This is mostly a problem with
    # the cylindrical grid we are using.

    if   [ $ncell -ge 8192 ]; then
        gravity_abs_tol="1.e-8"
    elif [ $ncell -eq 4096 ]; then
        gravity_abs_tol="1.e-9"
    elif [ $ncell -eq 2048 ]; then
        gravity_abs_tol="5.e-9"
    else
        gravity_abs_tol="1.e-10"
    fi

    stellar_refinement_list=""

    if   [ $ncell -eq 256 ]; then
        stellar_refinement_list="1"
    elif [ $ncell -eq 512 ]; then
        stellar_refinement_list="1"
    elif [ $ncell -eq 1024 ]; then
        stellar_refinement_list="1 4 8 16"
    elif [ $ncell -eq 2048 ]; then
        stellar_refinement_list="1 16"
    elif [ $ncell -eq 4096 ]; then
        stellar_refinement_list="1"
    else
        stellar_refinement_list="1"
    fi

    for stellar_refinement in $stellar_refinement_list
    do

        base_dir=$results_dir/collision_2D/mass_P_$mass_P/mass_S_$mass_S/n$ncell/r$stellar_refinement
        start_dir=$base_dir/start

        start_done="0"

        if [ -d $start_dir ]; then
            dir=$start_dir
            start_done=$(is_dir_done)
        fi

        if [ $start_done -ne 1 ]; then

            refinement=1

            dir=$start_dir
            set_run_opts

            # First, we need to do the initial run, up to the point
            # when burning starts.

            castro_ts_te_stopping_criterion="1.0e-6"

            if [ $to_run -eq 1 ]; then
                run
            fi

            unset castro_ts_te_stopping_criterion
            unset refinement

        else

            # At this point we should have our starting checkpoint in $start_dir.
            # Now we can do the runs by copying the checkpoint.

            burning_mode_list="self-heat suppressed"

            for burning_mode_str in $burning_mode_list
            do

                if [ $burning_mode_str == "self-heat" ]; then

                    burning_mode="1"

                    rtol_spec="1.d-6"
                    atol_spec="1.d-6"
                    rtol_temp="1.d-6"
                    atol_temp="1.d-6"
                    rtol_enuc="1.d-6"
                    atol_enuc="1.d-6"

                elif [ $burning_mode_str == "suppressed" ]; then

                    burning_mode="3"

                    # Use tighter tolerances for the suppressed burn;
                    # this seems to prevent integration failures.

                    rtol_spec="1.d-8"
                    atol_spec="1.d-8"
                    rtol_temp="1.d-8"
                    atol_temp="1.d-8"
                    rtol_enuc="1.d-8"
                    atol_enuc="1.d-8"

                fi

                # Only do the suppressed burn at certain resolutions.

                if [ $burning_mode_str == "suppressed" ]; then
                    if [ $ncell -gt 4096 ]; then
                        continue
                    fi
                    if [ $stellar_refinement -gt 1 ]; then
                        continue
                    fi
                fi

                # Complete the run with no special options.

                refinement=1

                dir=$base_dir/$burning_mode_str/finish
                set_run_opts

                copy_checkpoint

                if [ $to_run -eq 1 ]; then
                    run
                fi



                # The above is the only run we want to do with the suppressed burn.

                if [ $burning_mode_str == "suppressed" ]; then
                    continue
                fi



                # Test the effect of the burning timestep limiter parameter values.

                dtnuc_e_list=""
                dtnuc_X_list=""
                dtnuc_eX_list=""

                if [ $ncell -eq 256 ] && [ $stellar_refinement -eq 1 ]; then

                    dtnuc_e_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 5.0e0 2.0e0 1.0e0 5.0e-1 2.0e-1 1.0e-1 5.0e-2 2.0e-2 1.0e-2"
                    dtnuc_X_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 5.0e0 2.0e0 1.0e0 5.0e-1 2.0e-1 1.0e-1"
                    dtnuc_eX_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 5.0e0 2.0e0 1.0e0 5.0e-1 2.0e-1 1.0e-1"

                elif [ $ncell -eq 512 ] && [ $stellar_refinement -eq 1 ]; then

                    dtnuc_e_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0 5.0e-1 2.0e-1 1.0e-1"
                    dtnuc_X_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0 5.0e-1 2.0e-1 1.0e-1"
                    dtnuc_eX_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0 5.0e-1 2.0e-1 1.0e-1"

                elif [ $ncell -eq 1024 ] && [ $stellar_refinement -eq 1 ]; then

                    dtnuc_e_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0"
                    dtnuc_X_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0"
                    dtnuc_eX_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0"

                elif [ $ncell -eq 2048 ] && [ $stellar_refinement -eq 1 ]; then

                    dtnuc_e_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2"
                    dtnuc_X_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2"
                    dtnuc_eX_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2"

                fi

                for dtnuc in $dtnuc_e_list
                do

                    castro_dtnuc_e=$dtnuc

                    # For these tests only, disable timestep limiting for the plotting,
                    # since the whole point is to test the effect of the timestep.
                    # We'll still get enough plotfiles to examine the evolution if
                    # we choose to do so.

                    castro_plot_per_is_exact="0"
                    castro_small_plot_per_is_exact="0"

                    dir=$base_dir/$burning_mode_str/burning_limiter_e/dt$dtnuc
                    set_run_opts

                    copy_checkpoint

                    if [ $to_run -eq 1 ]; then
                        run
                    fi

                    castro_plot_per_is_exact="1"
                    castro_small_plot_per_is_exact="1"

                done

                castro_dtnuc_e=$dtnuc_e_default



                for dtnuc in $dtnuc_X_list
                do

                    castro_dtnuc_X=$dtnuc

                    castro_plot_per_is_exact="0"
                    castro_small_plot_per_is_exact="0"

                    dir=$base_dir/$burning_mode_str/burning_limiter_X/dt$dtnuc
                    set_run_opts

                    copy_checkpoint

                    if [ $to_run -eq 1 ]; then
                        run
                    fi

                    castro_plot_per_is_exact="1"
                    castro_small_plot_per_is_exact="1"

                done

                castro_dtnuc_X=$dtnuc_X_default



                for dtnuc in $dtnuc_eX_list
                do

                    castro_dtnuc_e=$dtnuc
                    castro_dtnuc_X=$dtnuc

                    castro_plot_per_is_exact="0"
                    castro_small_plot_per_is_exact="0"

                    dir=$base_dir/$burning_mode_str/burning_limiter_eX/dt$dtnuc
                    set_run_opts

                    copy_checkpoint

                    if [ $to_run -eq 1 ]; then
                        run
                    fi

                    castro_plot_per_is_exact="1"
                    castro_small_plot_per_is_exact="1"

                done

                castro_dtnuc_e=$dtnuc_e_default
                castro_dtnuc_X=$dtnuc_X_default



                # Do runs with stellar refinement up to a selected level.

                refinement_list=""

                if [ $ncell -eq 256 ]; then
                    if [ $stellar_refinement -eq 1 ]; then
                        refinement_list="2 4 8 16"
                    fi
                fi

                for refinement in $refinement_list
                do

                    full_stellar_refinement=1

                    dir=$base_dir/$burning_mode_str/stellar/r$refinement
                    set_run_opts

                    copy_checkpoint

                    if [ $to_run -eq 1 ]; then
                        run
                    fi

                    unset full_stellar_refinement

                done



                # Do the runs with temperature-based refinement.

                temp_list="1.0d9"

                for temperature_tagging_threshold in $temp_list
                do

                    refinement_list=""

                    if [ $ncell -eq 1024 ]; then
                        if [ $stellar_refinement -eq 16 ]; then
                            refinement_list="4"
                        fi
                    fi

                    for refinement in $refinement_list
                    do

                        dir=$base_dir/$burning_mode_str/temperature/r$refinement
                        set_run_opts

                        max_temperature_tagging_level=$amr_max_level

                        copy_checkpoint

                        if [ $to_run -eq 1 ]; then
                            run
                        fi

                        unset max_temperature_tagging_level

                    done

                done

                max_temperature_tagging_level="0"



                # Do runs with burning-based AMR up to a selected level.

                dxnuc_list="1.0e-2"

                for castro_dxnuc in $dxnuc_list
                do

                    refinement_list=""

                    if [ $ncell -eq 256 ]; then
                        if [ $stellar_refinement -eq 1 ]; then
                            refinement_list="2 4 8 16 32"
                        fi
                    elif [ $ncell -eq 512 ]; then
                        if [ $stellar_refinement -eq 1 ]; then
                            refinement_list="2 4 8 16 32"
                        fi
                    elif [ $ncell -eq 1024 ]; then
                        if [ $stellar_refinement -eq 1 ]; then
                            refinement_list="2 4 8 16 32 64 128"
                        elif [ $stellar_refinement -eq 4 ]; then
                            refinement_list="2 4 8 16 32 64"
                        elif [ $stellar_refinement -eq 16 ]; then
                            refinement_list="2 4 8 16 32"
                        fi
                    elif [ $ncell -eq 2048 ]; then
                        if [ $stellar_refinement -eq 1 ]; then
                            refinement_list="2 4 8 16 32 64"
                        fi
                    fi

                    for refinement in $refinement_list
                    do

                        dir=$base_dir/$burning_mode_str/dxnuc/f$castro_dxnuc/r$refinement/
                        set_run_opts

                        copy_checkpoint

                        if [ $to_run -eq 1 ]; then
                            run
                        fi

                    done

                done

                castro_dxnuc=$dxnuc_default



                # Do runs with static refinement in the center up to a selected level.

                center_tagging_radius="2.0d8"

                refinement_list=""

                if [ $ncell -eq 1024 ]; then
                    if [ $stellar_refinement -eq 16 ]; then
                        refinement_list="2"
                    fi
                fi

                for refinement in $refinement_list
                do

                    max_center_tagging_level=$amr_max_level

                    dir=$base_dir/$burning_mode_str/center/d$center_tagging_radius/r$refinement/
                    set_run_opts

                    copy_checkpoint

                    if [ $to_run -eq 1 ]; then
                        run
                    fi

                    unset max_center_tagging_level

                done

                unset center_tagging_radius

            done

            burning_mode=$burning_mode_default

            rtol_spec=$spec_tol_default
            atol_spec=$spec_tol_default
            rtol_temp=$temp_tol_default
            atol_temp=$temp_tol_default
            rtol_enuc=$enuc_tol_default
            atol_enuc=$enuc_tol_default

        fi

    done

done
