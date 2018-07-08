#!/bin/bash

source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

if [ -z $ncell ]; then
    echoerr "ncell not set; exiting."
    exit
fi

# The following assumes we are using Titan.

if [ $MACHINE == "TITAN" ]; then

    queue="batch"

    threads_per_task=4

    nprocs="16"
    walltime="2:00:00"

else

    echoerr "This machine is not set up for this job."

fi

# Set up the gridding.

amr_n_cell="$ncell"
amr_ref_ratio="2"

if [ ! -z $refinement ]; then

    if [ $refinement -eq 1 ]; then
        amr_max_level=0
    elif [ $refinement -eq 2 ]; then
        amr_max_level=1
    elif [ $refinement -eq 4 ]; then
        amr_max_level=2
    elif [ $refinement -eq 8 ]; then
        amr_max_level=3
    elif [ $refinement -eq 16 ]; then
        amr_max_level=4
    elif [ $refinement -eq 32 ]; then
        amr_max_level=5
    elif [ $refinement -eq 64 ]; then
        amr_max_level=6
    elif [ $refinement -eq 128 ]; then
        amr_max_level=7
    elif [ $refinement -eq 256 ]; then
        amr_max_level=8
    elif [ $refinement -eq 512 ]; then
        amr_max_level=9
    elif [ $refinement -eq 1024 ]; then
        amr_max_level=10
    elif [ $refinement -eq 2048 ]; then
        amr_max_level=11
    elif [ $refinement -eq 4096 ]; then
        amr_max_level=12
    elif [ $refinement -eq 8192 ]; then
        amr_max_level=13
    elif [ $refinement -eq 16384 ]; then
        amr_max_level=14
    elif [ $refinement -eq 32768 ]; then
        amr_max_level=15
    elif [ $refinement -eq 65536 ]; then
        amr_max_level=16
    elif [ $refinement -eq 131072 ]; then
        amr_max_level=17
    elif [ $refinement -eq 262144 ]; then
        amr_max_level=18
    elif [ $refinement -eq 524288 ]; then
        amr_max_level=19
    elif [ $refinement -eq 1048576 ]; then
        amr_max_level=20
    elif [ $refinement -eq 2097152 ]; then
        amr_max_level=21
    elif [ $refinement -eq 4194304 ]; then
        amr_max_level=22
    elif [ $refinement -eq 8388608 ]; then
        amr_max_level=23
    else
        echo "Unknown refinement factor: "$refinement"; exiting."
        exit
    fi

    max_temperr_lev=$amr_max_level
    max_tempgrad_rel_lev=$amr_max_level

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

problem_dir=$CASTRO_HOME/Exec/science/Detonation

use_first_castro_ex="1"

# Use the aprox13 network for all following tests.

NETWORK_DIR="aprox13"

# Get the right inputs and probin files.

inputs="inputs-collision"
probin="probin-collision"

# Set plotfile details.

amr_plot_per="-1.0"
amr_plot_vars="ALL"
amr_derive_plot_vars="ALL"

amr_small_plot_per="1.0e-2"
amr_small_plot_vars="density Temp rho_e rho_c12 rho_o16 rho_si28 rho_ni56 enuc"
amr_derive_small_plot_vars="pressure soundspeed x_velocity y_velocity t_sound_t_enuc"

# Ensure that plotting intervals are hit exactly.
# This helps in resolving the detonation point.

castro_plot_per_is_exact="1"
castro_small_plot_per_is_exact="1"

# Checkpoint save rate.

amr_check_per="-1.0"

# Checkpoints should not require plotfiles.

amr_write_plotfile_with_checkpoint="0"

# Initial timestep shortening factor.

castro_init_shrink="0.1"

# Set the interval for doing diagnostic global sums.

castro_sum_interval="1"

# Enable reactions.

castro_do_react="1"

# The timesteps can get quite small if you're fully 
# resolving the burning, so allow for this.

castro_dt_cutoff="1.0e-20"

# Allow the timestep to change by up to 5% per advance.

castro_change_max="1.05"

# Some defaults.

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

castro_output_at_completion="1"

castro_use_stopping_criterion="1"
castro_T_stopping_criterion="5.0e9"
castro_ts_te_stopping_criterion="1.0e200"

amr_subcycling_mode="None"

stop_time="10.0"
max_step="10000000"

size_list="8.192e8"

dens_list="5.d6"

g_list="1.1e8"

v_list="1.0d8"

T_l="1.0d7"
T_r="1.0d7"

ncell_list="32 64 128 256 512 1024 2048 4096 8192"

dtnuc_list="1.0e200"

burning_mode_list="self-heat suppressed"

dxnuc_list="1.0e200 1.0e-1"
temperr_list="1.0d20"
tempgrad_rel_list="1.0d20 0.5"

amr_n_error_buf="2"

# This loop nest covers all of the runs we will do,
# and is split into two basic parts. The outer loops
# govern the high level problem setup choices --
# the size of the domain, the initial velocity
# and gravitational acceleration, and the type
# of nuclear burning we are doing. We'll also cover
# the temporal resolution in this section.

# The inner part of the loop nest covers the
# spatial resolution options.

# We split it up this way because for all of
# the options that don't involve AMR, we can
# group the runs together in a single base
# case with uniform resolution, rather than
# redundantly calculate this uniform grid
# case for each one of those options.

for size in $size_list
do

    geometry_prob_lo="0.0"
    geometry_prob_hi="$size"

    for dens in $dens_list
    do

        for g in $g_list
        do

            gravity_const_grav=-$g

            for vel in $v_list
            do

                for burning_mode_str in $burning_mode_list
                do

                    if [ $burning_mode_str == "self-heat" ]; then
                        burning_mode="1"
                    elif [ $burning_mode_str == "suppressed" ]; then
                        burning_mode="3"
                    fi

                    if [ $burning_mode_str == "suppressed" ]; then
                        continue
                    fi

                    # Now we'll do the time resolution options.

                    for dtnuc in $dtnuc_list
                    do

                        castro_dtnuc_e=$dtnuc
                        castro_dtnuc_X=$dtnuc

                        # At this point we've completed all the non-spatial resolution
                        # options, and we move on to the number of zones in the coarse
                        # grid, and what AMR we will be doing.

                        base_dir=$results_dir/size$size/dens$dens/g$g/vel$vel/$burning_mode_str/dtnuc$dtnuc

                        for ncell in $ncell_list
                        do

                            # Before we reach the AMR section, complete a run
                            # using only the coarse grid.

                            refinement=1

                            dir=$base_dir/n$ncell/base
                            set_run_opts

                            if [ $to_run -eq 1 ]; then
                                run
                            fi

                            # Keep track of whether any of our inner loops
                            # actually enable AMR. If not, we don't need to
                            # do the run.

                            to_run=0

                            for castro_dxnuc in $dxnuc_list
                            do

                                if [ "$castro_dxnuc" != "1.0e200" ]; then
                                    to_run=1
                                fi

                                for temperr in $temperr_list
                                do

                                    if [ "$temperr" != "1.0d20" ]; then
                                        to_run=1
                                    fi

                                    for tempgrad_rel in $tempgrad_rel_list
                                    do

                                        if [ "$tempgrad_rel" != "1.0d20" ]; then
                                            to_run=1
                                        fi

                                        refinement_list=""

                                        if   [ $ncell -eq 32 ]; then
                                            refinement_list="2 4 8 16 32 64 128 256"
                                        elif [ $ncell -eq 64 ]; then
                                            refinement_list="2 4 8 16 32 64 128"
                                        elif [ $ncell -eq 128 ]; then
                                            refinement_list="2 4 8 16 32 64"
                                        elif [ $ncell -eq 256 ]; then
                                            refinement_list="2 4 8 16 32"
                                        elif [ $ncell -eq 512 ]; then
                                            refinement_list="2 4 8 16"
                                        elif [ $ncell -eq 1024 ]; then
                                            refinement_list="2 4 8"
                                        elif [ $ncell -eq 2048 ]; then
                                            refinement_list="2 4"
                                        elif [ $ncell -eq 4096 ]; then
                                            refinement_list="2 4"
                                        elif [ $ncell -eq 8192 ]; then
                                            refinement_list="2 4"
                                        fi

                                        for refinement in $refinement_list
                                        do

                                            # Now we've reached the point where we'll actually launch the AMR runs.

                                            # We should not be doing any r == 1 cases here; dummy check against that.

                                            if [ $refinement -eq 1 ]; then
                                                continue
                                            fi

                                            # Only do high temporal resolution runs for certain parameter combinations.

                                            if [ "$dtnuc" != "1.0e200" ]; then

                                                if [ $refinement -gt 1 ]; then
                                                    continue
                                                fi

                                                if [ $ncell -gt 1024 ]; then
                                                    continue
                                                fi

                                            fi

                                            dir=$base_dir/n$ncell/dxnuc$castro_dxnuc/temperr$temperr/tempgrad$tempgrad_rel/r$refinement
                                            set_run_opts

                                            if [ $to_run -eq 1 ]; then
                                                run
                                            fi

                                        done # refinement

                                    done # tempgrad

                                done # temperr

                            done # dxnuc

                            # Reset the run flag for the next series of iterations.

                            to_run=1

                        done # ncell

                    done # dtnuc

                done # burning mode

            done # vel

        done # g

    done # dens

done # size
