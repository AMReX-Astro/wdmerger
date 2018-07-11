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

    if [ ! -z $dxnuc_r ]; then
        castro_max_dxnuc_lev=$(log_base_2 $dxnuc_r)
    else
        castro_max_dxnuc_lev=0
    fi

    if [ ! -z $tempgrad_r ]; then
        max_tempgrad_rel_lev=$(log_base_2 $tempgrad_r)
    else
        max_tempgrad_lev=0
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

amr_small_plot_per="-0.01"
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

geometry_prob_lo="0.0"
geometry_prob_hi="8.192e8"

ofrac_list="0.0d0"

dens_list="5.d6"

g_list="1.1e8"

v_list="1.0d8"

T_l="1.0d7"
T_r="1.0d7"

ncell_list="16 24 32 48 64 96 128 256 512 1024 2048 4096"

burning_mode_list="self-heat suppressed"

dxnuc_list="1.0e-1"
tempgrad_rel_list="0.5"

amr_n_error_buf="2"

# This loop nest covers all of the runs we will do,
# and is split into two basic parts. The outer loops
# govern the high level problem setup choices --
# the oxygen fraction, the initial velocity
# and gravitational acceleration, and the type
# of nuclear burning we are doing.

# The inner part of the loop nest covers the
# spatial resolution options.

# We split it up this way because for all of
# the options that don't involve AMR, we can
# group the runs together in a single base
# case with uniform resolution, rather than
# redundantly calculate this uniform grid
# case for each one of those options.

for ofrac in $ofrac_list
do

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

                    # At this point we've completed all the non-spatial resolution
                    # options, and we move on to the number of zones in the coarse
                    # grid, and what AMR we will be doing.
                        
                    base_dir=$results_dir/ofrac$ofrac/dens$dens/g$g/vel$vel/$burning_mode_str

                    for ncell in $ncell_list
                    do

                        # Low resolution runs tend to skip right through the burning
                        # and have a very rough time with it as a result. Severely
                        # limit the CFL number for these runs to compensate.

                        if [ $ncell -lt 128 ]; then
                            castro_cfl="0.01"
                        else
                            castro_cfl="0.5"
                        fi

                        # Before we reach the AMR section, complete a run
                        # using only the coarse grid.

                        refinement=1

                        dir=$base_dir/n$ncell/base
                        set_run_opts

                        if [ $to_run -eq 1 ]; then
                            run
                        fi

                        tempgrad_r_list="1"
                        dxnuc_r_list="1"

                        if   [ $ncell -eq 256 ]; then
                            tempgrad_r_list="1 4 16"
                            dxnuc_r_list="1 4 16 64 256"
                        fi

                        for tempgrad_rel in $tempgrad_rel_list
                        do

                            for dxnuc in $dxnuc_list
                            do

                                castro_dxnuc=$dxnuc

                                for tempgrad_r in $tempgrad_r_list
                                do

                                    for dxnuc_r in $dxnuc_r_list
                                    do

                                        # Now we've reached the point where we'll actually launch the AMR runs.

                                        refinement=$dxnuc_r

                                        # We should not be doing any r == 1 cases here; dummy check against that.

                                        if [ $refinement -eq 1 ]; then
                                            continue
                                        fi

                                        # Skip runs where dxnuc refinement doesn't add on top of tempgrad.

                                        if [ $dxnuc_r -lt $tempgrad_r ]; then
                                            continue
                                        fi

                                        dir=$base_dir/n$ncell/tempgrad$tempgrad_rel/dxnuc$dxnuc/tempgrad_r$tempgrad_r/dxnuc_r$dxnuc_r
                                        set_run_opts

                                        if [ $to_run -eq 1 ]; then
                                            run
                                        fi

                                        to_run=1

                                    done # dxnuc_r

                                done # tempgrad_r

                            done # dxnuc

                        done # tempgrad

                    done # ncell

                done # burning mode

            done # vel

        done # g

    done # dens

done # ofrac
