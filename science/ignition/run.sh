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

    nprocs="16"
    walltime="2:00:00"

else

    echoerr "This machine is not set up for this job."

fi

# Set up the gridding.

amr_n_cell="$ncell"

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
    elif [ $refinement -eq 2048 ]; then
        amr_max_level=6
        amr_ref_ratio="4 4 4 4 4 2"
    elif [ $refinement -eq 4096 ]; then
        amr_max_level=6
        amr_ref_ratio="4 4 4 4 4 4"
    elif [ $refinement -eq 8192 ]; then
        amr_max_level=7
        amr_ref_ratio="4 4 4 4 4 4 2"
    elif [ $refinement -eq 16384 ]; then
        amr_max_level=7
        amr_ref_ratio="4 4 4 4 4 4 4"
    elif [ $refinement -eq 32768 ]; then
        amr_max_level=8
        amr_ref_ratio="4 4 4 4 4 4 4 2"
    elif [ $refinement -eq 65536 ]; then
        amr_max_level=8
        amr_ref_ratio="4 4 4 4 4 4 4 4"
    elif [ $refinement -eq 131072 ]; then
        amr_max_level=9
        amr_ref_ratio="4 4 4 4 4 4 4 4 2"
    elif [ $refinement -eq 262144 ]; then
        amr_max_level=9
        amr_ref_ratio="4 4 4 4 4 4 4 4 4"
    elif [ $refinement -eq 524288 ]; then
        amr_max_level=10
        amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
    elif [ $refinement -eq 1048576 ]; then
        amr_max_level=10
        amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
    elif [ $refinement -eq 2097152 ]; then
        amr_max_level=11
        amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
    elif [ $refinement -eq 4194304 ]; then
        amr_max_level=11
        amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 4"
    elif [ $refinement -eq 8388608 ]; then
        amr_max_level=12
        amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 4 2"
    else
        echo "Unknown refinement factor: "$refinement"; exiting."
        exit
    fi

    max_tempgrad_rel_lev=$amr_max_level

    probin_tagging_tempgrad=$tempgrad
    probin_tagging_max_tempgrad_rel_lev=$max_tempgrad_rel_lev

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

amr_small_plot_per="-1.0"
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

ncell_list="32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536"
ncell_list="32 64 128 256 512 1024 2048 4096"

dens_list="5.d6"

g_list="1.1e8"

v_list="0.0d0"

T_l="1.0d7"
T_r="1.0d7"

dtnuc_list="1.0e200 1.0e10 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0 1.0e-1"
dtnuc_list="1.0e200"

dxnuc_list="1.0e-1"

burning_mode_list="self-heat suppressed"

size_list="8.192e8"

tempgrad_list="1.d20 2.0 0.5"

amr_n_error_buf="2"

for size in $size_list
do

   geometry_prob_lo="0.0"
   geometry_prob_hi="$size"

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

       for dtnuc in $dtnuc_list
       do

           castro_dtnuc_e=$dtnuc
           castro_dtnuc_X=$dtnuc

           for castro_dxnuc in $dxnuc_list
           do

               for tempgrad in $tempgrad_list
               do

                   for dens in $dens_list
                   do

                       for g in $g_list
                       do

                           gravity_const_grav=-$g

                           for vel in $v_list
                           do

                               for ncell in $ncell_list
                               do

                                   refinement_list="1"

                                   if   [ $ncell -eq 32 ]; then
                                       refinement_list="1 8388608"
                                       refinement_list="1 4 16 64 256 1024 4096"
                                   elif [ $ncell -eq 64 ]; then
                                       refinement_list="1 4194304"
                                       refinement_list="1 4 16 64 256 1024 4096"
                                   elif [ $ncell -eq 128 ]; then
                                       refinement_list="1 2097152"
                                       refinement_list="1 4 16 64 256 1024 4096"
                                   elif [ $ncell -eq 256 ]; then
                                       refinement_list="1 4 16 64 256 1024 4096 16384 65536 262144"
                                       refinement_list="1 4 16 64 256 1024 4096"
                                   elif [ $ncell -eq 512 ]; then
                                       refinement_list="1 524288"
                                       refinement_list="1 4 16 64 256 1024 4096 16384"
                                       refinement_list="1"
                                   elif [ $ncell -eq 1024 ]; then
                                       refinement_list="1 4 16 262144"
                                       refinement_list="1 4 16 64 256 1024"
                                       refinement_list="1"
                                   elif [ $ncell -eq 2048 ]; then
                                       refinement_list="1 131072"
                                       refinement_list="1 4 16 64 256 1024"
                                       refinement_list="1"
                                   elif [ $ncell -eq 4096 ]; then
                                       refinement_list="1 65536"
                                       refinement_list="1"
                                   elif [ $ncell -eq 8192 ]; then
                                       refinement_list="1 32768"
                                       refinement_list="1"
                                   elif [ $ncell -eq 16384 ]; then
                                       refinement_list="1 16384"
                                       refinement_list="1"
                                   elif [ $ncell -eq 32768 ]; then
                                       refinement_list="1"
                                   elif [ $ncell -eq 65536 ]; then
                                       refinement_list="1"
                                   fi

                                   for refinement in $refinement_list
                                   do

                                       dir=$results_dir/size$size/$burning_mode_str/dtnuc$dtnuc/dxnuc$castro_dxnuc/tempgrad$tempgrad/dens$dens/g$g/v$vel/n$ncell/r$refinement

                                       # Only do high temporal resolution runs for certain parameter combinations.

                                       if [ "$dtnuc" != "1.0e200" ]; then

                                           if [ $refinement -gt 1 ]; then
                                               continue
                                           fi

                                           if [ $ncell -gt 1024 ]; then
                                               continue
                                           fi

                                       fi

                                       # Only do temperature gradient refinement for a sufficiently large number of zones;
                                       # for low zone numbers, it basically results in refining the entire grid.
                                       # Also avoid it for the high resolution runs to avoid duplication.

                                       if [ "$tempgrad" != "1.d20" ]; then

                                           if [ $ncell -lt 256 ]; then
                                               continue
                                           fi

                                           if [ $ncell -gt 256 ]; then
                                               continue
                                           fi

                                       fi

                                       set_run_opts

                                       if [ $to_run -eq 1 ]; then
                                           run
                                       fi

                                       to_run="1"

                                   done

                               done

                           done

                       done

                   done

               done

           done

       done

   done

done
