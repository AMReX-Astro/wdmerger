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

        if   [ $ncell -eq 4096 ]; then
            nprocs="2048"
            walltime="6:00:00"
        elif [ $ncell -eq 2048 ]; then
            nprocs="2048"
            walltime="6:00:00"
        elif [ $ncell -eq 1024 ]; then
            nprocs="512"
            walltime="2:00:00"
        elif [ $ncell -eq 512 ]; then
            nprocs="128"
            walltime="2:00:00"
        elif [ $ncell -eq 256 ]; then
            nprocs="32"
            walltime="2:00:00"
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

    # Defaults for the coarse uniform grid.

    amr_max_level="0"
    amr_ref_ratio="2"
    amr_blocking_factor="8"
    amr_max_grid_size="32"

    # If we specify a refinement parameter,
    # we can use that to control amr_ref_ratio.

    if [ ! -z $refinement ]; then

        # Disable refinement based on non-burning criteria.

        max_stellar_tagging_level="0"
        max_temperature_tagging_level="0"
        max_center_tagging_level="0"

	if [ $ncell -eq 256 ]; then

            amr_blocking_factor="8"

	    if [ $refinement -eq 1 ]; then
		amr_max_level=0
		amr_max_grid_size="32"
	    elif [ $refinement -eq 2 ]; then
		amr_max_level=1
		amr_ref_ratio="2"
		amr_max_grid_size="32 32"
	    elif [ $refinement -eq 4 ]; then
		amr_max_level=1
		amr_ref_ratio="4"
		amr_max_grid_size="32 32"
	    elif [ $refinement -eq 8 ]; then
		amr_max_level=2
		amr_ref_ratio="4 2"
		amr_max_grid_size="32 32 64"
	    elif [ $refinement -eq 16 ]; then
		amr_max_level=2
		amr_ref_ratio="4 4"
		amr_max_grid_size="32 32 64"
	    elif [ $refinement -eq 32 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 2"
		amr_max_grid_size="32 32 64 64"
	    elif [ $refinement -eq 64 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 4"
		amr_max_grid_size="32 32 64 64"
	    elif [ $refinement -eq 128 ]; then
		amr_max_level=4
		amr_ref_ratio="4 4 4 2"
		amr_max_grid_size="32 32 64 64 128"
	    elif [ $refinement -eq 256 ]; then
		amr_max_level=4
		amr_ref_ratio="4 4 4 4"
		amr_max_grid_size="32 32 64 64 128"
	    elif [ $refinement -eq 512 ]; then
		amr_max_level=5
		amr_ref_ratio="4 4 4 4 2"
		amr_max_grid_size="32 32 64 64 128 128"
	    elif [ $refinement -eq 1024 ]; then
		amr_max_level=5
		amr_ref_ratio="4 4 4 4 4"
		amr_max_grid_size="32 32 64 64 128 128"
	    fi

        elif [ $ncell -eq 512 ]; then

	    amr_blocking_factor="16"

	    if [ $refinement -eq 1 ]; then
		amr_max_level=0
		amr_max_grid_size="64"
	    elif [ $refinement -eq 2 ]; then
		amr_max_level=1
		amr_ref_ratio="2"
		amr_max_grid_size="64 64"
	    elif [ $refinement -eq 4 ]; then
		amr_max_level=1
		amr_ref_ratio="4"
		amr_max_grid_size="64 64"
	    elif [ $refinement -eq 8 ]; then
		amr_max_level=2
		amr_ref_ratio="4 2"
		amr_max_grid_size="64 64 128"
	    elif [ $refinement -eq 16 ]; then
		amr_max_level=2
		amr_ref_ratio="4 4"
		amr_max_grid_size="64 64 128"
	    elif [ $refinement -eq 32 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 2"
		amr_max_grid_size="64 64 128 128"
	    elif [ $refinement -eq 64 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 4"
		amr_max_grid_size="64 64 128 128"
            elif [ $refinement -eq 128 ]; then
                amr_max_level=4
                amr_ref_ratio="4 4 4 2"
                amr_max_grid_size="64 64 128 128 256"
            elif [ $refinement -eq 256 ]; then
                amr_max_level=4
                amr_ref_ratio="4 4 4 4"
                amr_max_grid_size="64 64 128 128 256"
	    fi

	elif [ $ncell -eq 1024 ]; then

            amr_blocking_factor="32"

	    if [ $refinement -eq 1 ]; then
		amr_max_level=0
		amr_max_grid_size="128"
	    elif [ $refinement -eq 2 ]; then
		amr_max_level=1
		amr_ref_ratio="2"
		amr_max_grid_size="128 128"
	    elif [ $refinement -eq 4 ]; then
		amr_max_level=1
		amr_ref_ratio="4"
		amr_max_grid_size="128 128"
	    elif [ $refinement -eq 8 ]; then
		amr_max_level=2
		amr_ref_ratio="4 2"
		amr_max_grid_size="128 128 256"
	    elif [ $refinement -eq 16 ]; then
		amr_max_level=2
		amr_ref_ratio="4 4"
		amr_max_grid_size="128 128 256"
	    elif [ $refinement -eq 32 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 2"
		amr_max_grid_size="128 128 256 256"
	    elif [ $refinement -eq 64 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 4"
		amr_max_grid_size="128 128 256 256"
            elif [ $refinement -eq 128 ]; then
                amr_max_level=4
                amr_ref_ratio="4 4 4 2"
                amr_max_grid_size="128 128 256 256 512"
	    fi

        elif [ $ncell -eq 2048 ]; then

            amr_blocking_factor="64"

	    if [ $refinement -eq 1 ]; then
		amr_max_level=0
		amr_max_grid_size="256"
	    elif [ $refinement -eq 2 ]; then
		amr_max_level=1
		amr_ref_ratio="2"
		amr_max_grid_size="256 256"
	    elif [ $refinement -eq 4 ]; then
		amr_max_level=1
		amr_ref_ratio="4"
		amr_max_grid_size="256 256"
	    elif [ $refinement -eq 8 ]; then
		amr_max_level=2
		amr_ref_ratio="4 2"
		amr_max_grid_size="256 256 512"
	    elif [ $refinement -eq 16 ]; then
		amr_max_level=2
		amr_ref_ratio="4 4"
		amr_max_grid_size="256 256 512"
	    elif [ $refinement -eq 32 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 2"
		amr_max_grid_size="256 256 512 512"
	    elif [ $refinement -eq 64 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 4"
		amr_max_grid_size="256 256 512 512"
	    fi

        elif [ $ncell -eq 4096 ]; then

            amr_blocking_factor="128"

	    if [ $refinement -eq 1 ]; then
		amr_max_level=0
		amr_max_grid_size="512"
	    elif [ $refinement -eq 2 ]; then
		amr_max_level=1
		amr_ref_ratio="2"
		amr_max_grid_size="512 512"
	    elif [ $refinement -eq 4 ]; then
		amr_max_level=1
		amr_ref_ratio="4"
		amr_max_grid_size="512 512"
	    elif [ $refinement -eq 8 ]; then
		amr_max_level=2
		amr_ref_ratio="4 2"
		amr_max_grid_size="512 512 1024"
	    elif [ $refinement -eq 16 ]; then
		amr_max_level=2
		amr_ref_ratio="4 4"
		amr_max_grid_size="512 512 1024"
	    elif [ $refinement -eq 32 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 2"
		amr_max_grid_size="512 512 1024 1024"
	    elif [ $refinement -eq 64 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 4"
		amr_max_grid_size="512 512 1024 1024"
	    fi

        elif [ $ncell -eq 8192 ]; then

            amr_blocking_factor="256"

	    if [ $refinement -eq 1 ]; then
		amr_max_level=0
		amr_max_grid_size="1024"
	    elif [ $refinement -eq 2 ]; then
		amr_max_level=1
		amr_ref_ratio="2"
		amr_max_grid_size="1024 1024"
	    elif [ $refinement -eq 4 ]; then
		amr_max_level=1
		amr_ref_ratio="4"
		amr_max_grid_size="1024 1024"
	    elif [ $refinement -eq 8 ]; then
		amr_max_level=2
		amr_ref_ratio="4 2"
		amr_max_grid_size="1024 1024 2048"
	    elif [ $refinement -eq 16 ]; then
		amr_max_level=2
		amr_ref_ratio="4 4"
		amr_max_grid_size="1024 1024 2048"
	    fi

        elif [ $ncell -eq 16384 ]; then

            amr_blocking_factor="512"

	    if [ $refinement -eq 1 ]; then
		amr_max_level=0
		amr_max_grid_size="2048"
	    elif [ $refinement -eq 2 ]; then
		amr_max_level=1
		amr_ref_ratio="2"
		amr_max_grid_size="2048 2048"
	    elif [ $refinement -eq 4 ]; then
		amr_max_level=1
		amr_ref_ratio="4"
		amr_max_grid_size="2048 2048"
	    elif [ $refinement -eq 8 ]; then
		amr_max_level=2
		amr_ref_ratio="4 2"
		amr_max_grid_size="2048 2048 4096"
	    elif [ $refinement -eq 16 ]; then
		amr_max_level=2
		amr_ref_ratio="4 4"
		amr_max_grid_size="2048 2048 4096"
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

# Make a full plotfile every second.

amr_plot_per="1.0"
amr_plot_vars="ALL"
amr_derive_plot_vars="ALL"

# Make small plotfiles rapidly.

amr_small_plot_per="0.01"
amr_small_plot_vars="density Temp"
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
mass_P_default="0.64"
mass_S_default="0.64"
limiter_mode_default="1"
burning_mode_default="1"
T_min_default="1.0e8"
rho_min_default="1.0e6"
small_temp_default="1.0e7"
spec_tol_default="1.0e-6"
enuc_tol_default="1.0e-6"
temp_tol_default="1.0e-6"
grav_tol_default="1.0e-10"
react_T_min_default="1.0e8"
react_rho_min_default="1.0e6"

ncell=$ncell_default
mass_P=$mass_P_default
mass_S=$mass_S_default
castro_dtnuc_e=$dtnuc_e_default
castro_dtnuc_X=$dtnuc_X_default
castro_dxnuc=$dxnuc_default
castro_react_T_min=$react_T_min_default
castro_react_rho_min=$react_rho_min_default
castro_small_temp=$small_temp_default



# The flag we will use to determine whether to run the job.

to_run=1



# Test the effect of the burning resolution limiter.

amr_max_level=9

castro_dxnuc="1.0e-1"

ncell_list="256 512 1024 2048 4096"

for ncell in $ncell_list
do

    # For the higher resolution runs, the tolerance
    # on the gravity solve needs to be loosened a little.

    if   [ $ncell -ge 8192 ]; then
        gravity_abs_tol="1.e-8"
    elif [ $ncell -eq 4096 ]; then
        gravity_abs_tol="1.e-9"
    elif [ $ncell -eq 2048 ]; then
        gravity_abs_tol="5.e-9"
    else
        gravity_abs_tol=$grav_tol_default
    fi

    if   [ $ncell -eq 256 ]; then
        refinement_list="1 2 4 8 16 32 64 128"
    elif [ $ncell -eq 512 ]; then
        refinement_list="1 2 4 8 16 32 64"
    elif [ $ncell -eq 1024 ]; then
        refinement_list="1 2 4 8 16"
    elif [ $ncell -eq 2048 ]; then
        refinement_list="1 2 4 8"
    elif [ $ncell -eq 4096 ]; then
        refinement_list="1"
    fi

    base_dir=$results_dir/amr/n$ncell
    start_dir=$base_dir/start
    save_dir=$base_dir/save

    start_done="0"

    if [ -d $start_dir ]; then
        dir=$start_dir
        start_done=$(is_dir_done)
    fi

    save_done="0"

    if [ -d "$save_dir/chk00000" ] || [ -d "$save_dir/output/chk00000" ]; then
        save_done="1"
    fi

    if [ $start_done -ne 1 ]; then

        # First, we need to do the initial run, up to a point
        # just prior to the detonation. We'll run up until the
        # density reaches 10^7 g/cc.

        castro_density_stopping_criterion="1.e7"

        dir=$start_dir
        set_run_opts
        if [ $to_run -eq 1 ]; then
            run
        fi

        unset castro_density_stopping_criterion

    elif [ $save_done -ne 1 ]; then

        # Now, if we have completed this step, the next step is
        # to relabel this checkpoint as step 0 at time 0.

        checkpoint=$(get_last_checkpoint $start_dir)

        old_stop_time=$stop_time

        amr_regrid_on_restart="1"
        amr_checkpoint_on_restart="1"
        castro_reset_checkpoint_time="0.0"
        castro_reset_checkpoint_step="0"
        stop_time="0.0"
        max_step="0"
        amr_check_int="1" # So no checkpoints get deleted by the run scripts

        dir=$save_dir
        set_run_opts

        # First, create the directory without submitting the job.

        no_submit=1

        if [ $to_run -eq 1 ]; then
            run
        fi

        unset no_submit

        # Now, with the directory created, copy the source checkpoint and do the run.

        if [ ! -d "$save_dir/$checkpoint" ]; then
            echo "Copying final checkpoint to" $save_dir "for relabelling."
            cp -r "$start_dir/$checkpoint" "$save_dir"
            rm -f "$save_dir/$checkpoint/jobIsDone"
        fi

        submit_even_if_done="1"
        if [ $to_run -eq 1 ]; then
            run
        fi
        unset submit_even_if_done

        unset amr_regrid_on_restart
        unset amr_checkpoint_on_restart
        unset castro_reset_checkpoint_time
        unset castro_reset_checkpoint_step
        unset max_step
        unset amr_check_int
        stop_time=$old_stop_time

    else

        # At this point we should have a chk00000 in $save_dir. The final step is to
        # create the run directory (without submitting the job), and copy the checkpoint
        # into it.

        for refinement in $refinement_list
        do

            dir=$base_dir/r$refinement

            job_flag=$(is_job_running $dir)

            if [ $job_flag -ne 1 ]; then

                done_flag=$(is_dir_done)

                if [ $done_flag -ne 1 ]; then

                    no_submit=1

                    set_run_opts
                    if [ $to_run -eq 1 ]; then
                        run
                    fi

                    unset no_submit

                    if [ ! -d "$dir/chk00000" ] && [ ! -d "$dir/output/chk00000" ]; then

                        echo "Copying initial checkpoint to" $dir"."

                        if [ -d "$save_dir/chk00000" ]; then
                            cp -r "$save_dir/chk00000" $dir
                        elif [ -d "$save_dir/output/chk00000" ]; then
                            cp -r "$save_dir/output/chk00000" $dir
                        else
                            echoerr "No initial checkpoint available, exiting."
                            exit
                        fi

                    fi

                    # Before we submit the job, we need to know how long to run.
                    # To figure this out, we'll subtract from the final stop time
                    # the amount already completed.

                    checkpoint=$(get_last_checkpoint $save_dir)
                    chk_time=$(awk 'NR==3' $save_dir/$checkpoint/Header)

                    old_stop_time=$stop_time
                    stop_time=$(echo "$stop_time - $chk_time" | bc -l)

                    replace_inputs_var "stop_time"

                    stop_time=$old_stop_time

                    # We also want to make sure that we get a final plotfile,
                    # since the new stop_time might not be a regular multiple
                    # of the plotting interval.

                    castro_output_at_completion=1

                    replace_inputs_var "castro_output_at_completion"

                    unset castro_output_at_completion

                fi

            fi

            # Now we can submit the job.

            if [ $to_run -eq 1 ]; then
                run
            fi

        done

    fi

done

ncell=$ncell_default
castro_dxnuc=$dxnuc_default
gravity_abs_tol=$grav_tol_default

unset refinement
unset castro_max_dt
unset amr_max_level



# Fix the resolution for the remainder of the runs.

ncell=$ncell_default
amr_max_level=0



# Test the effect of the different timestep limiter methods. For this we
# will only limit the timestep based on internal energy, for comparison
# to previous work.

castro_dtnuc_e="0.3"
dtnuc_mode_list="1 2 3 4"

for castro_dtnuc_mode in $dtnuc_mode_list
do

    dir=$results_dir/burning_limiter_mode/mode$castro_dtnuc_mode

    set_run_opts
    if [ $to_run -eq 1 ]; then
        run
    fi

done

castro_dtnuc_e=$dtnuc_e_default
castro_dtnuc_mode=$limiter_mode_default



# Test the effect of the burning timestep limiter parameter values.

dtnuc_list="10000.0 1000.0 100.0 10.0 1.0 0.5 0.4 0.3 0.2 0.1 0.01 0.001"

for dtnuc in $dtnuc_list
do

    castro_dtnuc_e=$dtnuc

    dir=$results_dir/burning_limiter_e/dt$castro$dtnuc

    set_run_opts
    if [ $to_run -eq 1 ]; then
        run
    fi

done

castro_dtnuc_e=$dtnuc_e_default



for dtnuc in $dtnuc_list
do

    castro_dtnuc_X=$dtnuc

    dir=$results_dir/burning_limiter_X/dt$castro$dtnuc

    set_run_opts
    if [ $to_run -eq 1 ]; then
        run
    fi

done

castro_dtnuc_X=$dtnuc_X_default



# Test the effect of the various burning modes.

burning_mode_list="0 1 2 3"

for burning_mode in $burning_mode_list
do

    dir=$results_dir/burning_mode/$burning_mode

    set_run_opts
    if [ $to_run -eq 1 ]; then
        run
    fi

done

burning_mode=$burning_mode_default



# Test the dependence on the minimum temperature for reactions.

T_min_list="1.0e7 1.0e8"

for castro_react_T_min in $T_min_list
do

    dir=$results_dir/T_min/T$castro_react_T_min

    set_run_opts
    if [ $to_run -eq 1 ]; then
        run
    fi

done

castro_react_T_min=$T_min_default



# Test the dependence on the minimum density for reactions.

rho_min_list="1.0e0 1.0e6"

for castro_react_rho_min in $rho_min_list
do

    dir=$results_dir/rho_min/rho$castro_react_rho_min

    set_run_opts
    if [ $to_run -eq 1 ]; then
        run
    fi

done

castro_react_rho_min=$rho_min_default



# Test the effect of ODE solver tolerance.

tol_list="1.d-12 1.d-10 1.d-8 1.d-6 1.d-4"

for enuc_tol in $tol_list
do

    atol_enuc=$enuc_tol
    rtol_enuc=$enuc_tol

    dir=$results_dir/enuc_tol/tol$enuc_tol

    set_run_opts
    if [ $to_run -eq 1 ]; then
        run
    fi

done

atol_enuc=$enuc_tol_default
rtol_enuc=$enuc_tol_default



for temp_tol in $tol_list
do

    atol_temp=$temp_tol
    rtol_temp=$temp_tol

    dir=$results_dir/temp_tol/tol$temp_tol

    set_run_opts
    if [ $to_run -eq 1 ]; then
        run
    fi

done

atol_temp=$temp_tol_default
rtol_temp=$temp_tol_default



for spec_tol in $tol_list
do

    atol_spec=$spec_tol
    rtol_spec=$spec_tol

    dir=$results_dir/spec_tol/tol$spec_tol

    set_run_opts
    if [ $to_run -eq 1 ]; then
        run
    fi

done

atol_spec=$spec_tol_default
rtol_spec=$spec_tol_default
