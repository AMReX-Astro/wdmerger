#!/bin/bash

source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

    if [ -z $ncell ]; then
	echo "ncell not set; exiting."
	exit
    fi

    # Set up the 2D geometry appropriately.

    amr_n_cell="$(echo "$ncell / 2" | bc) $ncell"
    geometry_prob_lo="0.0 $prob_lo"
    geometry_prob_hi="$prob_hi $prob_hi"

    # Define a maximum level for refinement of the 
    # stars and the high-temperature regions.

    amr_max_level="0"

    # If we specify a refinement parameter,
    # we can use that to control amr_ref_ratio.

    if [ ! -z $refinement ]; then

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

    if [ $MACHINE == "LIRED" ]; then

	queue="extended"
        nprocs="24"
        walltime="24:00:00"
        OMP_NUM_THREADS=1

    elif [ $MACHINE == "SEAWULF" ]; then

        queue="long"
        nprocs="28"
        walltime="24:00:00"
        OMP_NUM_THREADS=1

        if [ $ncell -gt 4096 ]; then
            queue="short"
            nprocs="224"
            walltime="4:00:00"
        fi

    fi

}

# Specify the problem directory.

problem_dir=$CASTRO_HOME/Exec/science/wdmerger

# Needed for the makefile: we want to compile in 2D for these tests.

DIM="2"

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
amr_derive_plot_vars="pressure x_velocity y_velocity z_velocity soundspeed"

# Take small plotfiles rapidly.

amr_small_plot_per="0.05"
amr_small_plot_vars="density Temp"

# Save checkpoints every second.

amr_check_per="1.0"

# Initial timestep shortening factor.

castro_init_shrink="0.1"

# Set the interval for doing diagnostic global sums.

castro_sum_interval="1"

# The interesting part of the evolution in this problem finishes by t = 10s.
# Since we know this in advance, we do not need to use the stopping criterion.

stop_time="10.0"
castro_use_stopping_criterion="0"

# Enable reactions.

castro_do_react="1"

# Disable rotation.

castro_do_rotation="0"

# The collision papers in the literature all use an equal C/O ratio 
# by mass in the initial white dwarfs. We will do this too for 
# comparison purposes.

co_wd_c_frac="0.5d0"
co_wd_o_frac="0.5d0"

# Turn off the sponge, it doesn't matter for collisions.

castro_do_sponge="0"

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
network_default="aprox13"
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



# Test the effect of the burning resolution limiter.
# For this test, set the maximum timestep and small_plot_per
# to be somewhat small so that we can make nice movies 
# of temperature and density.

castro_max_dt="0.05"

amr_max_level=9

castro_dxnuc="1.0e-1"

ncell_list="256 512 1024 2048 4096 8192 16384"

for ncell in $ncell_list
do

    # For the higher resolution runs, the tolerance
    # on the gravity solve needs to be loosened a little.

    if   [ $ncell -gt 8192 ]; then
        gravity_abs_tol="1.e-8"
    elif [ $ncell -gt 2048 ]; then
        gravity_abs_tol="1.e-9"
    else
        gravity_abs_tol=$grav_tol_default
    fi

    if   [ $ncell -eq 256 ]; then
        refinement_list="1 2 4 8 16 32 64 128 256"
    elif [ $ncell -eq 512 ]; then
        refinement_list="1 2 4 8 16 32 64 128"
    elif [ $ncell -eq 1024 ]; then
        refinement_list="1 2 4 8 16 32 64"
    elif [ $ncell -eq 2048 ]; then
        refinement_list="1 2 4 8 16 32"
    elif [ $ncell -eq 4096 ]; then
        refinement_list="1 2 4 8 16"
    elif [ $ncell -eq 8192 ]; then
        refinement_list="1 2 4 8"
    elif [ $ncell -eq 16384 ]; then
        refinement_list="1 2"
    fi

    for refinement in $refinement_list
    do

        dir=$results_dir/amr/n$ncell/r$refinement

        set_run_opts
        run

    done
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
    run

done

castro_dtnuc_e=$dtnuc_e_default
castro_dtnuc_mode=$limiter_mode_default



# Test the effect of the burning timestep limiter parameter values.

dtnuc_list="10000.0 5000.0 2000.0 1000.0 500.0 200.0 100.0 50.0 20.0 10.0 5.0 2.0 1.0 0.5 0.4 0.3 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0005 0.0002 0.0001"

for dtnuc in $dtnuc_list
do

    castro_dtnuc_e=$dtnuc

    dir=$results_dir/burning_limiter_e/dt$castro$dtnuc

    set_run_opts
    run

done

castro_dtnuc_e=$dtnuc_e_default



dtnuc_list="10000.0 5000.0 2000.0 1000.0 500.0 200.0 100.0 50.0 20.0 10.0 5.0 2.0 1.0 0.5 0.4 0.3 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001"

for dtnuc in $dtnuc_list
do

    castro_dtnuc_X=$dtnuc

    dir=$results_dir/burning_limiter_X/dt$castro$dtnuc

    set_run_opts
    run

done

castro_dtnuc_X=$dtnuc_X_default



# Test the effect of the various burning modes.

burning_mode_list="0 1 2 3"

for burning_mode in $burning_mode_list
do

    dir=$results_dir/burning_mode/$burning_mode

    set_run_opts
    run

done

burning_mode=$burning_mode_default



# Test the dependence on the minimum temperature for reactions.

T_min_list="2.0e7 4.0e7 6.0e7 8.0e7 1.0e8 2.0e8 3.0e8 4.0e8"

for castro_react_T_min in $T_min_list
do

    dir=$results_dir/T_min/T$castro_react_T_min

    set_run_opts
    run

done

castro_react_T_min=$T_min_default



# Test the dependence on the minimum density for reactions.

rho_min_list="1.0e0 1.0e1 1.0e2 1.0e3 1.0e4 1.0e5 1.0e6"

for castro_react_rho_min in $rho_min_list
do

    dir=$results_dir/rho_min/rho$castro_react_rho_min

    set_run_opts
    run

done

castro_react_rho_min=$rho_min_default



# Test the dependence on the temperature floor.

small_temp_list="1.0e5 1.0e6 1.0e7"

for castro_small_temp in $small_temp_list
do

    dir=$results_dir/small_temp/T$castro_small_temp

    set_run_opts
    run

done

castro_small_temp=$small_temp_default



# Test the effect of the various networks.

network_list="iso7 aprox13 aprox19 aprox21"
local_compile=1

for NETWORK_DIR in $network_list
do

    dir=$results_dir/networks/$NETWORK_DIR

    set_run_opts
    run

done

NETWORK_DIR=$network_default

unset local_compile



# Test the effect of ODE solver tolerance.

for enuc_tol in 1.d-8 1.d-7 1.d-6 1.d-5 1.d-4 1.d-3
do

    atol_enuc=$enuc_tol
    rtol_enuc=$enuc_tol

    dir=$results_dir/enuc_tol/tol$enuc_tol

    set_run_opts
    run

done

atol_enuc=$enuc_tol_default
rtol_enuc=$enuc_tol_default



for temp_tol in 1.d-8 1.d-7 1.d-6 1.d-5 1.d-4 1.d-3
do

    atol_temp=$temp_tol
    rtol_temp=$temp_tol

    dir=$results_dir/temp_tol/tol$temp_tol

    set_run_opts
    run

done

atol_temp=$temp_tol_default
rtol_temp=$temp_tol_default



for spec_tol in 1.d-12 1.d-11 1.d-10 1.d-9 1.d-8 1.d-7 1.d-6 1.d-5 1.d-4
do

    atol_spec=$spec_tol
    rtol_spec=$spec_tol

    dir=$results_dir/spec_tol/tol$spec_tol

    set_run_opts
    run

done

atol_spec=$spec_tol_default
rtol_spec=$spec_tol_default
