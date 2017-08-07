source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

    if [ -z $ncell ]; then
	echo "ncell not set; exiting."
	exit
    fi

    if [ $DIM -eq 2 ]; then
        amr_n_cell="$(echo "$ncell / 2" | bc) $ncell"
        geometry_prob_lo="0.0 $prob_lo"
        geometry_prob_hi="$prob_hi $prob_hi"
    else
        amr_n_cell="$ncell $ncell $ncell"
        geometry_prob_lo="$prob_lo $prob_lo $prob_lo"
        geometry_prob_hi="$prob_hi $prob_hi $prob_hi"
    fi

    # Define a maximum level for refinement of the 
    # stars and the high-temperature regions.

    max_level="0"

    if [ $ncell -eq 256 ]; then
	if [ -z $refinement ]; then
	    amr_max_level="0"
	fi
    elif [ $ncell -eq 512 ]; then
	if [ -z $refinement ]; then
	    amr_max_level="1"
	fi
	max_level="1"
	amr_ref_ratio="2 4 4 4 4 4 4 4 4"
    elif [ $ncell -eq 1024 ]; then
	if [ -z $refinement ]; then
	    amr_max_level="1"
	fi
	max_level="1"
	amr_ref_ratio="4 4 4 4 4 4 4 4 4"
    elif [ $ncell -eq 2048 ]; then
	if [ -z $refinement ]; then
	    amr_max_level="2"
	fi
	max_level="2"
	amr_ref_ratio="4 2 4 4 4 4 4 4 4"
    elif [ $ncell -eq 4096 ]; then
	if [ -z $refinement ]; then
	    amr_max_level="2"
	fi
	max_level="2"
	amr_ref_ratio="4 4 4 4 4 4 4 4 4"
    elif [ $ncell -eq 8192 ]; then
	if [ -z $refinement ]; then
	    amr_max_level="3"
	fi
	max_level="3"
	amr_ref_ratio="4 4 2 4 4 4 4 4 4"
    elif [ $ncell -eq 16384 ]; then
	if [ -z $refinement ]; then
	    amr_max_level="3"
	fi
	max_level="3"
	amr_ref_ratio="4 4 4 4 4 4 4 4 4"
    fi

    max_stellar_tagging_level=$max_level
    max_temperature_tagging_level=$max_level

    # If we specify a refinement parameter,
    # we can use that to control amr_ref_ratio.

    if [ ! -z $refinement ]; then

	if [ $ncell -eq 256 ]; then

	    if [ $refinement -eq 1 ]; then
		amr_max_level=0
		amr_blocking_factor="8"
		amr_max_grid_size="32"
	    elif [ $refinement -eq 2 ]; then
		amr_max_level=1
		amr_ref_ratio="2"
		amr_blocking_factor="8 8"
		amr_max_grid_size="32 32"
	    elif [ $refinement -eq 4 ]; then
		amr_max_level=1
		amr_ref_ratio="4"
		amr_blocking_factor="8 8"
		amr_max_grid_size="32 32"
	    elif [ $refinement -eq 8 ]; then
		amr_max_level=2
		amr_ref_ratio="4 2"
		amr_blocking_factor="8 8 8"
		amr_max_grid_size="32 32 48"
	    elif [ $refinement -eq 16 ]; then
		amr_max_level=2
		amr_ref_ratio="4 4"
		amr_blocking_factor="8 8 8"
		amr_max_grid_size="32 32 48"
	    elif [ $refinement -eq 32 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 2"
		amr_blocking_factor="8 8 8 8"
		amr_max_grid_size="32 32 48 64"
	    elif [ $refinement -eq 64 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 4"
		amr_blocking_factor="8 8 8 8"
		amr_max_grid_size="32 32 48 64"
	    elif [ $refinement -eq 128 ]; then
		amr_max_level=4
		amr_ref_ratio="4 4 4 2"
		amr_blocking_factor="8 8 8 8 8"
		amr_max_grid_size="32 32 48 64 64"
	    elif [ $refinement -eq 256 ]; then
		amr_max_level=4
		amr_ref_ratio="4 4 4 4"
		amr_blocking_factor="8 8 8 8 8"
		amr_max_grid_size="32 32 48 64 64"
	    elif [ $refinement -eq 512 ]; then
		amr_max_level=5
		amr_ref_ratio="4 4 4 4 2"
		amr_blocking_factor="8 8 8 8 8 8"
		amr_max_grid_size="32 32 48 64 64 128"
	    elif [ $refinement -eq 1024 ]; then
		amr_max_level=5
		amr_ref_ratio="4 4 4 4 2"
		amr_blocking_factor="8 8 8 8 8 8"
		amr_max_grid_size="32 32 48 64 64 128"
	    fi

	elif [ $ncell -eq 1024 ]; then

	    if [ $refinement -eq 1 ]; then
		amr_max_level=1
		amr_ref_ratio="4"
		amr_blocking_factor="8 8"
		amr_max_grid_size="32 32"
	    elif [ $refinement -eq 2 ]; then
		amr_max_level=2
		amr_ref_ratio="4 2"
		amr_blocking_factor="8 8 8"
		amr_max_grid_size="32 32 48"
	    elif [ $refinement -eq 4 ]; then
		amr_max_level=2
		amr_ref_ratio="4 4"
		amr_blocking_factor="8 8 8"
		amr_max_grid_size="32 32 48"
	    elif [ $refinement -eq 8 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 2"
		amr_blocking_factor="8 8 8 8"
		amr_max_grid_size="32 32 48 64"
	    elif [ $refinement -eq 16 ]; then
		amr_max_level=3
		amr_ref_ratio="4 4 4"
		amr_blocking_factor="8 8 8 8"
		amr_max_grid_size="32 32 48 64"
	    elif [ $refinement -eq 32 ]; then
		amr_max_level=4
		amr_ref_ratio="4 4 4 2"
		amr_blocking_factor="8 8 8 8 8"
		amr_max_grid_size="32 32 48 64 128"
	    elif [ $refinement -eq 64 ]; then
		amr_max_level=4
		amr_ref_ratio="4 4 4 4"
		amr_blocking_factor="8 8 8 8 8"
		amr_max_grid_size="32 32 48 64 128"
	    fi

	fi

    fi

    if [ $MACHINE == "LIRED" ]; then

	queue="extended"
        nprocs="24"
        walltime="24:00:00"
        OMP_NUM_THREADS=1

        if [ ! -z $refinement ]; then
            if [ $refinement -gt 256 ]; then
                queue="medium"
                nprocs="192"
                walltime="12:00:00"
                OMP_NUM_THREADS=8
            fi
        fi

        if [ $ncell -gt 1024 ]; then
            queue="medium"
            nprocs="192"
            walltime="12:00:00"
            OMP_NUM_THREADS=8
        fi

    fi

}

# Specify the problem directory.

problem_dir=$CASTRO_HOME/Exec/science/wdmerger

# Needed for the makefile: we want to compile in 2D for most of the tests.

DIM=2

# Get the right inputs and probin files.

inputs=inputs_2d
probin=probin

# Use a narrower domain than the usual default.
# This problem is insensitive to accuracy in the boundary conditions.

prob_lo=-4e9
prob_hi=4e9

# Variables we need to set up the collision.

problem=0
collision_separation=4.0
collision_impact_parameter=0.0
castro_do_rotation=0

# Make a full plotfile every tenth of a second.

amr_plot_per=0.1

# Save checkpoints every second.

amr_check_per=1.0

# Derive the sound speed, velocity, and pressure for plot files.
# These are useful variables in understanding what is happening 
# at the contact region between the stars.

amr_derive_plot_vars="pressure x_velocity y_velocity z_velocity soundspeed"

# Initial timestep shortening factor.

castro_init_shrink=0.01

# Set the interval for doing diagnostic global sums.

castro_sum_interval=1

# This problem self-terminates, so just pick a reasonably safe number here.

stop_time=20.0

# Enable reactions, but for efficiency disable all burning
# for T < 10^8 K and rho < 10^6 g/cc.

castro_do_react=1
castro_react_T_min=1.0e8
castro_react_rho_min=1.0e6

# We want a relatively high small_temp for this problem.

castro_small_temp=1.0e7

# The collision papers in the literature all use an equal C/O ratio 
# by mass in the initial white dwarfs. We will do this too for 
# comparison purposes.

co_wd_c_frac=0.5d0
co_wd_o_frac=0.5d0

# Turn off the sponge, it doesn't matter for collisions.

castro_do_sponge=0

# The timesteps can get quite small if you're fully 
# resolving the burning, so allow for this.

castro_dt_cutoff=1.0e-12

# Allow the timestep to change by up to 5% per advance.
# We don't need to be super strict on the timestep control
# for this problem, we have empirically found.

castro_change_max=1.05

# Some defaults.

ncell_default="256"

# Empirically we have found that the answer is qualitatively
# converged for dtnuc_e < 100. So let's set it here to get 
# enough timestep control to get a reasonable answer, but not
# so low that it takes forever for the runs with multiple levels
# of refinement to finish.

dtnuc_e_default="1.e2"
dtnuc_X_default="1.e200"
dxnuc_default="1.0e200"
mass_P_default="0.64"
mass_S_default="0.64"
network_default="aprox13"
limiter_mode_default="1"
burning_mode_default="1"
co_wd_c_frac_default="0.5"
co_wd_o_frac_default="0.5"
T_min_default="1.0e8"
rho_min_default="1.0e6"
small_temp_default="1.0e7"
spec_tol_default="1.0e-8"
enuc_tol_default="1.0e-6"
temp_tol_default="1.0e-6"

ncell=$ncell_default
mass_P=$mass_P_default
mass_S=$mass_S_default
castro_dtnuc_e=$dtnuc_e_default
castro_dtnuc_X=$dtnuc_X_default
castro_dxnuc=$dxnuc_default



# Test the effect of the burning resolution limiter.
# For this test, set the maximum timestep and small_plot_per
# to be somewhat small so that we can make nice movies 
# of temperature and density.

amr_small_plot_per="0.05"
amr_small_plot_vars="density Temp"
castro_max_dt="0.05"

# For this problem, we want to allow the refinement
# based on burning to go as deep as it needs to,
# at least within reason. For amr.max_level = 4,
# the maximum resolution will be about 1.5 km.

amr_max_level=9

castro_dxnuc="1.0e-1"
castro_dtnuc_e="1.e200"
castro_dtnuc_X="1.e200"

for refinement in 1 2 4 8 16 32 64 128 256 512 1024
do

    dir=$results_dir/2D/dxnuc/r$refinement

    set_run_opts
    run

done

castro_dxnuc=$dxnuc_default



# Test the effect of the resolution based on distance from the center.
# We can use the same scheme as above, but turn off the dxnuc tagging.

center_tagging_radius=2.0d8

for refinement in 1 2 4 8 16 32 64 128 256
do

    dir=$results_dir/2D/center_tagging/r$refinement

    set_run_opts
    run

done

castro_dtnuc_e=$dtnuc_e_default
castro_dtnuc_X=$dtnuc_X_default

unset center_tagging_radius

unset refinement
unset castro_max_dt
unset amr_max_level



# Do full refinement on the stars.

for ncell in 256 512 1024 2048 4096
do

  dir=$results_dir/2D/stellar_tagging/n$ncell

  set_run_opts
  run

done

ncell=$ncell_default



# Fix the resolution for the remainder of the runs.

ncell=$ncell_default
amr_max_level=0

# For the remainder of the runs we do not need to plot frequently.

amr_plot_per=1.0
amr_small_plot_per=-1.0
amr_check_per=1.0



# Test the effect of the different timestep limiter methods. For this we
# will only limit the timestep based on internal energy, for comparison
# to previous work.

castro_dtnuc_e="0.3"
dtnuc_mode_list="1 2 3 4"

for castro_dtnuc_mode in $dtnuc_mode_list
do

    dir=$results_dir/2D/burning_limiter_mode/mode$castro_dtnuc_mode

    set_run_opts
    run

done

castro_dtnuc_e=$dtnuc_e_default
castro_dtnuc_mode=$limiter_mode_default



# Test the effect of the burning timestep limiter parameter values.

castro_dtnuc_X="1.e200"

dtnuc_list="10000.0 1000.0 100.0 10.0 5.0 2.0 1.0 0.5 0.4 0.3 0.2 0.1 0.01 0.001"

for dtnuc in $dtnuc_list
do

    castro_dtnuc_e=$dtnuc

    dir=$results_dir/2D/burning_limiter_e/dt$castro$dtnuc

    set_run_opts
    run

done

castro_dtnuc_e=$dtnuc_e_default
castro_dtnuc_X=$dtnuc_X_default



castro_dtnuc_e="1.e200"

dtnuc_list="10000.0 1000.0 100.0 10.0 5.0 2.0 1.0 0.5 0.4 0.3 0.2 0.1 0.01 0.001"

for dtnuc in $dtnuc_list
do

    castro_dtnuc_X=$dtnuc

    dir=$results_dir/2D/burning_limiter_X/dt$castro$dtnuc

    set_run_opts
    run

done

castro_dtnuc_X=$dtnuc_X_default
castro_dtnuc_e=$dtnuc_e_default



# Test the effect of the various burning modes.

burning_mode_list="0 1 2 3"

for burning_mode in $burning_mode_list
do

    dir=$results_dir/2D/burning_mode/$burning_mode

    set_run_opts
    run

done

burning_mode=$burning_mode_default



# Test the dependence on the C/O mass fraction.

c_frac_list=(0.30 0.40 0.50 0.60 0.70)
o_frac_list=(0.70 0.60 0.50 0.40 0.30)

list_length=${#c_frac_list[@]}

for index in $(seq 0 $(($list_length-1)))
do

    co_wd_c_frac=${c_frac_list[$index]}
    co_wd_o_frac=${o_frac_list[$index]}

    dir=$results_dir/2D/co/c"$co_wd_c_frac"o"$co_wd_o_frac"

    set_run_opts
    run

done

co_wd_c_frac=$co_wd_c_frac_default
co_wd_o_frac=$co_wd_o_frac_default



# Test the dependence on the minimum temperature for reactions.

T_min_list="2.0e7 4.0e7 6.0e7 8.0e7 1.0e8 2.0e8 3.0e8 4.0e8"

for castro_react_T_min in $T_min_list
do

    dir=$results_dir/2D/T_min/T$castro_react_T_min

    set_run_opts
    run

done

castro_react_T_min=$T_min_default



# Test the dependence on the minimum density for reactions.

rho_min_list="1.0e0 1.0e1 1.0e2 1.0e3 1.0e4 1.0e5 1.0e6"

for castro_react_rho_min in $rho_min_list
do

    dir=$results_dir/2D/rho_min/rho$castro_react_rho_min

    set_run_opts
    run

done

castro_react_rho_min=$rho_min_default



# Test the dependence on the temperature floor.

small_temp_list="1.0e5 1.0e6 1.0e7"

for castro_small_temp in $small_temp_list
do

    dir=$results_dir/2D/small_temp/T$castro_small_temp

    set_run_opts
    run

done

castro_small_temp=$small_temp_default



# Test the dependence on the dual energy criteria.

eta2_list="0.00 0.0001 0.001 0.01 0.05 0.10"

for castro_dual_energy_eta2 in $eta2_list
do

    dir=$results_dir/2D/eta2/eta$castro_dual_energy_eta2

    set_run_opts
    run

done

unset castro_dual_energy_eta2

eta3_list="0.00 0.0001 0.001 0.01 0.05 0.10"

for castro_dual_energy_eta3 in $eta3_list
do

    dir=$results_dir/2D/eta3/eta$castro_dual_energy_eta3

    set_run_opts
    run

done

unset castro_dual_energy_eta3



# Test the effect of the various networks.

network_list="iso7 aprox19 aprox21 aprox13"
local_compile=1

for Network_dir in $network_list
do

    dir=$results_dir/2D/networks/$Network_dir

    set_run_opts
    run

done

Network_dir=$network_default

unset local_compile



# Test the effect of ODE solver tolerance.

for enuc_tol in 1.d-8 1.d-7 1.d-6 1.d-5 1.d-4 1.d-3
do

    atol_enuc=$enuc_tol
    rtol_enuc=$enuc_tol

    dir=$results_dir/2D/enuc_tol/tol$enuc_tol

    set_run_opts
    run

done

atol_enuc=$enuc_tol_default
rtol_enuc=$enuc_tol_default



for temp_tol in 1.d-8 1.d-7 1.d-6 1.d-5 1.d-4 1.d-3
do

    atol_temp=$temp_tol
    rtol_temp=$temp_tol

    dir=$results_dir/2D/temp_tol/tol$temp_tol

    set_run_opts
    run

done

atol_temp=$temp_tol_default
rtol_temp=$temp_tol_default



for spec_tol in 1.d-12 1.d-11 1.d-10 1.d-9 1.d-8 1.d-7 1.d-6 1.d-5 1.d-4
do

    atol_spec=$spec_tol
    rtol_spec=$spec_tol

    dir=$results_dir/2D/spec_tol/tol$spec_tol

    set_run_opts
    run

done

atol_spec=$spec_tol_default
rtol_spec=$spec_tol_default



# Do a run without any burning. Let this one run for a while
# so that we know what happens to the collision on longer timescales.

castro_use_stopping_criterion=0

for castro_do_react in 0 1
do

    dir=$results_dir/2D/reactions/react$castro_do_react

    set_run_opts
    run

done

unset castro_use_stopping_criterion
castro_do_react=1



# Switch to a 3D test.

DIM=3
inputs=inputs_3d

castro_use_stopping_criterion=0

for castro_do_react in 0 1
do

    dir=$results_dir/3D/reactions/react$castro_do_react

    set_run_opts
    run

done

unset castro_use_stopping_criterion
castro_do_react=1
