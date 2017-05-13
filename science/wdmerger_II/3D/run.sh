source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

    if [ -z $ncell ]; then
	echo "ncell not set; exiting."
	exit
    fi

    amr_n_cell="$ncell $ncell $ncell"

    # Define a maximum level for refinement of the 
    # stars and the high-temperature regions.

    max_level="0"

    max_stellar_tagging_level=$max_level
    max_temperature_tagging_level=$max_level

    # If we specify a refinement parameter,
    # we can use that to control amr_ref_ratio.

    # This problem requires some careful thought
    # for load balancing. When we initially start the
    # burn, the burning will be confined to a relatively
    # small region. If we use a too-large max_grid_size,
    # then the work will not be even distributed. The burn
    # strongly dominates the run time cost when we're on 
    # refined levels, and we're intentionally disabling
    # the elliptic solve on those levels, so we can use
    # very small boxes on the refined levels in an effort
    # to distribute the work better.

    if [ ! -z $refinement ]; then
	if [ $refinement -eq 1 ]; then
	    amr_max_level=0
	    amr_max_grid_size="32"
	    amr_blocking_factor="32"
	elif [ $refinement -eq 2 ]; then
	    amr_max_level=1
	    amr_ref_ratio="2"
	    amr_max_grid_size="32 32"
	    amr_blocking_factor="32 8"
	elif [ $refinement -eq 4 ]; then
	    amr_max_level=1
	    amr_ref_ratio="4"
	    amr_max_grid_size="32 32"
	    amr_blocking_factor="32 8"
	elif [ $refinement -eq 8 ]; then
	    amr_max_level=2
	    amr_ref_ratio="4 2"
	    amr_max_grid_size="32 32 48"
	    amr_blocking_factor="32 8 8"
	elif [ $refinement -eq 16 ]; then
	    amr_max_level=2
	    amr_ref_ratio="4 4"
	    amr_max_grid_size="32 32 48"
	    amr_blocking_factor="32 8 8"
	elif [ $refinement -eq 32 ]; then
	    amr_max_level=3
	    amr_ref_ratio="4 4 2"
	    amr_max_grid_size="32 32 48 48"
	    amr_blocking_factor="32 8 8 8"
	elif [ $refinement -eq 64 ]; then
	    amr_max_level=3
	    amr_ref_ratio="4 4 4"
	    amr_max_grid_size="32 32 48 48"
	    amr_blocking_factor="32 8 8 8"
	elif [ $refinement -eq 128 ]; then
	    amr_max_level=4
	    amr_ref_ratio="4 4 4 2"
	    amr_max_grid_size="32 32 48 48 64"
	    amr_blocking_factor="32 8 8 8 8"
	elif [ $refinement -eq 256 ]; then
	    amr_max_level=4
	    amr_ref_ratio="4 4 4 4"
	    amr_max_grid_size="32 32 48 48 64"
	    amr_blocking_factor="32 8 8 8 8"
	elif [ $refinement -eq 512 ]; then
	    amr_max_level=5
	    amr_ref_ratio="4 4 4 4 2"
	    amr_max_grid_size="32 32 48 48 64 128"
	    amr_blocking_factor="32 8 8 8 8 16"
	elif [ $refinement -eq 1024 ]; then
	    amr_max_level=5
	    amr_ref_ratio="4 4 4 4 4"
	    amr_max_grid_size="32 32 48 48 64 128"
	    amr_blocking_factor="32 8 8 8 8 16"
	fi
    fi

    if [ $MACHINE == "TITAN" ]; then

	if [ $ncell -eq 256 ]; then
	    if [ ! -z $refinement ]; then
		if [ $refinement -le 8 ]; then
		    nprocs="512"
		    walltime="2:00:00"
		elif [ $refinement -le 64 ]; then
		    nprocs="2048"
		    walltime="6:00:00"
		elif [ $refinement -le 256 ]; then
		    nprocs="512"
		    walltime="2:00:00"
		else
		    nprocs="2048"
		    walltime="6:00:00"
		fi
	    else
		nprocs="512"
		walltime="2:00:00"
	    fi
	fi

    fi

}

# Needed for the makefile: we want to compile in 3D.

DIM=3

# Grab the right inputs and probin files.

inputs=inputs_3d
probin=probin

# Variables we need to set up the collision.

problem=0
collision_separation=4.0
collision_impact_parameter=0.0
castro_do_rotation=0

# Set the maximum AMR level allowed.
# For this problem, we want to allow the refinement
# based on burning to go as deep as it needs to,
# at least within reason. For amr.max_level = 4,
# the maximum resolution will be about 1.5 km.

max_level_default=9
amr_max_level=$max_level_default

# Make a full plotfile every tenth of a second.

amr_plot_per=0.1

# Save checkpoints every second.

amr_check_per=1.0

# Derive the sound speed, velocity, and pressure for plot files.
# These are useful variables in understanding what is happening 
# at the contact region between the stars.

amr_derive_plot_vars="pressure x_velocity y_velocity z_velocity soundspeed"

# Rapidly print out plots of the density and temperature, for visualization purposes.

amr_small_plot_per="0.05"
amr_small_plot_vars="density Temp"

# Don't allow the timestep to be larger than the small_plot_per value.
# This ensures that we don't have any jumps in the iteration.

castro_max_dt="0.05"

# Initial timestep shortening factor.

castro_init_shrink=1.0

# Set the interval for doing diagnostic global sums.

castro_sum_interval=-1
castro_sum_per=0.001

# This problem determines its own stopping condition,
# so just set a safe upper limit.

stop_time=100.0

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

co_wd_c_frac=0.5
co_wd_o_frac=0.5

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

ncell="256"
castro_dtnuc_e="1.e200" #"1.0e2"
castro_dtnuc_X="1.e200"

mass_P="0.64"
mass_S="0.64"



# Vary the impact parameter and level of refinement.

castro_dxnuc="1.0e-1"

for collision_impact_parameter in 0.0
do

  for refinement in 1 2 4 8 16
  do

      dir=$results_dir/3D/mass_P_$mass_P/mass_S_$mass_S/c$co_wd_c_frac/o$co_wd_o_frac/b$collision_impact_parameter/dxnuc/r$refinement/

      set_run_opts
      run

  done

done

castro_dxnuc="1.0e200"



center_tagging_radius="2.0d8"

for collision_impact_parameter in 0.0
do

  for refinement in 1 2 4 8 16
  do

      dir=$results_dir/3D/mass_P_$mass_P/mass_S_$mass_S/c$co_wd_c_frac/o$co_wd_o_frac/b$collision_impact_parameter/center_tagging/r$refinement/

      set_run_opts
      run

  done

done

center_tagging_radius="0.0d0"

unset refinement



# Test the effect of the impact parameter. We'll use low resolution runs for this.

b_list="0.0 0.01 0.02 0.03 0.04 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75"

for collision_impact_parameter in $b_list
do

    dir=$results_dir/3D/mass_P_$mass_P/mass_S_$mass_S/c$co_wd_c_frac/o$co_wd_o_frac/impact_parameter/b$collision_impact_parameter/

    set_run_opts
    run

done

collision_impact_parameter="0.0"



# Do an unequal mass collision for the purpose of generating a gravitational wave signal.

mass_P=0.80
mass_S=0.60

collision_impact_parameter="0.8"

castro_use_stopping_criterion="0"
stop_time="20.0"

dir=$results_dir/3D/mass_P_$mass_P/mass_S_$mass_S/c$co_wd_c_frac/o$co_wd_o_frac/gw_signal/b$collision_impact_parameter/

set_run_opts
run

unset castro_use_stopping_criterion
stop_time="100.0"
