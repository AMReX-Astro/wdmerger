source $WDMERGER_HOME/job_scripts/run_utils.sh

# Functions we'll need

function set_run_opts {

  if [ $MACHINE == "BLUE_WATERS" ]; then
      if   [ $ncell == 128 ]; then
	  nprocs=128
	  walltime=12:00:00
      elif [ $ncell == 256 ]; then
	  nprocs=1024
	  walltime=24:00:00
      elif [ $ncell == 512 ]; then
	  nprocs=1024
	  walltime=24:00:00
      elif [ $ncell == 1024 ]; then
	  nprocs=2048
	  walltime=24:00:00
      elif [ $ncell == 2048 ]; then
	  nprocs=4096
	  walltime=24:00:00
      fi
  fi

  if [ ! -z $ncell ]; then

    if [ $ncell -eq 256 ]; then
	amr_n_cell="256 256 256"
	amr_max_level="0"
	amr_max_grid_size="32"
    elif [ $ncell -eq 512 ]; then
	amr_n_cell="256 256 256"
	amr_max_level="1"
	amr_ref_ratio="2"
	amr_max_grid_size="32 32"
    elif [ $ncell -eq 1024 ]; then
	amr_n_cell="256 256 256"
	amr_max_level="1"
	amr_ref_ratio="4"
	amr_max_grid_size="32 48"
    elif [ $ncell -eq 2048 ]; then
	amr_n_cell="256 256 256"
	amr_max_level="2"
	amr_ref_ratio="4 2"
	amr_max_grid_size="32 32 64"
    fi

  fi

}


# Some global variables we'll need

castro_rotational_period=100.0

amr_check_per=20.0
amr_plot_per=20.0

# Main test of equal and unequal mass binaries for multiple orders.

ncell=256
num_periods=25
stop_time=$(echo "$num_periods * $castro_rotational_period" | bc -l)

castro_grav_source_type=4
castro_rot_source_type=4

N_iters=$num_periods

for ratio in equal unequal
do

  if [ $ratio == "equal" ]; then
      mass_P=0.90
      mass_S=0.90
  elif [ $ratio == "unequal" ]; then
      mass_P=0.90
      mass_S=0.60
  fi

  for castro_do_rotation in 0 1
  do

    dir=$results_dir/$ratio/rot$castro_do_rotation

    set_run_opts
    chain

  done
done

# Fixed parameters for the remainder of the tests

mass_S=0.90
mass_P=0.60
ncell=256

stop_time=$castro_rotational_period

# Now do a single orbit test of the various gravity and rotation source terms.

for castro_grav_source_type in 1 2 3 4
do
  for castro_rot_source_type in 0 1 2 3 4
  do

    dir=$results_dir/sources/gs$castro_grav_source_type/rs$castro_rot_source_type

    if [ $castro_rot_source_type == 0 ]; then
        castro_do_rotation=0
    else
        castro_do_rotation=1
    fi

    set_run_opts
    run

  done
done

castro_grav_source_type=4
castro_rot_source_type=4

# A test of the gravity boundary conditions

for castro_do_rotation in 0 1
do
  for gravity_max_multipole_order in 0 2 4 6 8 10 12 14 16
  do

    dir=$results_dir/gravity_bcs/rot$castro_do_rotation/order$gravity_max_multipole_order

    set_run_opts
    run
  
  done
done

# Do a spatial resolution test in both reference frames

for castro_do_rotation in 0 1
do
  for ncell in 256 512 1024
  do

    dir=$results_dir/spatial_convergence/rot$castro_do_rotation/n$ncell

    set_run_opts
    run

  done
done

# Do a time resolution test

ncell=256
castro_init_shrink=1.0

for castro_do_rotation in 0 1
do
  for castro_fixed_dt in 0.01 0.005 0.0025 0.00125
  do

    dir=$results_dir/time_convergence/rot$castro_do_rotation/dt$castro_fixed_dt

    set_run_opts
    run

  done
done
