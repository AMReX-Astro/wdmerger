source $WDMERGER_HOME/job_scripts/run_utils.sh

# Functions we'll need

function set_run_opts {

  if [ $MACHINE == "LIRED" ]; then
      if   [ $ncell == 128 ]; then
	  nprocs=24
	  walltime=12:00:00
      elif [ $ncell == 256 ]; then
	  nprocs=72
	  walltime=12:00:00
      elif [ $ncell == 512 ]; then
	  nprocs=288
	  walltime=12:00:00
      elif [ $ncell == 1024 ]; then
	  nprocs=576
	  walltime=12:00:00
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

amr_check_per=100.0
amr_plot_per=100.0

amr_small_plot_per=1.0
amr_small_plot_vars="density"

ncell_list="256 512"

castro_do_react=1
castro_react_T_min=1.0e8

mass_P=1.20
mass_S=0.90

stop_time=10000.0

roche_radius_factor=2.0d0

castro_riemann_solver=0

# So that our diagnostic files are not inundated with information,
# print out the global diagnostic sums somewhat sparsely.

castro_sum_interval=-1
castro_sum_per=1.0

# First set up the problem using the approximate ICs.

problem=2
castro_do_rotation=0
castro_do_sponge=0

for ncell in $ncell_list
do

    dir=$results_dir/approximate/n$ncell

    set_run_opts
    run

done


# Now set up the problem using the accurate ICs.

problem=3
castro_do_rotation=1
castro_do_sponge=1
frame_list="1 2"
castro_rotational_dPdt=-1.0

for ncell in $ncell_list
do
    for accurate_IC_frame in $frame_list
    do

	dir=$results_dir/accurate/n$ncell/frame$accurate_IC_frame

	set_run_opts
	run

    done
done
