source $WDMERGER_HOME/job_scripts/run_utils.sh

# Problem-specific variables

geometry_is_periodic="1 1"

geometry_prob_lo="0.0 0.0"
geometry_prob_hi="1.0 1.0"

castro_lo_bc="0 0"
castro_hi_bc="0 0"

# Loop over problem choices

for p in 1 2 3
do

  # Loop over possible bulk flow velocities

  for vel in 0 1 3 10 30 100
  do

    # Loop over the resolutions in question

    for ncell in 64 128 256 1024 2048 4096
    do

      # We only want to do the high-resolution convergence test
      # for the zero-velocity problem.

      if [ $ncell -ge 1024 ] && [ $vel -gt 0 ]; then
	  continue
      fi

      dir=$results_dir/problem$p/velocity$vel/n$ncell

      problem=$p
      bulk_velocity=$vel
      amr_n_cell="$ncell $ncell $ncell"
      
      # Determine stopping time based on problem of interest

      if [ $p -eq 1 ]; then
	  stop_time=2.0
	  amr_plot_per=0.05
      elif [ $p -eq 2 ]; then
	  stop_time=2.0
	  amr_plot_per=0.05
      elif [ $p -eq 3 ]; then
	  stop_time=10.0
	  amr_plot_per=0.1
      fi

      amr_check_per=$amr_plot_per

      # Set number of processors based on amount of work

      if [ $MACHINE == "BLUE_WATERS" ]; then

	if [ $ncell -eq 64 ]; then
	    nprocs=16
	    walltime=00:30:00
	elif [ $ncell -eq 128 ]; then
	    nprocs=128
	    walltime=01:00:00
	elif [ $ncell -eq 256 ]; then
	    nprocs=1024
	    walltime=04:00:00
	elif [ $ncell -eq 512 ]; then
	    nprocs=2048
	    walltime=12:00:00
	elif [ $ncell -eq 1024 ]; then
	    nprocs=128
	    walltime=02:00:00
	elif [ $ncell -eq 2048 ]; then
	   nprocs=1024
	   walltime=04:00:00
	elif [ $ncell -eq 4096 ]; then
	   nprocs=1024
	   walltime=20:00:00
	fi

      fi

      run $dir $nprocs $walltime

    done
  done
done
