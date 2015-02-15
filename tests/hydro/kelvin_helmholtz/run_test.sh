source $WDMERGER_HOME/job_scripts/run_utils.sh

# Loop over problem choices

for problem in 1 2 3
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

      dir=$results_dir/problem$problem/velocity$vel/n$ncell

      sed -i "/problem/c problem = $problem" $compile_dir/$probin
      sed -i "/bulk_velocity/c bulk_velocity = $vel" $compile_dir/$probin

      sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell" $compile_dir/$inputs

      # Determine stopping time based on problem of interest

      if [ $problem -eq 1 ]; then
	  stop_time=2.0
	  plot_per=0.05
      elif [ $problem -eq 2 ]; then
	  stop_time=2.0
	  plot_per=0.05
      elif [ $problem -eq 3 ]; then
	  stop_time=10.0
	  plot_per=0.1
      fi

      sed -i "/amr.plot_per/c amr.plot_per = $plot_per" $compile_dir/$inputs
      sed -i "/amr.check_per/c amr.check_per = $plot_per" $compile_dir/$inputs
      sed -i "/stop_time/c stop_time = $stop_time" $compile_dir/$inputs

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
