source $WDMERGER_HOME/job_scripts/run_utils.sh

# Loop over the resolutions in question

for problem in 1 2 3
do

  for ncell in 16 32 64 128 256
  do

    dir=$results_dir/p$problem/n$ncell
    sed -i "/problem/c problem = $problem" $compile_dir/$probin
    sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell" $compile_dir/$inputs

    if [ $MACHINE == "BLUE_WATERS" ]; then
	if   [ $ncell -eq 16  ]; then
	    nprocs=1
	    walltime=1:00:00
	elif [ $ncell -eq 32  ]; then
	    nprocs=8
	    walltime=1:00:00
	elif [ $ncell -eq 64  ]; then
	    nprocs=16
	    walltime=1:00:00
	elif [ $ncell -eq 128 ]; then
	    nprocs=64
	    walltime=1:00:00
	elif [ $ncell -eq 256 ]; then
	    nprocs=512
	    walltime=1:00:00
	fi
    fi

    run $dir $nprocs $walltime

  done
done
