source $WDMERGER_HOME/job_scripts/run_utils.sh

# Problem inputs choices

prob_lo="-1.6 -1.6 -1.6"
prob_hi=" 1.6  1.6  1.6"

do_hydro=0
do_rotation=0

small_dens="1.e-10"

max_step=0

direct_sum_bcs=1

# Loop over problem choices

for p in 1 2 3
do

  # Loop over resolutions

  for ncell in 16 32 64 128 256
  do

    dir=$results_dir/p$problem/n$ncell
    problem=$p
    n_cell="$ncell $ncell $ncell"

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
