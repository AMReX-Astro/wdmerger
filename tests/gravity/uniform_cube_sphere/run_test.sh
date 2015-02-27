source $WDMERGER_HOME/job_scripts/run_utils.sh

# Problem inputs choices

geometry_prob_lo="-1.6 -1.6 -1.6"
geometry_prob_hi=" 1.6  1.6  1.6"

castro_do_hydro=0
castro_do_rotation=0

castro_small_dens="1.e-10"

max_step=0

gravity_direct_sum_bcs=1

# Loop over problem choices

for p in 1 2 3
do

  # Loop over resolutions

  for ncell in 16 32 64 128 256
  do

    dir=$results_dir/p$problem/n$ncell
    problem=$p
    amr_n_cell="$ncell $ncell $ncell"

    if [ $MACHINE == "BLUE_WATERS" ]; then
	if   [ $ncell -eq 16  ]; then
	    nprocs=32
	    walltime=1:00:00
	elif [ $ncell -eq 32  ]; then
	    nprocs=32
	    walltime=1:00:00
	elif [ $ncell -eq 64  ]; then
	    nprocs=32
	    walltime=1:00:00
	elif [ $ncell -eq 128 ]; then
	    nprocs=128
	    walltime=1:00:00
	elif [ $ncell -eq 256 ]; then
	    nprocs=1024
	    walltime=1:00:00
	fi
    fi

    run $dir $nprocs $walltime

  done
done
