source $WDMERGER_HOME/job_scripts/run_utils.sh

TEST_DIR=$CASTRO_DIR/Exec/uniform_cube_sphere

cp $TEST_DIR/inputs source/
cp $TEST_DIR/probin source/

# Loop over problem choices

for problem in 1 2 3
do

  # Loop over resolutions

  for ncell in 16 32 64 128 256
  do

    dir=$results_dir/p$problem/n$ncell
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
