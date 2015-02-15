source $WDMERGER_HOME/job_scripts/run_utils.sh

# Loop over the resolutions in question

for ncell in 64 128 256
do
  dir=$results_dir/n$ncell
  sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell" $compile_dir/$inputs

  if [ $MACHINE == "BLUE_WATERS" ]; then

    if   [ $ncell -eq 64  ]; then
	nprocs="32"
	walltime="2:00:00"
    elif [ $ncell -eq 128 ]; then
	nprocs="128"
	walltime="2:00:00"
    elif [ $ncell -eq 256 ]; then
	nprocs="1024"
	walltime="2:00:00"
    fi

  fi

  run $dir $nprocs $walltime
done
