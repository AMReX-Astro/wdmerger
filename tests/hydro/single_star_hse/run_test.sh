source $WDMERGER_HOME/job_scripts/run_utils.sh

# Loop over the resolutions in question

for ncell in 32 64 128
do
  dir=$results_dir/n$ncell
  sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell" $compile_dir/$inputs

  if [ $MACHINE == "BLUE_WATERS" ]; then

    if [ $ncell -eq 32 ]; then
	nprocs=16
	walltime=1:00:00
    elif [ $ncell -eq 64 ]; then
	nprocs=16
	walltime=12:00:00
    elif [ $ncell -eq 128 ]; then
	nprocs=64
	walltime=24:00:00
    fi

  fi

  run $dir $nprocs $walltime
done
