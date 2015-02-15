source $WDMERGER_HOME/job_scripts/run_utils.sh

# Loop over the resolutions in question

for gst in 1 2 3 4
do
  for ncell in 64 128
  do
    dir=$results_dir/gs$gst/n$ncell

    sed -i "/castro.grav_source_type/c castro.grav_source_type = $gst" $compile_dir/$inputs
    sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell" $compile_dir/$inputs

    if [ $MACHINE == "BLUE_WATERS" ]; then

      if [ $ncell -eq 64 ]; then
	  nprocs=16
	  walltime=1:00:00
      elif [ $ncell -eq 128 ]; then
	  nprocs=64
	  walltime=2:00:00
      fi

    fi

    run $dir $nprocs $walltime
  done
done
