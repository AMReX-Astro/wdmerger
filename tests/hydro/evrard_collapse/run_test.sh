source $WDMERGER_HOME/job_scripts/run_utils.sh

# Problem-specific variables

stop_time=15.0

geometry_prob_lo=" -2.56e9 -2.56e9 -2.56e9"
geometry_prob_hi="  2.56e9  2.56e9  2.56e9"

castro_do_rotation=0

# Loop over the resolutions in question

for gs in 1 2 3 4
do
  for ncell in 64 128
  do
    dir=$results_dir/gs$gs/n$ncell

    castro_grav_source_type=$gs
    amr_n_cell="$ncell $ncell $ncell"

    if [ $MACHINE == "BLUE_WATERS" ]; then

      if [ $ncell -eq 64 ]; then
	  nprocs=32
	  walltime=1:00:00
      elif [ $ncell -eq 128 ]; then
	  nprocs=128
	  walltime=2:00:00
      fi

    fi

    run $dir $nprocs $walltime
  done
done
