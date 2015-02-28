source $WDMERGER_HOME/job_scripts/run_utils.sh

# Problem-specific variables

mass_P=" 0.90"
mass_S="-1.00"

# Loop over resolutions

for ncell in 64 128 256
do
  dir=$results_dir/n$ncell
  amr_n_cell="$ncell $ncell $ncell"

  if [ $MACHINE == "BLUE_WATERS" ]; then

    if [ $ncell -eq 32 ]; then
	nprocs=32
	walltime=1:00:00
    elif [ $ncell -eq 64 ]; then
	nprocs=32
	walltime=12:00:00
    elif [ $ncell -eq 128 ]; then
	nprocs=128
	walltime=24:00:00
    fi

  fi

  run $dir $nprocs $walltime
done
