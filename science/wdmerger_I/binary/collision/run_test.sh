source $WDMERGER_HOME/job_scripts/run_utils.sh

# Problem-specific variables

mass_P=0.64
mass_S=0.64

castro_do_rotation=0

collision=T
collision_separation=4.0

amr_plot_per=0.1
amr_check_per=0.1

stop_time=10.0

# Loop over the resolutions in question

for ncell in 256
do
  dir=$results_dir/n$ncell
  amr_n_cell="$ncell $ncell $ncell"

  if [ $MACHINE == "BLUE_WATERS" ]; then

    if [ $ncell -eq 256 ]; then
	nprocs="2048"
	walltime="2:00:00"
    fi

  fi

  run $dir $nprocs $walltime
done
