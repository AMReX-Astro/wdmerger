source $WDMERGER_HOME/job_scripts/run_utils.sh

# Problem-specific variables

mass_P=0.90
mass_S=0.81

castro_do_rotation=0
castro_rotational_period=100.0

no_orbital_kick=T

amr_plot_per=1.0
amr_check_per=1.0

# We can work out the stopping time using the formula
# t_freefall = rotational_period / (4 * sqrt(2)).
# We'll stop 90% of the way there because that's about
# when the stars start coming into contact, and the 
# assumption of spherically symmetric stars breaks down.

stop_time=$(echo "0.90 * $(rotational_period) / (4.0 * sqrt(2))" | bc -l)

# Loop over the resolutions in question

for ncell in 64 128 256
do
  dir=$results_dir/n$ncell
  amr_n_cell="$ncell $ncell $ncell"

  if [ $MACHINE == "BLUE_WATERS" ]; then

    if   [ $ncell -eq 64  ]; then
	nprocs="64"
	walltime="2:00:00"
    elif [ $ncell -eq 128 ]; then
	nprocs="256"
	walltime="2:00:00"
    elif [ $ncell -eq 256 ]; then
	nprocs="2048"
	walltime="2:00:00"
    fi

  fi

  run $dir $nprocs $walltime
done
