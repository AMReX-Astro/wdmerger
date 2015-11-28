source $WDMERGER_HOME/job_scripts/run_utils.sh

# Problem-specific variables

mass_P=0.90
mass_S=0.75

castro_do_rotation=0
castro_rotational_period=100.0

no_orbital_kick=T

amr_plot_per=1.0
amr_check_per=1.0

castro_do_react=0

# We can work out the stopping time using the formula
# t_freefall = rotational_period / (4 * sqrt(2)).
# We'll stop 90% of the way there because that's about
# when the stars start coming into contact, and the 
# assumption of spherically symmetric stars breaks down.

stop_time=$(echo "0.90 * $castro_rotational_period / (4.0 * sqrt(2))" | bc -l)

dir=$results_dir

if [ $MACHINE == "BLUE_WATERS" ]; then

  nprocs="2048"
  walltime="2:00:00"

elif [ $MACHINE == "LIRED" ]; then

  nprocs="144"
  walltime="12:00:00"

fi

run
