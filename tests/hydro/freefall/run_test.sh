source $WDMERGER_HOME/job_scripts/run_utils.sh

# Problem-specific variables

mass_P=0.90
mass_S=0.81

do_rotation=0
rotational_period=100.0

plot_per=1.0
check_per=1.0

# We can work out the stopping time using the formula
# t_freefall = rotational_period / (4 * sqrt(2))
# bc doesn't have a native way to evaluate fractional exponents like square roots,
# so we perform the trick described here:
# http://www.linuxquestions.org/questions/programming-9/bc-and-exponents-containing-decimals-and-fractions-755260/
# The argument is that a^b = exp(b*ln(a)).
# By calling bc with the -l option we include the math library, 
# which knows how to evaluate exponents with e() and logarithms with l().
# The $$() tells the makefile that we want to expand the following in the shell.
# Of course, I could have just used sqrt(2) = 1.414, but where's the fun in that?

stop_time=$(echo "0.90 * $(rotational_period) / 4.0 / e(0.5*l(2))" | bc -l)

# Loop over the resolutions in question

for ncell in 64 128 256
do
  dir=$results_dir/n$ncell
  n_cell="$ncell $ncell $ncell"

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
