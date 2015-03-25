# During the run, we'll check whether we're about to run out of time.
# If so, we create a dump_and_stop file. At the end of the next timestep,
# BoxLib will check whether this file exists and terminate the run.
# The argument gives us how long we expect the run to last, in seconds.

if [ ! -z $1 ]; then
    dir=$1
else
    dir=.
fi

# Bring in main helper functions
source $WDMERGER_HOME/job_scripts/run_utils.sh

# Default to cycling every second; we'll 
# update with a smarter choice as we go.

cycle_time=1

# Safety factor: end the run if we're within
# this many timesteps of the total walltime.

safety=20

# Get the name of the output file. There should only be one running.

while true; do

  sleep $cycle_time

  filename=$(find $dir -maxdepth 1 -name "*$run_ext")

  if [ ! -z $filename ]; then
      break
  fi

done

# Extract the job number from the filename.

job_number=${filename#$dir/}
job_number=${job_number%%.*}

# Store the current wall time, in seconds.

start_time=$(date +%s.%N)

# Determine how much time the job has, in seconds.

total_time=$(get_remaining_walltime $job_number)

# Do a continuous loop, but give ourselves an escape hatch
# by creating a particular stop_all file. This may be helpful
# because we're going to background and disown this process
# so we want some way to kill the process in case things go wrong.

while true; do

  if [ -e $WDMERGER_HOME/job_scripts/stop_all ]; then
      exit
  fi

  if [ -e $dir/stop_all ]; then
      exit
  fi

  # Get the median timestep wall time using the last 10 timesteps.

  timestep=$(get_median_timestep $filename 10)

  if [ -z $timestep ]; then
      sleep 1
      continue
  fi

  # Determine the total time remaining in the job, in seconds.

  curr_time=$(date +%s.%N)
  time_elapsed=$(echo "$curr_time - $start_time" | bc -l)
  time_remaining=$(echo "$total_time - $time_elapsed" | bc -l)

  # If we're within $safety steps of the end, kill the run.
  # If not, use the current timestep to determine how long
  # we want to wait before checking again. It should be
  # quite a bit less than the safety factor so that we 
  # guarantee we check before the run terminates.

  to_stop=$(echo "$time_remaining < $safety * $timestep" | bc -l)

  if [ $to_stop -eq 1 ]; then
      touch $dir/dump_and_stop
      exit
  else
      cycle_time=$(echo "$timestep * $safety / 2.0" | bc -l)
  fi

  sleep $cycle_time

done
