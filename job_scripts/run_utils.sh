# run_utils.sh: helper functions for job submission scripts.

# There are several scripts that this file relies upon. The 
# home location for these is in $WDMERGER_HOME/job_scripts.
# However, we make local copies of these files for run 
# directories, so test first if we have a local job_scripts 
# subdirectory, and check there first.

if [ -d "job_scripts" ]; then
  script_dir="job_scripts"
else
  script_dir="$WDMERGER_HOME/job_scripts"
fi

# Functions for working with inputs and probin files.
source $script_dir/inputs.sh
source $script_dir/probin.sh

# Basic mathematical routines.
source $script_dir/math.sh

# System information.
source $script_dir/machines.sh



# Return a string to append to the 'make' command
# if we have various makefile options specified.

function compile_options {

  compile_opts=''

  if [ ! -z $CASTRO_HOME ]; then
      compile_opts=$compile_opts' CASTRO_HOME='$CASTRO_HOME
  fi

  if [ ! -z $DIM ]; then
      compile_opts=$compile_opts' DIM='$DIM
  fi

  if [ ! -z $Network_dir ]; then
      compile_opts=$compile_opts' Network_dir='$Network_dir
  fi

  echo $compile_opts

}



# This uses the functionality built into the CASTRO makefile setup,
# where make print-$VAR finds the variable VAR in the makefile
# variable list and prints it out to stdout. It is the last word
# on the last line of the make output.

function get_make_var {

  make print-$1 -C $compile_dir $(compile_options) &> temp_compile.out 
  cat temp_compile.out | tail -2 | head -1 | awk '{ print $NF }'
  rm -f temp_compile.out

}



# Use the batch submission system to obtain the set of currently submitted jobs.
# It fills the variable job_list

function get_submitted_jobs {

  job_file=$WDMERGER_HOME/job_scripts/jobs.txt

  # Store the result of the queue information command, 
  # but only keep the lines that store jobs with our 
  # username so that we remove all formatting output.
  # We are only interested in storing the job
  # numbers that are somewhere in the queue.

  if [ $batch_system == "PBS" ]; then

      showq -u $USER | grep $USER | awk '{ print $1 }' > $job_file

  elif [ $batch_system == "COBALT" ]; then

      qstat -u $USER | grep $USER | awk '{ print $1 }' > $job_file

  fi

  # Store it as an array.

  if [ -e $job_file ]; then

      num_jobs=$(cat $job_file | wc -l)
      
      for i in $(seq 0 $(($num_jobs-1)))
      do
	  line=$(awk "NR == $i+1" $job_file)
	  job_arr[$i]=$(echo $line | awk '{ print $1 }')
      done

  fi

}



# Given a directory as the first argument, return the numerically last output file.
# Assumes that all output files have the same number of digits; this should work 
# in general except in rare cases.

function get_last_output {

  if [ -z $1 ]; then
      echo "No directory passed to get_last_output; exiting."
      return
  else
      dir=$1
  fi

  # First try looking for jobs that are currently running.

  output=$(find $dir -name "*$run_ext" | sort | tail -1)

  # If that fails, look through completed jobs.

  if [ -z $output ]; then
      output=$(find $dir -name "$job_name*" | sort | tail -1)
  fi

  # Extract out the search directory from the result.

  output=$(echo ${output#$dir/})
  output=$(echo ${output#$dir})

  echo $output

}



# Given a directory as the first argument, return the numerically last checkpoint file.

function get_last_checkpoint {

  if [ -z $1 ]; then
      echo "No directory passed to get_last_checkpoint; exiting."
      return
  else
      dir=$1
  fi

  # Optionally pass in a second argument corresponding to how far back we want to go.
  # By default we return the very latest one.

  if [ ! -z $2 ]; then
      num=$2
  else
      num=1
  fi

  # Doing a search this way will treat first any checkpoint files 
  # with seven digits, and then will fall back to ones with six and then five digits.
  # We want to be smart about this and list the ones in the current directory first,
  # before checking any output directories where the data is archived, because
  # the former are the most likely to be recently created checkpoints.
  # Bash can be pretty slow when looping through a large number of variables,
  # so the added complication is justified.

  checkpointList=""
  checkpointNums=""

  checkpointList+=" $(find $dir -maxdepth 1 -type d -name "*chk???????" | sort -r)"
  checkpointList+=" $(find $dir -mindepth 2 -type d -name "*chk???????" | sort -r)"
  checkpointNums+=" $(find $dir -maxdepth 1 -type d -name "*chk???????" | awk -F/ '{ print $NF }' | sort -r)"
  checkpointNums+=" $(find $dir -mindepth 2 -type d -name "*chk???????" | awk -F/ '{ print $NF }' | sort -r)"

  checkpointList+=" $(find $dir -maxdepth 1 -type d -name "*chk??????" | sort -r)"
  checkpointList+=" $(find $dir -mindepth 2 -type d -name "*chk??????" | sort -r)"
  checkpointNums+=" $(find $dir -maxdepth 1 -type d -name "*chk??????" | awk -F/ '{ print $NF }' | sort -r)"
  checkpointNums+=" $(find $dir -mindepth 2 -type d -name "*chk??????" | awk -F/ '{ print $NF }' | sort -r)"

  checkpointList+=" $(find $dir -maxdepth 1 -type d -name "*chk?????" | sort -r)"
  checkpointList+=" $(find $dir -mindepth 2 -type d -name "*chk?????" | sort -r)"
  checkpointNums+=" $(find $dir -maxdepth 1 -type d -name "*chk?????" | awk -F/ '{ print $NF }' | sort -r)"
  checkpointNums+=" $(find $dir -mindepth 2 -type d -name "*chk?????" | awk -F/ '{ print $NF }' | sort -r)"

  if [ -z "$checkpointList" ]; then

      echo ""
      return

  fi

  # Match up the last checkpoint number with the actual file path location.

  numSkipped=1

  for chkNum in $checkpointNums
  do

      if [ $numSkipped -lt $num ]; then
	  numSkipped=$((numSkipped + 1))
	  continue
      fi

      for chkFile in $checkpointList
      do
	  currBaseName=$(echo $chkFile | awk -F/ '{ print $NF }')

	  if [ "$currBaseName" == "$chkNum" ]; then

	      # The Header is the last thing written -- check if it's there, otherwise,
	      # we can skip this iteration, because it means the latest checkpoint file 
	      # is still being written.

	      if [ -f ${chkFile}/Header ]; then
		  checkpoint=$chkFile
		  break
	      fi
	  fi
      done

      if [ ! -z $checkpoint ]; then
	  break
      fi

  done

  # Extract out the search directory from the result.

  checkpoint=$(echo ${checkpoint#$dir/})
  checkpoint=$(echo ${checkpoint#$dir})

  echo $checkpoint

}



# Obtain the median coarse timestep length from a given output file.

function get_median_timestep {

  # First argument is the name of the file.

  if [ -z $1 ]; then
      # No file was passed to the function, so exit.
      return
  else
      file=$1
  fi

  # Second argument is the number of most recent timesteps to use.
  # If it doesn't exist, default to using all timesteps.

  if [ -z $2 ]; then
      nsteps=-1
  else
      nsteps=$2
  fi

  # Use grep to get all lines containing the coarse timestep time;
  # then, use awk to extract the actual times.

  if [ $nsteps -gt 0 ]; then
      timesteps=$(grep "Coarse" $file | awk -F "Coarse TimeStep time: " '{ print $2 }' | tail -$nsteps)
  else
      timesteps=$(grep "Coarse" $file | awk -F "Coarse TimeStep time: " '{ print $2 }')
  fi

  # Calculate the median.
    
  median_timestep=$(median "$timesteps")

  echo $median_timestep

}



# Predict the next timestep wall time based on the trend of recent steps.

function predict_next_timestep {

  # First argument is the name of the file.

  if [ -z $1 ]; then

      # No file was passed to the function, so exit.

      echo "-1"
      return

  else

      file=$1

  fi

  # Second argument is the number of most recent timesteps to use.
  # If it doesn't exist, default to using the last ten timesteps.

  if [ -z $2 ]; then
      nsteps=10
  else
      nsteps=$2
  fi

  # Use grep to get all lines containing the coarse timestep time;
  # then, use awk to extract the actual times.

  timesteps=$(grep "Coarse" $file | awk -F "Coarse TimeStep time: " '{ print $2 }' | tail -$nsteps)

  if [ -z "$timesteps" ]; then
      echo "-1"
      return
  fi

  # Sort them numerically.

  timesteps=$(echo $timesteps | tr " " "\n" | sort -n)
  last_timestep=$(echo $timesteps | tr " " "\n" | tail -1)

  # Calculate the minimum and maximum.

  maximum_timestep=$(maximum "$timesteps")
  minimum_timestep=$(minimum "$timesteps")

  # Fit a line that matches these extrema.

  pos1=0
  pos2=0

  for timestep in $timesteps
  do
      pos1=$((pos1+1))
      if [ "$minimum_timestep" == "$timestep" ]; then
	  break
      fi
  done

  for timestep in $timesteps
  do
      pos2=$((pos2+1))
      if [ "$maximum_timestep" == "$timestep" ]; then
	  break
      fi
  done

  # Get a "slope". We want a conservative choice here.

  if [ "$maximum_timestep" == "$minimum_timestep" ]; then
      slope=0.0
  else
      slope=$(echo "($maximum_timestep - $minimum_timestep) / ($pos2 - $pos1)" | bc -l)
  fi

  # Number of timesteps to predict forward by (one, in this case).

  dx=1

  next_timestep=$(echo "$last_timestep + $slope * $dx" | bc -l)

  echo $next_timestep

}



# Obtain the length of walltime remaining on the current job.

function get_remaining_walltime {

  if [ ! -z $1 ]; then
      job_number=$1
  else
      echo "Job number not given to get_remaining_walltime" >&2
      return
  fi

  total_time=''

  if [ $batch_system == "PBS" ]; then

      # For PBS we can get the remaining time by doing 
      # showq and then grepping for the line that contains
      # the relevant job number.

      total_time=$(showq -u $USER | grep $job_number | awk '{ print $5 }')

      if [ -z $total_time ]; then
	  echo "Unable to get total_time from job submission system" >&2
      fi

      if [ ! -z $total_time ]; then
	  total_time=$(hours_to_seconds $total_time)
      fi

  elif [ $batch_system == "COBALT" ]; then

      export QSTAT_HEADER=JobId:User:WallTime:RunTime:Nodes:Mode:State:Queue

      # For Cobalt we need to subtract the run time from the total allotted time.

      total_time=$(qstat -u $USER | grep $job_number | awk '{ print $3 }')
      total_time=$(hours_to_seconds $total_time)

      time_used=$(qstat -u $USER | grep $job_number | awk '{ print $4 }')
      time_used=$(hours_to_seconds $time_used)

      if [ ! -z $total_time ] && [ ! -z $time_used ]; then
	  total_time=$(echo "$total_time - $time_used" | bc -l)
      fi

  else

      echo "Unknown job submission system in get_remaining_walltime" >&2

  fi

  echo $total_time

}



# Return a string that is used for restarting from the latest checkpoint.
# Optionally you can hand this a directory, otherwise it will default 
# to whatever get_last_checkpoint determines.

function get_restart_string {

  if [ -z $1 ]; then
      dir="./"
  else
      dir=$1
  fi

  checkpoint=$(get_last_checkpoint $dir)

  # restartString will be empty if no chk files are found -- i.e. this is a new run.

  if [ ! -n "$checkpoint" ]; then
      restartString=""
  else
      restartString="amr.restart="$checkpoint
  fi

  echo $restartString

}



# Return 1 if a job is running in the directory given in the first argument.

function is_job_running {

  if [ -z $1 ]; then
      return 0;
  else
      dir=$1
  fi

  job_status=0

  # Check if a file with the appropriate run extension exists in that directory.
  # If not, check if there is an active job yet to be completed or started 
  # in the directory using the results of the queue information (showq, qstat, etc.).

  if [ -e "$dir/*$run_ext" ]; then

      job_status=1

  elif [ -e $dir/jobs_submitted.txt ] && [ ! -z $num_jobs ]; then

      num_jobs_in_dir=$(cat $dir/jobs_submitted.txt | wc -l)
      jobs_in_directory=$(cat $dir/jobs_submitted.txt | awk '{print $1}')

      for job1 in ${job_arr[@]}
      do
	  for job2 in $jobs_in_directory
	  do
	      if [ $job1 == $job2 ]; then
		  job_status=1
	      fi
          done
      done

  fi

  echo $job_status

}



# Determine if directory dir has reached the desired stopping time or max_step.

function is_dir_done {

  if [ -z $dir ]; then
      directory="./"
  else
      directory=$dir
  fi

  # If the directory already exists, check to see if we've reached the desired stopping point.
  # There are two places we can look: in the last checkpoint file, or in the last stdout file. 

  checkpoint=$(get_last_checkpoint $directory)
  last_output=$(get_last_output $directory)

  # Get the desired stopping time and max step from the inputs file in the directory.
  # Alternatively, we may have set this from the calling script, so prefer that.

  if [ -z $stop_time ] && [ -e $directory/$inputs ]; then
      stop_time=$(get_inputs_var "stop_time" $directory)
  fi

  if [ -z $max_step ] && [ -e $directory/$inputs ]; then
      max_step=$(get_inputs_var "max_step" $directory)
  fi

  # Assume we're not done, by default.

  time_flag=""
  step_flag=""

  done_status=0

  if [ -e "$directory/jobIsDone" ]; then

      # The user has explicitly signalled that the simulation is complete; we can stop here.

      done_status=1

  elif [ -e "$directory/$checkpoint/jobIsDone" ]; then

      # The problem has explicitly signalled that the simulation is complete; we can stop here.

      done_status=1

  elif [ -e "$directory/$checkpoint/jobIsNotDone" ]; then

      # The problem has explicitly signalled that the simulation is NOT complete; again, we can stop here.

      done_status=0

  elif [ -e "$directory/$checkpoint/Header" ]; then

      # Extract the checkpoint time. It is stored in row 3 of the Header file.

      chk_time=$(awk 'NR==3' $directory/$checkpoint/Header)

      # Convert to floating point, since it might be in exponential format.

      chk_time=$(printf "%f" $chk_time)

      # Extract the current timestep. We can get it from the 
      # name of the checkpoint file. cut will do the trick;
      # just capture everything after the 'k' of 'chk'.

      chk_step=$(echo $checkpoint | cut -d"k" -f2)

      if [ ! -z $stop_time ]; then
	  time_flag=$(echo "$chk_time >= $stop_time" | bc -l)
      fi

      if [ ! -z $max_step ]; then
	  step_flag=$(echo "$chk_step >= $max_step" | bc)
      fi

  elif [ ! -z "$last_output" ] && [ -e "$directory/$last_output" ]; then

      output_time=$(grep "STEP =" $directory/$last_output | tail -1 | awk '{print $6}')
      output_step=$(grep "STEP =" $directory/$last_output | tail -1 | awk '{print $3}')

      # bc can't handle numbers in scientific notation, so use printf to convert it to floating point.

      output_time=$(printf "%f" $output_time)

      time_flag=$(echo "$output_time >= $stop_time" | bc -l)
      step_flag=$(echo "$output_step >= $max_step" | bc)

  fi

  # If we don't have valid variables for checking against the timestep and max_time
  # criteria, we assume that we're not done because we just haven't run the job yet.

  if [ ! -z $time_flag ] && [ ! -z $step_flag ]; then

    # If the variables are valid, check if either one indicates that we are done.

    if [ $time_flag -eq 1 ] || [ $step_flag -eq 1 ]; then
      done_status=1
    fi

  fi

  echo $done_status

}



# Archive the file or directory given in the first argument, 
# to the same path on the archive machine relative to the machine's $workdir.

function archive {

  if [ ! -z $1 ]; then
      if [ -d $1 ]; then
	  echo "Archiving contents of directory "$1"."
      else
	  echo "Archiving location "$1"."
      fi
  else
      echo "No file to archive; exiting."
      return
  fi

  # We may get a directory to archive, so call basename to make sure $file
  # doesn't appear with a trailing slash.

  f=$(basename $1)
  d=$(dirname $1)

  # Get the absolute path to this directory, and then 
  # remove everything from the directory up to the username.
  # The assumption here is that everything after that was 
  # created by the user, and that's the directory structure we want
  # to preserve when moving things over to the storage system.

  cd $d 
  abs_path=$(pwd)
  cd - > /dev/null

  local_path=$abs_path
  storage_path=$abs_path

  storage_path=${storage_path#*$USER/}
  storage_path=${storage_path#*$USER}

  # Archive based on the method chosen for this machine.

  if   [ $archive_method == "htar" ]; then

      # htar will give us the path to the file in the 
      # tar file if we do it from outside the local
      # directory, so let's jump in first and avoid that.

      cd $d

      src=$f
      dst=$storage_path/$f.tar

      $HTAR $dst $src

      cd - > /dev/null

  elif [ $archive_method == "globus" ]; then

      src=$globus_src_endpoint/$local_path/$f
      dst=$globus_dst_endpoint/$storage_path/$f

      if [ -d $d/$f ]; then
          # If we're transferring a directory, Globus needs to explicitly know
          # that it is recursive, and needs to have trailing slashes.
          $globus_archive -- $src/ $dst/ -r
      else
          # We're copying a normal file.
	  $globus_archive -- $src $dst
      fi

  fi

}



# Archive all the output files in the directory given in the first argument.
# The directory must be an absolute path. The strategy will be to determine 
# all files we want to archive, then create a list and pass that list to the
# main archive function.

function archive_all {

  if [ ! -z $1 ]; then
      directory=$1
  else
      directory="./"
  fi

  if [ ! -d $directory/output/ ]; then
      mkdir $directory/output/
  fi

  archivelist=""

  # Archive the plotfiles and checkpoint files.
  # Make sure that they have been completed by checking if
  # the Header file exists, which is the last thing created.

  pltlist=$(find $directory -maxdepth 1 -type d -name "*plt*" | sort)
  chklist=$(find $directory -maxdepth 1 -type d -name "*chk*" | sort)

  # Move all completed plotfiles and checkpoints to the output
  # directory, and add them to the list of things to archive.
  # It is possible that a plotfile or checkpoint of the same name
  # will be duplicated if we started from an earlier checkpoint;
  # in this case, delete the old one and replace it with the new.

  for file in $pltlist
  do
      if [ -e $file/Header ]; then
	  if [ -e output/$file ]; then
	      rm -rf output/$file
	  fi
	  mv $file $directory/output/
	  f=$(basename $file)
	  archivelist=$archivelist" "$f
      fi
  done

  lastCheckpoint=$(get_last_checkpoint $directory)

  for file in $chklist
  do
      # We need to be careful with checkpoints because we do not want 
      # to move a checkpoint that is being read to restart the next 
      # run. So we will skip the last checkpoint under the assumption 
      # that this is the only one needed to perform a restart.

      f=$(basename $file)

      if [ "$f" != "$lastCheckpoint" ]; then
	  if [ -e $file/Header ]; then
	      if [ -e output/$file ]; then
		  rm -rf output/$file
	      fi
	      mv $file $directory/output/
	      archivelist=$archivelist" "$f
	  fi
      fi
  done

  diaglist=$(find $directory -maxdepth 1 -name "*diag*.out")

  # For the diagnostic files, we just want to make a copy and move it to the 
  # output directory; we can't move it, since the same file needs to be there
  # for the duration of the simulation if we want a continuous record. But 
  # we want to avoid archiving the files again if the run has already been
  # completed, so we check the timestamps and only move the file to the output
  # directory if the archived version is older than the main version.

  for file in $diaglist
  do
      f=$(basename $file)
      if [ -e $directory/output/$f ]; then
	  if [ $directory/output/$f -nt $file ]; then
	      continue
	  fi
      fi
      cp $file $directory/output/

      archivelist=$archivelist" "$f
  done

  # Same thing for the runtime stdout files.

  outlist=$(find $directory -maxdepth 1 -name "*$job_name*")

  for file in $outlist
  do
      f=$(basename $file)
      if [ -e $directory/output/$f ]; then
	  if [ $directory/output/$f -nt $file ]; then
	      continue
	  fi
      fi
      cp $file $directory/output/
      archivelist=$archivelist" "$f
  done

  # Same strategy for the inputs and probin files.

  inputs_list=$(find $directory -maxdepth 1 -name "$inputs")

  for file in $inputs_list
  do
      f=$(basename $file)
      if [ -e $directory/output/$f ]; then
	  if [ $directory/output/$f -nt $directory/$f ]; then	  
	      continue
	  fi
      fi
      cp $file $directory/output/
      archivelist=$archivelist" "$f
  done

  probin_list=$(find $directory -maxdepth 1 -name "*probin*")

  for file in $probin_list
  do
      f=$(basename $file)
      if [ -e $directory/output/$f ]; then
	  if [ $directory/output/$f -nt $directory/$f ]; then
	      continue
	  fi
     fi
     cp $file $directory/output/
     archivelist=$archivelist" "$f
  done

  # If there is nothing to archive,
  # then assume we have completed the run and exit.

  if [[ -z $archivelist ]]; then
      return
  fi

  # Now we'll do the archiving for all files in $archivelist.
  # Determine the archiving method based on machine.

  if [ $do_storage -eq 1 ]; then

    if   [ $MACHINE == "TITAN"       ]; then

	# For Titan, just loop over every file we're archiving and htar it.

	for file in $archivelist
	do
	    archive $directory/output/$file
	done

    elif [ $MACHINE == "BLUE_WATERS" ]; then

	# For Blue Waters, we're using Globus Online, which has a cap on the number 
	# of simultaneous transfers a user can have. Therefore our strategy is
	# to sync the entire output directory of this location rather than 
	# transferring the files independently.

	archive $directory/output/

    fi

  fi

}



function get_safety_factor {

  if [ -z $1 ]; then
      echo "No walltime total passed to get_safety_factor; exiting."
      return
  else
      tot_time=$1
  fi

  twoHours=$(hours_to_seconds 2:00:00)

  # For small enough jobs we need to introduce a little extra buffer 
  # if we archiving during this run.

  if [ $tot_time -le $twoHours ] && [ -z $archive_queue ]; then
      safety_factor=0.2
  else
      safety_factor=0.1
  fi

  echo $safety_factor

}



function check_to_stop {

  # Get the job number if we are not on a system that
  # provides the job ID through a runtime variable.

  if [ -z "$job_number" ]; then
      job_number=$(get_last_submitted_job)

      # Clean out any extraneous information.

      job_number=${job_number%%.*}
  fi

  # Get the current UNIX time in seconds.

  start_wall_time=$(date +%s)
  curr_wall_time=$start_wall_time

  # Determine how much time the job has, in seconds.
  
  total_time=$(get_remaining_walltime $job_number)

  # Account for the possibility that we don't yet have
  # this information because the job is just starting;
  # cycle until we do.

  while [ -z $total_time ]
  do
      sleep 1
      total_time=$(get_remaining_walltime $job_number)
  done

  # Now we'll plan to stop when (1 - safety_factor) of the time has been used up.

  safety_factor=$(get_safety_factor $total_time)
  time_remaining=$(echo "(1.0 - $safety_factor) * $total_time" | bc -l)

  end_wall_time=$(echo "$start_wall_time + $time_remaining" | bc)

  # Round to nearest integer.

  end_wall_time=$(printf "%.0f" $end_wall_time)

  # We'll subdivide the remaining interval into a given number of chunks,
  # and periodically wake up to check if we're past the time limit. This
  # is intended to deal with potential issues where the function doesn't
  # wake up on time, which I have seen in the past when a system is overloaded.

  numSleepIntervals=1000

  intervalsElapsed=0
  sleepInterval=$(echo "$time_remaining / $numSleepIntervals" | bc -l)

  numCheckpointIntervals=10
  checkpointInterval=$(echo "$time_remaining / $numCheckpointIntervals" | bc)

  nextCheckpointTime=$(echo "$curr_wall_time + $checkpointInterval" | bc)

  while [ $intervalsElapsed -lt $numSleepIntervals ] && [ $curr_wall_time -lt $end_wall_time ]
  do

      if [ ! -z $pid ]; then
	  if [ $(ps -p $pid | wc -l) -eq 1 ]; then
	      break
	  fi
      fi

      sleep $sleepInterval

      curr_wall_time=$(date +%s)

      intervalsElapsed=$(echo "$intervalsElapsed + 1" | bc)

      # We also use a safety check based on the length of the recent timesteps.
      # If we predict that we're close to the end, break out of here and stop the run.

      num_safety_timesteps=5
      file_to_check=$(get_last_output .)
      predicted_timestep=$(predict_next_timestep $file_to_check 10)

      # Make sure that we round up to the nearest second.
      # This relies on the trick that bc will truncate an integer
      # if you don't use the floating point extension with -l.

      if [ "$predicted_timestep" != "-1" ]; then

	  time_to_stop=$(echo "$curr_wall_time + ($num_safety_timesteps * $predicted_timestep + 1) / 1" | bc)

	  if [ "$time_to_stop" -ge $end_wall_time ]; then
	      echo ""
	      echo "Stopping run since we are possibly within "$num_safety_timesteps" timesteps of the end."
	      echo ""
	      break
	  fi

      fi

      # Periodically dump checkpoints as a safeguard against system crashes.
      # Obviously for this to work properly, we need sleepInterval << checkpointInterval.

      if [ $curr_wall_time -gt $nextCheckpointTime ]; then
	  touch "dump_and_continue"
	  nextCheckpointTime=$(echo "$curr_wall_time + $checkpointInterval" | bc)
      fi

      # Delete checkpoints as we go so that we don't fill up the disk with temporaries.
      # We will only delete the second-to-last, not the last, for safety.
      # Only do this if the checkpoint is one of our temporary checkpoints to
      # guard against a crash, not if it's a regularly spaced checkpoint printed
      # at the request of BoxLib. We can test this by checking if the timestep
      # number or simulation time is a multiple of amr.check_int or amr.check_per.

      checkpoint=$(get_last_checkpoint . 2)

      if [ ! -z $checkpoint ] && [ -d $checkpoint ]; then

	  # Extract the checkpoint time. It is stored in row 3 of the Header file.

	  chk_time=$(awk 'NR==3' $checkpoint/Header)

	  # Convert to floating point, since it might be in exponential format.

	  chk_time=$(printf "%f" $chk_time)

	  # Extract the current timestep. We can get it from the
	  # name of the checkpoint file. cut will do the trick;
	  # just capture everything after the 'k' of 'chk'.

	  chk_step=$(echo $checkpoint | cut -d"k" -f2)

	  # Now get the I/O parameters from the inputs file.

	  chk_int=$(get_inputs_var "amr_check_int" .)
	  chk_per=$(get_inputs_var "amr_check_per" .)

	  if [ $chk_int -le 0 ] && [ $(echo "$chk_per <= 0" | bc -l) -eq 1 ]; then

	      rm -rf $checkpoint

	  elif [ $chk_int -gt 0 ]; then

	      if [ $(( $chk_step % $chk_int )) -eq 0 ]; then

		  rm -rf $checkpoint

	      fi

	  else # chk_per > 0

	      if [ $(echo "$chk_per <= 0" | bc -l) -eq 1 ]; then

		  rm -rf $checkpoint

	      else

		  # Get the coarse timestep. There's a row in the Header that lists
		  # the timestep on each level but to make things tricky, the row is
		  # dependent on the number of levels. The bottom line is that there
		  # are five fixed lines of output, and then nlevs+1 lines after that
		  # describing the domain geometry (where nlevs == amr.max_level + 1),
		  # and then one more line with the refinement ratios, and then finally
		  # the line with the timesteps (where the level 0 timestep is first).
		  # See Amr::checkPoint if you don't believe me and/or this changes
		  # in the future. So the only variable we need to determine where to
		  # look for dt is the maximum number of levels, which conveniently
		  # is one of those fixed lines (the fourth, in particular).

		  nlevs=$((1 + $(awk 'NR==4' $checkpoint/Header)))

		  line=$((5 + $nlevs + 1 + 2))

		  coarseDt=$(awk -v line=$line 'FNR==line {print $1}' $checkpoint/Header)

		  # Convert to floating point, since it might be in exponential format.

		  coarseDt=$(printf "%f" $coarseDt)

		  # Use the same logic as in the Amr class for determining
		  # if we just wrote a checkpoint using amr.check_per.

		  num_per_old=$(printf "%.0f" $(echo "($chk_time - $coarseDt) / $chk_per" | bc -l))
		  num_per_new=$(printf "%.0f" $(echo "($chk_time            ) / $chk_per" | bc -l))

		  if [ "$num_per_old" == "$num_per_new" ]; then

		      rm -rf $checkpoint

		  fi

	      fi

	  fi

      fi

  done

  # BoxLib's framework requires a particular file name to exist in the local directory, 
  # to trigger a checkpoint and quit.

  touch "dump_and_stop"

  # Now handle the case where the job is running on one timestep so long that dump_and_stop
  # was not sufficient. What we will do is send a hard kill to the job if we get halfway to the
  # distance between the safety factor and the end.

  numSleepIntervals=1000

  end_wall_time=$(echo "$start_wall_time + (1.0 - $safety_factor / 2.0) * $total_time" | bc -l)
  end_wall_time=$(printf "%.0f" $end_wall_time)
  time_remaining=$(echo "$end_wall_time - $curr_wall_time" | bc -l)
  sleepInterval=$(echo "$time_remaining / $numSleepIntervals" | bc -l)

  while [ $curr_wall_time -lt $end_wall_time ]
  do

      if [ ! -z $pid ]; then
	  if [ $(ps -p $pid | wc -l) -eq 1 ]; then
	      break
	  fi
      fi

      sleep $sleepInterval

      curr_wall_time=$(date +%s)

  done

  # Check to make sure we are done, and if not, re-submit the job.

  if [ -z "$no_continue" ]; then

      if [ ! -e no_submit ]; then
	  dir_done_status=$(is_dir_done)
	  if [ "$dir_done_status" != "1" ]; then
	     submit_job
	  fi
      else
	  rm -f no_submit
      fi

  fi

  # Run the archive script at the end of the simulation.

  if [ ! -z $archive_queue ]; then

      # We want to submit the archive job from the output/ directory.
      # The point is that our job software gets confused if there are
      # multiple jobs running in a directory at once, so we arbitrarily
      # submit this job from a different location so that the other functions
      # here don't see multiple outputs from separate runs.

      if [ ! -d output/ ]; then
	  mkdir output/
      fi

      cd output/
      archive_job_number=`$exec $archive_script`
      cd - > /dev/null
      echo "Submitted an archive job with job number $archive_job_number."

  else

      if [ $do_storage -ne 0 ]; then
	  archive_all
      fi

  fi

  rm -f dump_and_stop
  rm -f dump_and_continue

  if [ ! -z "$job_number" ] && [ ! -z "$cancel_job" ]; then
      if [ ! -z "$pid" ] && [ $(ps -p $pid | wc -l) -ne 1 ]; then
	  echo "Implementing a hard termination of the job since we are too close to the time limit."
	  $cancel_job $job_number
      fi
  fi

}



# Copies all relevant files needed for a CASTRO run into the target directory,
# and updates the inputs and probin according to any shell variables we set.

function copy_files {

  if [ -z $dir ]; then
      echo "No directory passed to copy_files; exiting."
  fi

  if [ ! -d "$dir/output/" ]; then
      mkdir $dir/output/
  fi

  if [ ! -e $dir/$CASTRO ] && [ -z "$inputs_only" ]; then
      if [ ! -z "$local_compile" ] && [ "$local_compile" -eq "1" ]; then
	  if [ -e $dir/$compile_dir/$CASTRO ]; then
	      cp $dir/$compile_dir/$CASTRO $dir
	  fi
      elif [ -e $compile_dir/$CASTRO ]; then
	  cp $compile_dir/$CASTRO $dir
      fi
  fi

  if [ ! -e "$dir/helm_table.dat" ] && [ -z "$inputs_only" ]; then
      if [ -e "$compile_dir/helm_table.dat" ]; then
	  cp $compile_dir/helm_table.dat $dir
      fi
  fi

  new_inputs="F"

  if [ ! -e "$dir/$inputs" ]; then

      new_inputs="T"

      if [ -e "$compile_dir/$inputs" ]; then
          cp $compile_dir/$inputs $dir/$inputs
      else
	  cp $problem_dir/$inputs $dir/$inputs
      fi
  fi

  new_probin="F"

  if [ ! -e "$dir/$probin" ]; then

      new_probin="T"

      if [ -e "$compile_dir/$probin" ]; then
	  cp $compile_dir/$probin $dir/$probin
      else
	  cp $problem_dir/$probin $dir/$probin
      fi
  fi

  # Copy over all the helper scripts, so that these are 
  # fixed in time for this run and don't change if we update the repository.

  if [ ! -e "$dir/job_scripts/run_utils.sh" ] && [ -z "$inputs_only" ]; then
      mkdir -p "$dir/job_scripts"
      cp -r $WDMERGER_HOME/job_scripts/*.sh $dir/job_scripts/
  fi

  if [ -z "$inputs_only" ]; then
      touch "$dir/jobs_submitted.txt"
  fi

  # Now determine all the variables that have been added
  # since we started; then search for them in the inputs
  # file and do a replace as needed. This relies on the 
  # comm function, which when run with the -3 option 
  # returns only the strings that aren't common to 
  # both of two files. To get our variable list to 
  # play nice with it, we use tr to replace spaces
  # with newlines, so that comm thinks it's being 
  # handled a file in the same format as if you did ls.

  shell_list_new=$(compgen -v)

  input_vars=$(comm -3 <( echo $shell_list | tr " " "\n" | sort) <( echo $shell_list_new | tr " " "\n" | sort))

  # Loop through all new variables and call both replace_inputs_var and 
  # replace_probin_var. These will only take action if the variable exists
  # in the respective files, and there should not be any common variables,
  # so there is no harm in the redundancy.

  for var in $input_vars
  do

      if [ $new_inputs == "T" ] || (([ $var == "stop_time" ] || [ $var == "max_step" ]) && [ ! -z "$update_stopping_criteria" ]); then
	  replace_inputs_var $var
      fi

      if [ $new_probin == "T" ]; then
	  replace_probin_var $var
      fi

  done

}



# Submits a job from the job script, then captures the 
# job number output. This function is designed to be run 
# while we live in the job directory itself.

function submit_job {

  # Get the current date for printing to file so we know 
  # when we submitted. Since we only one want to add one column,
  # we'll use the +%s option, which is seconds since January 1, 1970.

  current_date=$(date +%s)

  # Sometimes the code crashes and we get into an endless cycle of 
  # resubmitting the job and then crashing again soon after,
  # which is liable to make system administrators mad at us.
  # Let's protect against this by putting in a safeguard.
  # If the job stops within the first 25% of its requested runtime,
  # it is a safe bet that we are crashing and we don't want to
  # submit a new job.

  old_date=$(tail -1 jobs_submitted.txt | awk '{print $2}')
  old_walltime=$(tail -1 jobs_submitted.txt | awk '{print $3}')
  old_nprocs=$(tail -1 jobs_submitted.txt | awk '{print $4}')

  if [ ! -z $old_date ] && [ ! -z $old_walltime ]; then

      date_diff=$(( $current_date - $old_date ))

      safety_factor=$(get_safety_factor $old_walltime)
      submit_flag=$(echo "$date_diff > 0.25 * $old_walltime" | bc -l)

      if [ $submit_flag -eq 0 ]; then
	  echo "Refusing to submit job because the last job ended too soon."
	  return 1
      fi

  fi

  # We can also test directly on whether the code crashed by looking
  # for a BoxLib Backtrace file.

  if ls Backtrace* 1> /dev/null 2>&1; then
      echo "Refusing to submit job because a Backtrace file was found."
      return 1
  fi

  # Determine the requested walltime, in seconds.

  if [ -z $walltime ]; then
      if [ ! -z $old_walltime ]; then
	  walltime=$old_walltime
      else
	  echo "Don't know what the walltime request is in submit_job; aborting."
	  return 1
      fi
  fi

  walltime_in_seconds=$(hours_to_seconds $walltime)
  walltime_in_minutes=$(hours_to_minutes $walltime)

  # Determine the number of nodes; some job submission systems
  # (e.g. Cobalt) need to know this at submission.

  if [ -z $nprocs ]; then
      if [ ! -z $old_nprocs ]; then
	  nprocs=$old_nprocs
      else
	  echo "Don't know how many processors this job needs; aborting."
	  return 1
      fi
  fi

  nodes=$(compute_num_nodes)

  # If we made it to this point, now actually submit the job.

  if [ $batch_system == "PBS" ]; then
      submitted_job_number=`$exec $job_script`
  elif [ $batch_system == "COBALT" ]; then
      submitted_job_number=`$exec -A $allocation -t $walltime_in_minutes -n $nodes --mode script run_script`
  fi

  # Some systems like Blue Waters include the system name
  # at the end of the number, so remove any appended text.

  submitted_job_number=${submitted_job_number%%.*}

  if [ ! -z "$submitted_job_number" ]; then
      echo "$submitted_job_number $current_date $walltime_in_seconds $nprocs" >> jobs_submitted.txt
  fi

}



# Get the last job number submitted by examining the jobs_submitted file.

function get_last_submitted_job {

  if [ -e jobs_submitted.txt ]; then

      job_number=$(tail -1 jobs_submitted.txt | awk '{print $1}')

  else

      echo "No jobs_submitted.txt file to obtain last job number from" >&2

      job_number=-1

  fi

  echo $job_number

}



# Generate a run script in the given directory.

function create_job_script {

  if [ ! -z $1 ]; then
      dir=$1
  else
      echo "No directory given to create_job_script; exiting."
      return
  fi

  if [ ! -z $2 ]; then
      nprocs=$2
  else
      echo "Number of processors not given to create_job_script; exiting."
      return
  fi

  if [ ! -z $3 ]; then
      walltime=$3
  else
      echo "Walltime not given to create_job_script; exiting."
      return
  fi

  nodes=$(compute_num_nodes)

  # Number of threads for OpenMP. This will be equal to 
  # what makes the most sense for the machine architecture 
  # by default. For example, the Titan XK7 and Blue Waters XE6
  # Cray nodes are composed of Interlagos boards which are composed 
  # of two NUMA nodes (each NUMA node has 8 integer cores and 4 
  # floating point cores). If the user doesn't set it,
  # we'll update the default with our experience from running on 
  # these machines with tiling. When the grids are small enough,
  # there isn't enough work to justify the OpenMP overhead. So 
  # we'll use two OpenMP threads for small problems and four
  # threads for bigger problems.

  if [ -z $OMP_NUM_THREADS ]; then

      if [ -z $threads_per_task ]; then
	  threads_per_task=1
      fi

      OMP_NUM_THREADS=$threads_per_task

      # First, get the maximum grid size. If this has been
      # set by the including script, we use that; otherwise, 
      # we read in the value from the main inputs file.

      if [ -z "$amr_max_grid_size" ]; then
	  amr_max_grid_size=$(get_inputs_var "amr_max_grid_size" $dir)
      fi

      max_level_grid_size=0

      for grid_size in $amr_max_grid_size
      do
	  if [ $grid_size -gt $max_level_grid_size ]; then
	      max_level_grid_size=$grid_size
	  fi
      done

      if [ $MACHINE == "TITAN" ]; then
	  if [ $max_level_grid_size -lt "64" ]; then
	      OMP_NUM_THREADS=2
	  elif [ $max_level_grid_size -lt "128" ]; then
	      OMP_NUM_THREADS=4
	  fi
      elif [ $MACHINE == "BLUE_WATERS" ]; then
	  if [ $max_level_grid_size -lt "64" ]; then
	      OMP_NUM_THREADS=2
	  elif [ $max_level_grid_size -lt "128" ]; then
	      OMP_NUM_THREADS=4
	  fi
      fi

      # Also, we want to make sure that OMP_NUM_THREADS is equal to one
      # if we didn't compile with OpenMP.

      do_omp=$(get_make_var USE_OMP)
      if [ $do_omp == "FALSE" ]; then
	  OMP_NUM_THREADS=1
      fi

  fi

  # If the number of processors is less than the number of processors per node,
  # there are scaling tests where this is necessary; we'll assume the user understands
  # what they are doing and set it up accordingly.

  old_ppn=$ppn

  if [ $nodes -eq 0 ]; then
      nodes="1"
      ppn=$nprocs
  fi

  num_mpi_tasks=$(echo "$nprocs / $OMP_NUM_THREADS" | bc)
  tasks_per_node=$(echo "$ppn / $OMP_NUM_THREADS" | bc)

  # Create the job script and make it executable.

  touch $dir/$job_script
  chmod u+x $dir/$job_script

  if [ $batch_system == "PBS" ]; then

      echo "#!/bin/bash" > $dir/$job_script

      # Select the project allocation we're charging this job to
      if [ ! -z $allocation ]; then
	  echo "#PBS -A $allocation" >> $dir/$job_script
      fi

      # Set the name of the job
      echo "#PBS -N $job_name" >> $dir/$job_script

      # Combine standard error into the standard out file
      echo "#PBS -j oe" >> $dir/$job_script

      # Amount of wall time for the simulation
      echo "#PBS -l walltime=$walltime" >> $dir/$job_script

      # Number of nodes, the number of MPI tasks per node, and the node type to use
      if [ $MACHINE == "BLUE_WATERS" ]; then
	  echo "#PBS -l nodes=$nodes:ppn=$ppn:$node_type" >> $dir/$job_script
      elif [ $MACHINE == "HOPPER" ]; then
	  echo "#PBS -l mppwidth=$nprocs" >> $dir/$job_script
      elif [ $MACHINE == "LIRED" ]; then
	  echo "#PBS -l nodes=$nodes:ppn=$ppn" >> $dir/$job_script
      else
	  echo "#PBS -l nodes=$nodes" >> $dir/$job_script
      fi

      # Queue to submit to. This is required for some systems.
      if [ ! -z $queue ]; then
	  echo "#PBS -q $queue" >> $dir/$job_script
      fi

      echo "" >> $dir/$job_script

      # Set name of problem directory, if applicable.

      if [ ! -z $problem_dir ]; then
	  echo "problem_dir=$problem_dir" >> $dir/$job_script
	  echo "" >> $dir/$job_script
      fi

      # Set the name of the inputs and probin files, in case they are 
      # unique for this problem.

      if [ $inputs != "inputs" ]; then
	  echo "inputs=$inputs" >> $dir/$job_script
	  echo "" >> $dir/$job_script
      fi

      if [ $probin != "probin" ]; then
	  echo "probin=$probin" >> $dir/$job_script
	  echo "" >> $dir/$job_script
      fi

      # Indicate if we don't want to continue the run.

      if [ ! -z $no_continue ]; then
	  echo "no_continue=$no_continue" >> $dir/$job_script
	  echo "" >> $dir/$job_script
      fi

      # Indicate if we don't want to archive the data.

      if [ ! -z $do_storage ]; then
	  echo "do_storage=$do_storage" >> $dir/$job_script
	  echo "" >> $dir/$job_script
      fi

      # Store the job number.

      echo "job_number=\$PBS_JOBID" >> $dir/$job_script
      echo "job_number=\${job_number%%.*}" >> $dir/$job_script
      echo "" >> $dir/$job_script

      # We assume that the directory we submitted from is eligible to 
      # work in, so cd to that directory.

      echo "cd \$PBS_O_WORKDIR" >> $dir/$job_script
      echo "" >> $dir/$job_script

      # Load up our helper functions.

      echo "source job_scripts/run_utils.sh" >> $dir/$job_script
      echo "" >> $dir/$job_script

      # Number of OpenMP threads

      echo "export OMP_NUM_THREADS=$OMP_NUM_THREADS" >> $dir/$job_script

      # Amount of memory allocated to each OpenMP thread.

      echo "export OMP_STACKSIZE=64M" >> $dir/$job_script
      echo "" >> $dir/$job_script

      # Set the aprun options.

      if [ $launcher == "aprun" ]; then
	  launcher_opts="-n $num_mpi_tasks -N $tasks_per_node -d $OMP_NUM_THREADS"
	  redirect=""
      elif [ $launcher == "mpirun" ]; then
	  launcher_opts="-np $num_mpi_tasks --map-by ppr:$threads_per_task"
	  redirect="> $job_name.OU"
      fi

      # Main job execution.

      echo "$launcher $launcher_opts $CASTRO $inputs \$(get_restart_string) $redirect &" >> $dir/$job_script
      echo "" >> $dir/$job_script

      echo "pid=\$!" >> $dir/$job_script
      echo "" >> $dir/$job_script

      # Call the function that determines when we're going to stop the run.
      # It should run in the background, to allow the main job to execute.

      echo "check_to_stop &" >> $dir/$job_script
      echo "" >> $dir/$job_script

      echo "wait" >> $dir/$job_script
      echo "" >> $dir/$job_script

      # With mpirun we redirect the output to a file; let's move that to a file 
      # named by the job number so that we retain it later.

      if [ $launcher == "mpirun" ]; then

	echo "mv $job_name.OU \$job_number.out" >> $dir/$job_script
	echo "" >> $dir/$job_script

      fi



      # Now write the archive script, if this is a system that 
      # has a data transfer queue.

      if [ ! -z $archive_queue ]; then

	  echo "#!/bin/bash" > $dir/output/$archive_script

          # Select the project allocation we're charging this job to
	  if [ ! -z $allocation ]; then
	      echo "#PBS -A $allocation" >> $dir/output/$archive_script
	  fi

           # Set the name of the job
	  echo "#PBS -N archive" >> $dir/output/$archive_script

           # Combine standard error into the standard out file
	  echo "#PBS -j oe" >> $dir/output/$archive_script

          # Amount of wall time for the simulation
	  echo "#PBS -l walltime=$archive_wclimit" >> $dir/output/$archive_script

          # Number of nodes, the number of MPI tasks per node, and the node type to use
	  echo "#PBS -l nodes=1" >> $dir/output/$archive_script

          # Queue to submit to.
	  echo "#PBS -q $archive_queue" >> $dir/output/$archive_script

	  echo "" >> $dir/output/$archive_script

	  # Set the location of some variables.

	  echo "WDMERGER_HOME=$WDMERGER_HOME" >> $dir/output/$archive_script
	  echo "MACHINE=$MACHINE" >> $dir/output/$archive_script
	  echo "" >> $dir/output/$archive_script

          # We assume that the directory we submitted from is eligible to 
          # work in, so cd to that directory. We submit from the output
	  # directory, so move up one level before actually executing the
	  # scripts.

	  echo "cd \$PBS_O_WORKDIR" >> $dir/output/$archive_script
	  echo "cd ../" >> $dir/output/$archive_script
	  echo "" >> $dir/output/$archive_script

          # Load up our helper functions.

	  echo "source job_scripts/run_utils.sh" >> $dir/output/$archive_script
	  echo "" >> $dir/output/$archive_script

          # Main job execution.

	  if [ $do_storage -ne 1 ]; then
	      echo "do_storage=$do_storage" >> $dir/output/$archive_script
	  fi
	  echo "archive_all >> archive_log.out" >> $dir/output/$archive_script
	  echo "" >> $dir/output/$archive_script

      fi

  elif [ $batch_system == "COBALT" ]; then

      echo "#!/bin/bash" > $dir/$job_script

      echo "" >> $dir/$job_script

      # Set name of problem directory, if applicable.

      if [ ! -z $problem_dir ]; then
	  echo "problem_dir=$problem_dir" >> $dir/$job_script
	  echo "" >> $dir/$job_script
      fi

      # Set the name of the inputs and probin files, in case they are 
      # unique for this problem.

      if [ $inputs != "inputs" ]; then
	  echo "inputs=$inputs" >> $dir/$job_script
	  echo "" >> $dir/$job_script
      fi

      if [ $probin != "probin" ]; then
	  echo "probin=$probin" >> $dir/$job_script
	  echo "" >> $dir/$job_script
      fi

      # Indicate if we don't want to continue the run.

      if [ ! -z $no_continue ]; then
	  echo "no_continue=$no_continue" >> $dir/$job_script
	  echo "" >> $dir/$job_script
      fi

      # Indicate if we don't want to archive the data.

      if [ ! -z $do_storage ]; then
	  echo "do_storage=$do_storage" >> $dir/$job_script
	  echo "" >> $dir/$job_script
      fi

      # Load up our helper functions.

      echo "source job_scripts/run_utils.sh" >> $dir/$job_script
      echo "" >> $dir/$job_script

      if [ $launcher == "runjob" ]; then
	  launcher_opts="--np $nprocs -p $tasks_per_node --block=\$COBALT_PARTNAME --verbose=INFO : "
      fi

      redirect="> $job_name.OU"

      # Main job execution.

      echo "$launcher $launcher_opts $CASTRO $inputs \$(get_restart_string) $redirect" >> $dir/$job_script
      echo "" >> $dir/$job_script

      echo "pid=\$!" >> $dir/$job_script
      echo "" >> $dir/$job_script

      # Call the function that determines when we're going to stop the run.
      # It should run in the background, to allow the main job to execute.

      echo "check_to_stop &" >> $dir/$job_script
      echo "" >> $dir/$job_script

      echo "wait" >> $dir/$job_script
      echo "" >> $dir/$job_script

   elif [ $batch_system == "batch" ]; then

      echo "echo \"mpiexec -n $nprocs $CASTRO $inputs > $job_name$run_ext\" | batch" > $dir/$job_script

   fi

   # Restore the number of processors per node in case we changed it.

   ppn=$old_ppn

}



# Delete the last submitted job in the directory.

function cancel {

  if [ -z $dir ]; then
      echo "No directory given to cancel; exiting."
      return
  fi

  if [ -d $dir ]; then

      cd $dir

      job_number=$(get_last_submitted_job)

      if [ $job_number -gt 0 ]; then

	  echo "Cancelling job number $job_number in directory $dir."

	  $cancel_job $job_number

      fi

      cd - > /dev/null

  fi

}



# Cancel the last submitted job, but do it gracefully using
# a checkpoint dump.

function soft_cancel {

  if [ -z $dir ]; then
      echo "No directory given to soft_cancel; exiting."
      return
  fi

  if [ -d $dir ]; then

      cd $dir

      job_number=$(get_last_submitted_job)

      if [ $job_number -gt 0 ]; then

	  echo "Cancelling job number $job_number in directory $dir."

	  touch "no_submit"
	  touch "dump_and_stop"

      fi

      cd - > /dev/null

  fi

}



# Pause the last submitted job in the directory.

function pause {

  if [ -z $dir ]; then
      echo "No directory given to pause; exiting."
      return
  fi

  if [ -d $dir ]; then

      cd $dir

      job_number=$(get_last_submitted_job)

      if [ $job_number -gt 0 ]; then

	  echo "Pausing job number $job_number in directory $dir."

	  $pause_job $job_number

      fi

      cd - > /dev/null

  fi

}



# Resume the last submitted job in the directory.

function resume {

  if [ -z $dir ]; then
      echo "No directory given to resume; exiting."
      return
  fi

  if [ -d $dir ]; then

      cd $dir

      job_number=$(get_last_submitted_job)

      if [ $job_number -gt 0 ]; then

	  echo "Resuming job number $job_number in directory $dir."

	  $resume_job $job_number

      fi

      cd - > /dev/null

  fi

}



# Main submission function. Checks which Linux variant we're on,
# and uses the relevant batch submission script. If you want to
# use a different machine, you'll need to include a run script
# for it in the job_scripts directory.
# The variable dir must exist in the calling script, and is the 
# name of the directory where we want to # submit this job from.
# Optionally, if nprocs is defined in the calling script,
# that is the number of processors to use; otherwise we default to
# using all processors available on one node.
# If walltime isn't defined, we run for one hour.

function run {

  # Require the user to tell us where the source files and inputs are coming from.

  if [ -z $problem_dir ]; then
      echo "No problem directory specified; exiting."
      return
  fi

  if [ -z $dir ]; then
      echo "No directory given to run; exiting."
      return
  fi

  if [ -z $DIM ]; then
      DIM="3"
  fi

  set_up_problem_dir

  if [ -z $nprocs ]; then
      nprocs=$ppn
  fi

  if [ -z $walltime ]; then
      walltime=1:00:00
  fi

  do_job=0

  if [ ! -d $dir ]; then

    if [ -z "$no_submit" ] && [ -z "$inputs_only" ]; then
  	echo "Creating directory "$dir" and submitting job."
    else
  	echo "Creating directory "$dir" without submitting the job."
    fi

    mkdir -p $dir

    do_job=1

  elif [ -z "$inputs_only" ] && [ -z "$no_submit" ]; then

    # First as a sanity check, make sure the desired job isn't already running.

    job_running_status=$(is_job_running $dir)

    if [ $job_running_status -eq 1 ]; then

  	echo "Job currently in process or queued in directory $dir."

    else

      # Remove the dump_and_stop file if it exists.

      rm -f $dir/dump_and_stop

      done_flag=$(is_dir_done)

      if [ $done_flag -eq 0 ]; then

  	  echo "Continuing job in directory $dir."

  	  do_job=1

      else

  	  # If we make it here, then we've already reached either stop_time
  	  # or max_step, so we should conclude that the run is done.

  	  echo "Job has already been completed in directory $dir."

      fi

    fi

  fi

  # Optionally, the user can force a recompile from the run script.

  if [ ! -z "$local_compile" ] && [ -z "$inputs_only" ]; then
      if [ "$local_compile" -eq "1" ]; then
	  compile_in_job_directory $dir
      fi
  fi

  copy_files $dir

  if [ -z "$inputs_only" ]; then
      if [ ! -e "$dir/$job_script" ]; then
	  create_job_script $dir $nprocs $walltime
      fi
  fi

  # If we are continuing or starting a job, change into the run directory, 
  # submit the job, then come back to the main directory.

  if [ $do_job -eq 1 ]; then

    # Sometimes we'll want to set up the run directories but not submit them,
    # e.g. for testing purposes, so if the user creates a no_submit variable,
    # we won't actually submit this job.

    if [ -z "$no_submit" ] && [ -z "$inputs_only" ]; then

  	cd $dir

        # Run the job, and capture the job number for output.

  	submit_job

	return_value=$?

	if [ $return_value -eq 0 ]; then
  	    echo "The job number is $(get_last_submitted_job)."
	fi

  	cd - > /dev/null

    fi

  fi

}



function set_up_problem_dir {

    # Upon initialization, store some variables and create results directory.
    # We only want to initialize these variables if we're currently in a root problem directory.

    if [ -d $compile_dir ]; then

        # Some variables we need for storing job information.

	num_jobs=
	job_arr=()
	state_arr=()
	walltime_arr=()

        # Build the main executable if we haven't yet.

	if [ ls $compile_dir/*"$DIM"d*.ex 1> /dev/null 2>&1 ]; then

	    echo "Detected that the executable doesn't exist yet; building executable now."

	    cd $compile_dir
	    make -j8 $(compile_options) &> compile.out
	    cd - > /dev/null

	    echo "Done building executable."

	fi

        # Fill these arrays.

	get_submitted_jobs

	if [ -e $compile_dir/GNUmakefile ]; then

	    CASTRO=$(get_make_var executable)

	fi

	if [ ! -d $plots_dir ]; then
	    mkdir $plots_dir
	fi


    fi

}



# Compile the executable locally in the job directory.
# This is used if you want to compile with non-standard
# make options, like a different network or EOS.

function compile_in_job_directory {

  if [ -z $1 ]; then

      echo "Error: no directory given to function compile_in_job_directory."

  fi

  if [ ! -d $1/$compile_dir ]; then
      mkdir $1/$compile_dir
  fi

  if [ ! -e $1/$compile_dir/GNUmakefile ]; then
      cp $problem_dir/GNUmakefile $1/$compile_dir
  fi

  if [ ! -e $1/$compile_dir/$CASTRO ]; then

      echo "Locally compiling the executable in directory "$1/$compile_dir"."

      cd $1/$compile_dir
      make -j8 $(compile_options) &> compile.out
      cd - > /dev/null

  fi

}



# The number of nodes is equal to the number of processors divided 
# by the number of processors per node. This will not be an integer 
# if the number of processors requested is not an integer multiple 
# of the number of processors per node. So we want to use a trick 
# that guarantees we round upward so that we have enough nodes:
# if we want the result of (A / B) to be always rounded upward 
# in integer arithmetic, we evaluate (A + B - 1) / B.

function compute_num_nodes {

  nodes=$(echo "($nprocs + $ppn - 1) / $ppn" | bc)

  echo $nodes

}



########################################################################

# Define variables

# Before we get started, save the list of current shell variables.
# We'll use this to be able to sort out only the ones that have been 
# set by this script or the including script.

shell_list=$(compgen -v)

job_name="wdmerger"
job_script="run_script"
archive_script="archive_script"

archive_method="none"

do_storage=1

# Directory to compile the executable in

compile_dir="compile"

# Directory for executing and storing results

results_dir="results"

# Directory for placing plots from analysis routines

plots_dir="plots"

if [ -z $inputs ]; then
    inputs=inputs
fi

if [ -z $probin ]; then
    probin=probin
fi

# Get current machine and set preferences accordingly.
# Note: workdir is the name of the directory you submit 
# jobs from (usually scratch directories).

MACHINE=$(get_machine)

set_machine_params

# Set parameters for our archiving scripts.
if   [ $archive_method == "htar" ]; then
    copies=2
    HTAR="htar -H copies=$copies -Pcvf"
elif [ $archive_method == "globus" ]; then
    globus_username="mkatz"
    globus_hostname="cli.globusonline.org"

    # Give Globus Online a one hour time limit.

    time_limit="1h"

    # If we're transferring a directory, tell Globus to only sync either new files or altered files.

    sync_level=2

    # Main archiving command

    globus_archive="ssh $globus_username@$globus_hostname transfer -d $time_limit -s $sync_level"
fi
