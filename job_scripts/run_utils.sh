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

  if [ ! -z $CASTRO_DIR ]; then
      compile_opts=$compile_opts' CASTRO_DIR='$CASTRO_DIR
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

  if [ $batch_system == "PBS" ]; then

      # Store the result of showq

      job_file=$WDMERGER_HOME/job_scripts/jobs.txt
      showq -u $USER > $job_file.temp

      # Remove the lines that don't store jobs.
      # A useful heuristic is that lines containing our username
      # are going to have jobs printed on that line.

      grep $USER $job_file.temp > $job_file
      rm -f $job_file.temp

      # Store as arrays the status of all jobs.

      num_jobs=$(cat $job_file | wc -l)
      
      for i in $(seq 0 $(($num_jobs-1)))
      do
	  line=$(awk "NR == $i+1" $job_file)
	  job_arr[$i]=$(echo $line | awk '{ print $1 }')
	  state_arr[$i]=$(echo $line | awk '{ print $3 }')
	  walltime_arr[$i]=$(echo $line | awk '{ print $5 }')
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

  # Doing a search this way will treat first any checkpoint files 
  # with seven digits, and then will fall back to ones with six and then five digits.
  # On recent versions of GNU sort, one can simplify this with sort -V.

  checkpoint=$(find $dir -type d -name "*chk???????" | sort | tail -1)

  if [ -z $checkpoint ]; then

      checkpoint=$(find $dir -type d -name "*chk??????" | sort | tail -1)

      if [ -z $checkpoint ]; then

          checkpoint=$(find $dir -type d -name "*chk?????"  | sort | tail -1)

      fi

  fi

  # The Header is the last thing written -- check if it's there, otherwise,
  # fall back to the second-to-last check file written, because it means 
  # the latest checkpoint file is still being written.

  if [ ! -f ${checkpoint}/Header ]; then

      # How many *chk?????? files are there? if only one, then skip,
      # as there are no valid and complete checkpoint files.

      nl=$(find $dir -type d -name "*chk??????" -print | sort | wc -l)

      if [ $nl -gt 1 ]; then
	  checkpoint=$(find $dir -type d -name "*chk??????" -print | sort | tail -2 | head -1)    
      else
	  checkpoint=""
      fi
      
  fi

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



# Obtain the length of walltime remaining on the current job.

function get_remaining_walltime {

  if [ ! -z $1 ]; then
      job_number=$1
  else
      return
  fi

  total_time=0.0

  # Extract the job number from the filename, then
  # use the relevant batch submission system.

  if [ $batch_system == "PBS" ]; then

      # For PBS we can get the remaining time by doing 
      # showq and then grepping for the line that contains
      # the relevant job number.

      total_time=$(showq -u $USER | grep $job_number | awk '{ print $5 }')
	
      total_time=$(hours_to_seconds $total_time)

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
  # in the directory using the results of showq.

  if [ -e $dir/*$run_ext ]; then

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

  time_flag=0
  step_flag=0

  done_status=0

  if [ -e $directory/$checkpoint/Header ]; then

      # Extract the checkpoint time. It is stored in row 3 of the Header file.

      chk_time=$(awk 'NR==3' $directory/$checkpoint/Header)

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

  elif [ ! -z $last_output ] && [ -e $directory/$last_output ]; then

      output_time=$(grep "STEP =" $directory/$last_output | tail -1 | awk '{print $6}')
      output_step=$(grep "STEP =" $directory/$last_output | tail -1 | awk '{print $3}')

      # bc can't handle numbers in scientific notation, so use printf to convert it to floating point.

      output_time=$(printf "%f" $output_time)

      time_flag=$(echo "$output_time >= $stop_time" | bc -l)
      step_flag=$(echo "$output_step >= $max_step" | bc)

  fi

  # If we don't have valid variables for checking against the timestep and max_time
  # criteria, we assume that we're not done because we just haven't run the job yet.

  if [ -z $time_flag ] || [ -z $step_flag ]; then
      done_status=0
  elif [ $time_flag -eq 1 ] || [ $step_flag -eq 1 ]; then
      done_status=1
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
	  echo "Archiving location " $1"."
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
  # remove everything from the directory up to the username. The assumption here is that
  # everything after that was created by the user, and that's the directory structure we want
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

      $HTAR ${d}/${f}.tar $d/$f

  elif [ $archive_method == "globus" ]; then

      src=$globus_src_endpoint/$local_path/$f
      dst=$globus_dst_endpoint/$storage_dir/$f

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

  for file in $chklist
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



function check_to_stop {

  job_number=$(get_last_submitted_job)

  # Determine how much time the job has, in seconds.
    
  total_time=$(get_remaining_walltime $job_number)

  # Now we'll plan to stop when (1 - safety_factor) of the time has been used up.

  time_remaining=$(echo "(1.0 - $safety_factor) * $total_time" | bc -l)

  sleep $time_remaining

  # BoxLib's framework requires a particular file name to exist in the local directory, 
  # to trigger a checkpoint and quit.

  touch "dump_and_stop"

}



# Copies all relevant files needed for a CASTRO run into the target directory,
# and updates the inputs and probin according to any shell variables we set.

function copy_files {

  if [ -z $dir ]; then
      echo "No directory passed to copy_files; exiting."
  fi
  
  if [ ! -e $dir/$CASTRO ]; then
      cp $compile_dir/$CASTRO $dir
  fi

  if [ ! -e "$dir/helm_table.dat" ]; then
      if [ -e "$compile_dir/helm_table.dat" ]; then
	  cp $compile_dir/helm_table.dat $dir
      fi
  fi
  
  if [ ! -e "$dir/$inputs" ]; then
      if [ -e "$compile_dir/$inputs" ]; then
          cp $compile_dir/inputs $dir/$inputs
      else
	  if [ ! -z $problem_dir ]; then
	      cp $problem_dir/$inputs $dir/$inputs
	  else
              cp $WDMERGER_HOME/source/inputs $dir/$inputs
	  fi
      fi
  fi

  if [ ! -e "$dir/$probin" ]; then
      if [ -e "$compile_dir/$probin" ]; then
	  cp $compile_dir/probin $dir/$probin
      else
	  if [ ! -z $problem_dir ]; then
	      cp $problem_dir/$probin $dir/$probin
	  else
	      cp $WDMERGER_HOME/source/probin $dir/$probin
	  fi
      fi
  fi

  # Copy over all the helper scripts, so that these are 
  # fixed in time for this run and don't change if we update the repository.

  if [ ! -e "$dir/job_scripts/run_utils.sh" ]; then
      mkdir -p "$dir/job_scripts"
      cp -r $WDMERGER_HOME/job_scripts/*.sh $dir/job_scripts/
  fi

  touch "$dir/jobs_submitted.txt"

  if [ $DIM -eq "2" ] && [ -z $problem_dir ]; then
      convert_to_2D
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

      replace_inputs_var $var
      replace_probin_var $var

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
  # Normally the job should never end before (1.0 - safety_factor) 
  # of the walltime, so if it has, we know that the job exited 
  # abnormally (or, say, it completed) and so we don't want to
  # submit a new job.

  old_date=$(tail -1 jobs_submitted.txt | awk '{print $2}')
  old_walltime=$(tail -1 jobs_submitted.txt | awk '{print $3}')

  if [ ! -z $old_date ] && [ ! -z $old_walltime ]; then

      date_diff=$(( $current_date - $old_date ))

      submit_flag=$( echo "$date_diff > (1.0 - $safety_factor) * $old_walltime" | bc -l )

      if [ $submit_flag -eq 0 ]; then
	  echo "Refusing to submit job because the last job ended too soon."
	  return
      fi

  fi

  # If we made it to this point, now actually submit the job.

  job_number=`$exec $job_script`

  # Some systems like Blue Waters include the system name
  # at the end of the number, so remove any appended text.

  job_number=${job_number%%.*}

  # Store the requested walltime, in seconds.

  if [ -z $walltime ]; then
      if [ ! -z $old_walltime ]; then
	  walltime=$old_walltime
      else
	  echo "Don't know what the walltime request is in submit_job; aborting."
	  return
      fi
  fi

  walltime_in_seconds=$(hours_to_seconds $walltime)

  echo "$job_number $current_date $walltime_in_seconds" >> jobs_submitted.txt

}



# Get the last job number submitted by examining the jobs_submitted file.

function get_last_submitted_job {

  if [ -e jobs_submitted.txt ]; then

      job_number=$(tail -1 jobs_submitted.txt | awk '{print $1}')

  else

      job_number=-1

  fi

  echo $job_number

}



# Convert some inputs variables from 3D into their 2D equivalents.

function convert_to_2D {

    # Cylindrical (R-Z) coordinate system.

    geometry_coord_sys=1

    if [ -z "$geometry_is_periodic" ]; then
	geometry_is_periodic=$(get_inputs_var "geometry_is_periodic" $dir)
    fi

    geometry_is_periodic=$(echo $geometry_is_periodic | awk '{print $1, $2}')

    # Set the radial coordinate to have lower boundary value = 0.

    if [ -z "$geometry_prob_lo" ]; then
	geometry_prob_lo=$(get_inputs_var "geometry_prob_lo" $dir)
    fi

    if [ -z "$geometry_prob_hi" ]; then
	geometry_prob_hi=$(get_inputs_var "geometry_prob_hi" $dir)
    fi

    if [ -z "$castro_center" ]; then
	castro_center=$(get_inputs_var "castro_center" $dir)
    fi

    geometry_prob_lo=$(echo $geometry_prob_lo | awk '{print "0.0e0", $2}')
    geometry_prob_hi=$(echo $geometry_prob_hi | awk '{print      $1, $2}')
    castro_center=$(echo $castro_center | awk '{print $1, $2}')

    # Use half as many radial points to keep dr = dz.

    if [ -z "$amr_n_cell" ]; then
	amr_n_cell=$(get_inputs_var "amr_n_cell" $dir)
    fi

    nr=$(echo $amr_n_cell | awk '{print $1}')
    nz=$(echo $amr_n_cell | awk '{print $2}')

    nr=$(echo "$nz / 2" | bc)	

    amr_n_cell="$nr $nz"

    # Use a symmetric lower boundary condition for radial coordinate.

    if [ -z "$castro_lo_bc" ]; then
	castro_lo_bc=$(get_inputs_var "castro_lo_bc" $dir)
    fi

    if [ -z "$castro_hi_bc" ]; then
	castro_hi_bc=$(get_inputs_var "castro_hi_bc" $dir)
    fi

    castro_lo_bc=$(echo $castro_lo_bc | awk '{print  3, $2}')
    castro_hi_bc=$(echo $castro_hi_bc | awk '{print $1, $2}')

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

  # The number of nodes is equal to the number of processors divided 
  # by the number of processors per node. This will not be an integer 
  # if the number of processors requested is not an integer multiple 
  # of the number of processors per node. So we want to use a trick 
  # that guarantees we round upward so that we have enough nodes:
  # if we want the result of (A / B) to be always rounded upward 
  # in integer arithmetic, we evaluate (A + B - 1) / B.

  nodes=$(echo "($nprocs + $ppn - 1) / $ppn" | bc)

  # Number of threads for OpenMP. This will be equal to 
  # what makes the most sense for the machine architecture 
  # by default. For example, the Titan XK7 and Blue Waters XE6
  # Cray nodes are composed of Interlagos boards which are composed 
  # of two NUMA nodes (each NUmA node has 8 integer cores and 4 
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

      if [ $MACHINE == "BLUE_WATERS" ]; then
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

      # We assume that the directory we submitted from is eligible to 
      # work in, so cd to that directory.

      echo "cd \$PBS_O_WORKDIR" >> $dir/$job_script
      echo "" >> $dir/$job_script

      # Load up our helper functions.

      echo "source job_scripts/run_utils.sh" >> $dir/$job_script
      echo "" >> $dir/$job_script

      # Call the function that determines when we're going to stop the run.
      # It should run in the background, to allow the main job to execute.

      echo "check_to_stop &" >> $dir/$job_script
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

      echo "$launcher $launcher_opts $CASTRO $inputs \$(get_restart_string) $redirect" >> $dir/$job_script
      echo "" >> $dir/$job_script

      # With mpirun we redirect the output to a file; let's move that to a file 
      # named by the job number so that we retain it later.

      if [ $launcher == "mpirun" ]; then

	echo "mv $job_name.OU \$(get_last_submitted_job).out" >> $dir/$job_script
	echo "" >> $dir/$job_script

      fi

      # Check to make sure we are done, and if not, re-submit the job.

      if [ -z $no_continue ]; then

	echo "if [ \$(is_dir_done) -ne 1 ]; then" >> $dir/$job_script
	echo "  submit_job" >> $dir/$job_script
	echo "fi" >> $dir/$job_script
	echo "" >> $dir/$job_script

      fi

      # Run the archive script at the end of the simulation.

      if [ $do_storage -ne 1 ]; then
	  echo "do_storage=$do_storage" >> $dir/$job_script
      fi
      echo "archive_all" >> $dir/$job_script
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

    echo "Submitting job in directory "$dir"."

    mkdir -p $dir

    do_job=1

  else

    # First as a sanity check, make sure the desired job isn't already running.

    job_running_status=$(is_job_running $dir)

    if [ $job_running_status -eq 1 ]; then

	echo "Job currently in process or queued in directory $dir."

    else

      # Remove the dump_and_stop file if it exists.

      rm -f $dir/dump_and_stop

      # Archive any files that have not yet been saved to the storage system.

      cwd=$(pwd)

      archive_all $cwd/$dir

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

  # If we are continuing or starting a job, change into the run directory, 
  # submit the job, then come back to the main directory.

  if [ $do_job -eq 1 ]; then

    copy_files $dir
    create_job_script $dir $nprocs $walltime

    # Sometimes we'll want to set up the run directories but not submit them,
    # e.g. for testing purposes, so if the user creates a no_submit variable,
    # we won't actually submit this job.

    if [ -z $no_submit ]; then

	cd $dir

        # Run the job, and capture the job number for output.

	submit_job

	echo "The job number is $(get_last_submitted_job)."

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

        # Build the executable if we haven't yet. Optionally,
	# the user can force a recompile from the run script.

	if [ ! -z $force_recompile ]; then

	    if [ $force_recompile -eq 1 ]; then

		echo "Re-compiling the executable at the user's request."

		cd $compile_dir
		make realclean &> compile.out
		make -j8 $(compile_options) &> compile.out
		cd - > /dev/null

	    fi

	fi

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


########################################################################

# Define variables

# Before we get started, save the list of current shell variables.
# We'll use this to be able to sort out only the ones that have been 
# set by this script or the including script.

shell_list=$(compgen -v)

# Get current machine and set preferences accordingly.
# Note: workdir is the name of the directory you submit 
# jobs from (usually scratch directories).

MACHINE=$(get_machine)

set_machine_params

job_name="wdmerger"
job_script="run_script"

archive_method="none"

do_storage=1

# Directory to compile the executable in

compile_dir="compile"

# Directory for executing and storing results

results_dir="results"

# Directory for placing plots from analysis routines

plots_dir="plots"

# This factor determines the amount of time we want to 
# use as a safety factor in ending jobs.

safety_factor=0.1

if [ -z $inputs ]; then
    inputs=inputs
fi

if [ -z $probin ]; then
    probin=probin
fi

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
