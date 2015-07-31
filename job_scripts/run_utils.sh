# run_utils.sh: helper functions for job submission scripts.

# Bring in scripts for working with inputs and probin files.
source $WDMERGER_HOME/job_scripts/inputs.sh
source $WDMERGER_HOME/job_scripts/probin.sh

# Basic mathematical routines.
source $WDMERGER_HOME/job_scripts/math.sh

# This uses the functionality built into the CASTRO makefile setup,
# where make print-$VAR finds the variable VAR in the makefile
# variable list and prints it out to stdout. It is the last word
# on the last line of the make output.

function get_make_var {

    make print-$1 -C $compile_dir | tail -2 | head -1 | awk '{ print $NF }'

}



# Return the name of the machine we're currently working on.

function get_machine {

  # Get the name of the machine by using uname;
  # then store it in a file in the wdmerger root.
  # This storage helps when we're on the compute nodes,
  # which often don't have the same system name as the login nodes.

  if [ ! -e $WDMERGER_HOME/job_scripts/machine ]; then

      UNAMEN=$(uname -n)

      if   [[ $UNAMEN == *"h2o"*    ]]; then
	MACHINE=BLUE_WATERS
      elif [[ $UNAMEN == *"titan"*  ]]; then
	MACHINE=TITAN
      elif [[ $UNAMEN == *"hopper"* ]]; then
	MACHINE=HOPPER
      else
	MACHINE=GENERICLINUX
      fi

      # Store the name 

      echo "$MACHINE" > $WDMERGER_HOME/job_scripts/machine

  else

      MACHINE=$(cat $WDMERGER_HOME/job_scripts/machine)

  fi

  echo $MACHINE

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

    output=$(find $dir -name "$job_name*" | sort | tail -1)

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
	jobs_in_directory=$(cat $dir/jobs_submitted.txt)

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
      return
  fi

  # If the directory already exists, check to see if we've reached the desired stopping point.
  # There are two places we can look: in the last checkpoint file, or in the last stdout file. 

  checkpoint=$(get_last_checkpoint $dir)
  last_output=$(get_last_output $dir)

  # Get the desired stopping time and max step from the inputs file in the directory.
  # Alternatively, we may have set this from the calling script, so prefer that.

  if [ -z $stop_time ]; then
      stop_time=$(get_inputs_var "stop_time" $dir)
  fi

  if [ -z $max_step ]; then
      max_step=$(get_inputs_var "max_step" $dir)
  fi

  # Assume we're done, by default.

  time_flag=1
  step_flag=1

  done_status=0

  if [ -e $dir/$checkpoint/Header ]; then

      # Extract the checkpoint time. It is stored in row 3 of the Header file.

      chk_time=$(awk 'NR==3' $dir/$checkpoint/Header)

      # Extract the current timestep. We can get it from the 
      # name of the checkpoint file. cut will do the trick;
      # just capture everything after the 'k' of 'chk'.

      chk_step=$(echo $checkpoint | cut -d"k" -f2)

      time_flag=$(echo "$chk_time >= $stop_time" | bc -l)
      step_flag=$(echo "$chk_step >= $max_step" | bc)

  elif [ -e $dir/$last_output ]; then

      output_time=$(grep "STEP =" $dir/$last_output | tail -1 | awk '{print $6}')
      output_step=$(grep "STEP =" $dir/$last_output | tail -1 | awk '{print $3}')

      # bc can't handle numbers in scientific notation, so use printf to convert it to floating point.

      output_time=$(printf "%f" $output_time)

      time_flag=$(echo "$output_time >= $stop_time" | bc -l)
      step_flag=$(echo "$output_step >= $max_step" | bc)

  fi

  if [ $time_flag -eq 1 ] || [ $step_flag -eq 1 ]; then
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
      echo "No directory passed to function archive_all; exiting."
      return
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

  inputs_list=$(find $directory -maxdepth 1 -name "*inputs*")

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
	    echo $directory/output/$file
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
            cp $WDMERGER_HOME/source/inputs $dir/$inputs
	fi
    fi

    if [ ! -e "$dir/$probin" ]; then
	if [ -e "$compile_dir/$probin" ]; then
	    cp $compile_dir/probin $dir/$probin
	else
	    cp $WDMERGER_HOME/source/probin $dir/$probin
	fi
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

  nodes=$(expr $nprocs / $ppn)

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

      OMP_NUM_THREADS=$threads_per_task

      # First, get the maximum grid size. If this has been
      # set by the including script, we use that; otherwise, 
      # we read in the value from the main inputs file.

      if [ -z "$amr_max_grid_size" ]; then
	  amr_max_grid_size=$(get_inputs_var "amr_max_grid_size")
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
      echo "#PBS -A $allocation" >> $dir/$job_script

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
      else
	  echo "#PBS -l nodes=$nodes" >> $dir/$job_script
      fi

      # Queue to submit to. This is required for some systems.
      if [ ! -z $queue ]; then
	  echo "#PBS -q $queue" >> $dir/$job_script
      fi

      echo "" >> $dir/$job_script

      # We assume that the directory we submitted from is eligible to 
      # work in, so cd to that directory.

      echo "cd \$PBS_O_WORKDIR" >> $dir/$job_script

      # Tell the script where the root wdmerger directory is, since 
      # environment variables are generally not loaded onto compute nodes.

      echo "export WDMERGER_HOME=$WDMERGER_HOME" >> $dir/$job_script
      echo "source $WDMERGER_HOME/job_scripts/run_utils.sh" >> $dir/$job_script
         
      # Set up the script that periodically checks whether we should
      # terminate the run. Only do this for jobs with large enough
      # processor counts. Run the process in the background.

      if [ $nodes -ge 8 ]; then
	  echo "bash $WDMERGER_HOME/job_scripts/check_to_stop.sh . &" >> $dir/$job_script
      fi

      # Number of OpenMP threads

      echo "export OMP_NUM_THREADS=$OMP_NUM_THREADS" >> $dir/$job_script

      # Set the aprun options.

      aprun_opts="-n $num_mpi_tasks -N $tasks_per_node -d $OMP_NUM_THREADS"

      echo "restartString=\$(get_restart_string .)" >> $dir/$job_script

      # Main job execution.

      echo "" >> $dir/$job_script
      echo "aprun $aprun_opts $CASTRO $inputs \$restartString" >> $dir/$job_script
      echo "" >> $dir/$job_script

      # Run the archive script at the end of the simulation.

      echo "archive_all ." >> $dir/$job_script

      # Check to make sure we are done, and if not, re-submit the job.

      if [ -z $do_chain ]; then

	echo "" >> $dir/$job_script
	echo "dir=." >> $dir/$job_script
	echo "done_flag=\$(is_dir_done)" >> $dir/$job_script
	echo "if [ \$done_flag -ne 1 ]; then" >> $dir/$job_script
	echo "  job_number=`$exec $job_script`" >> $dir/$job_script
	echo "  job_number=\${job_number%%.*}" >> $dir/$job_script
	echo "  echo $job_number >> jobs_submitted.txt" >> $dir/$job_script
	echo "fi" >> $dir/$job_script
	echo "" >> $dir/$job_script

     fi

   elif [ $batch_system == "batch" ]; then

      echo "echo \"mpiexec -n $nprocs $CASTRO $inputs > $job_name$run_ext\" | batch" > $dir/$job_script

   fi

   # Restore the number of processors per node in case we changed it.

   ppn=$old_ppn

}



# Wrapper script for run. The first argument is an integer N that
# causes us to divide $stop_time into N equal increments.

function chain {

  if [ -z $N_iters ]; then
      echo "N_iters not set in call to chain; exiting."
      return
  fi

  if [ -z $stop_time ]; then
      echo "stop_time not defined in call to chain; exiting."
  fi

  if [ -z $dir ]; then
      echo "No directory given to chain; exiting."
      return
  fi

  orig_stop_time=$stop_time
  orig_inputs=$inputs
  orig_probin=$probin

  do_chain=1

  # First check to see if we've completed the run yet.

  done_flag=0

  if [ -d $dir ]; then
      inputs=$(find $dir -maxdepth 1 -name "inputs_*" | sort -n | tail -1)

      # If there are no inputs files in the directory, we know we haven't yet started.

      if [ ! -z $inputs ]; then
          inputs=$(basename $inputs)
      fi

      done_flag=$(is_dir_done)
  else
      mkdir -p $dir
  fi

  # It is possible that we are continuing or extending a run.
  # Check the current time from the last check point, and 
  # only submit the jobs that are remaining to do.

  checkpoint=$(get_last_checkpoint $dir)

  if [ ! -z $checkpoint ]; then
      chk_time=$(awk 'NR==3' $dir/$checkpoint/Header)
  else
      chk_time=0.0
  fi

  job_running_status=$(is_job_running $dir)


  if [ $done_flag -ne 1 ] && [ $job_running_status -ne 1 ]; then

      for N in $(seq 1 $N_iters)
      do
	  inputs="inputs_"$N
	  stop_time=$(echo "$orig_stop_time * $N / $N_iters" | bc -l)
	  run_test=$(echo "$stop_time > $chk_time" | bc -l)
	  if [ $run_test -eq 1 ]; then
	      run
	      job_dependency=$job_number
          fi
      done

  else

      if [ $done_flag -eq 1 ]; then

	  echo "Chain completed in directory "$dir$"."

      elif [ $job_running_status -eq 1 ]; then

	  echo "Chain currently queued or in process in directory "$dir"."

      fi

  fi

  job_dependency=""
  stop_time=$orig_stop_time
  inputs=$orig_inputs
  probin=$orig_probin
  do_chain=""

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

  if [ -z $nprocs ]; then
      nprocs=$ppn
  fi

  if [ -z $walltime ]; then
      walltime=1:00:00
  fi

  if [ -z $inputs ]; then
      inputs=inputs
  fi

  if [ -z $probin ]; then
      probin=probin
  fi

  do_job=0

  if [ ! -d $dir ]; then

    echo "Submitting job in directory "$dir"."

    mkdir -p $dir

    do_job=1

  elif [ ! -z $do_chain ]; then

      echo "Continuing chain in directory $dir to time $stop_time."

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

    cd $dir

    exec_command=$exec

    if [ ! -z $job_dependency ]; then
	if [ $batch_system == "PBS" ]; then
	    exec_command="$exec -W depend=afterok:$job_dependency"
	fi
    fi

    # Capture the job number output.

    job_number=`$exec_command $job_script`

    # Some systems like Blue Waters include the system name
    # at the end of the number, so remove it.

    job_number=${job_number%%.*}

    echo "The job number is $job_number."

    echo "$job_number" >> jobs_submitted.txt

    cd - > /dev/null

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

job_name="wdmerger"
job_script="run_script"

# Determine the number of OpenMP threads per task to run 
# with, by default. We'll update this for the various machines
# given their configurations, and then in the run script 
# depending on the size of the problem, but the user 
# can overwrite these defaults by setting OMP_NUM_THREADS
# in the including script.

threads_per_task="1"

archive_method="none"

do_storage=1

if [ $MACHINE == "GENERICLINUX" ]; then

    exec="bash"
    ppn="16"
    batch_system="batch"
    run_ext=".OU"

elif [ $MACHINE == "BLUE_WATERS" ]; then

    allocation="jni"
    exec="qsub"
    ppn="32"
    threads_per_task="8"
    node_type="xe"
    run_ext=".OU"
    batch_system="PBS"
    archive_method="globus"
    globus_src_endpoint="ncsa#BlueWaters"
    globus_dst_endpoint="ncsa#Nearline/projects/sciteam/$allocation/$USER"

elif [ $MACHINE == "TITAN" ]; then

    allocation="ast106"
    exec="qsub"
    ppn="16"
    threads_per_task="8"
    run_ext=".OU"
    batch_system="PBS"
    archive_method="htar"

elif [ $MACHINE == "HOPPER" ]; then
    
    allocation="m1400"
    exec="qsub"
    ppn="24"
    run_ext=".OU"
    batch_system="PBS"
    queue="regular"

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


# Directory to compile the executable in

compile_dir="compile"

# Some variables we need for storing job information.

num_jobs=
job_arr=()
state_arr=()
walltime_arr=()

# Fill these arrays.

get_submitted_jobs

# Upon initialization, store some variables and create results directory.
# We only want to initialize these variables if we're currently in a root problem directory.

if [ -d $compile_dir ]; then

    if [ -e $compile_dir/GNUmakefile ]; then

	CASTRO=$(get_make_var executable)

    fi

    # Directory for executing and storing results

    results_dir="results"

    if [ ! -d $results_dir ]; then
      mkdir $results_dir
    fi

    # Directory for placing plots from analysis routines

    plots_dir="plots"

    if [ ! -d $plots_dir ]; then
	mkdir $plots_dir
    fi

fi
