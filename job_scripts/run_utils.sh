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

    if [ -e $dir/*$run_ext ]; then
	echo "1"
    else
	echo "0"
    fi

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

  # Remove everything from the directory up to the username. The assumption here is that
  # everything after that was created by the user, and that's the directory structure we want
  # to preserve when moving things over to the storage system.

  storage_dir=$d

  storage_dir=${storage_dir#*$USER/}
  storage_dir=${storage_dir#*$USER}

  # Archive based on the method chosen for this machine.

  if   [ $archive_method == "htar" ]; then

      $HTAR ${storage_dir}/${f}.tar $d/$f

  elif [ $archive_method == "globus" ]; then

      archive_dir=$allocation/$USER/$storage_dir

      src=$globus_src_endpoint/$d/$f
      dst=$globus_dst_endpoint/$archive_dir/$f

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
      dir=$1
  else
      echo "No directory passed to function archive_all; exiting."
      return
  fi

  if [ ! -d $dir/output/ ]; then
      mkdir $dir/output/
  fi

  archivelist=""

  # Archive the plotfiles and checkpoint files.
  # Make sure that they have been completed by checking if
  # the Header file exists, which is the last thing created.

  pltlist=$(find $dir -maxdepth 1 -type d -name "*plt*" | sort)
  chklist=$(find $dir -maxdepth 1 -type d -name "*chk*" | sort)

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
	  mv $file $dir/output/
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
	  mv $file $dir/output/
	  f=$(basename $file)
	  archivelist=$archivelist" "$f
      fi
  done

  diaglist=$(find $dir -maxdepth 1 -name "*diag*.out")

  # For the diagnostic files, we just want to make a copy and move it to the 
  # output directory; we can't move it, since the same file needs to be there
  # for the duration of the simulation if we want a continuous record. But 
  # we want to avoid archiving the files again if the run has already been
  # completed, so we check the timestamps and only move the file to the output
  # directory if the archived version is older than the main version.

  for file in $diaglist
  do
      diag_basename=$(basename $file)
      if [ -e $dir/output/$diag_basename ]; then
	  if [ $dir/output/$diag_basename -nt $file ]; then
	      continue
	  fi
      fi
      cp $file $dir/output/
      f=$(basename $file)
      archivelist=$archivelist" "$f
  done

  # Same thing for the runtime stdout files.

  outlist=$(find $dir -maxdepth 1 -name *$job_name*)

  for file in $outlist
  do
      output_basename=$(basename $file)
      if [ -e $dir/output/$output_basename ]; then
	  if [ $dir/output/$output_basename -nt $file ]; then
	      continue
	  fi
      fi
      cp $file $dir/output/
      f=$(basename $file)
      archivelist=$archivelist" "$f
  done

  # Same strategy for the inputs and probin files.

  if ([ ! -e $dir/output/inputs ] || [ $dir/output/inputs -ot $dir/inputs ]); then
      cp $dir/inputs $dir/output/
      archivelist=$archivelist" "inputs
  fi

  if ([ ! -e $dir/output/probin ] || [ $dir/output/probin -ot $dir/probin ]); then
      cp $dir/probin $dir/output/
      archivelist=$archivelist" "probin
  fi

  # If there is nothing to archive,
  # then assume we have completed the run and exit.

  if [[ -z $archivelist ]]; then
      return
  fi

  # Now we'll do the archiving for all files in $archivelist.
  # Determine the archiving method based on machine.

  if   [ $MACHINE == "TITAN"       ]; then

      # For Titan, just loop over every file we're archiving and htar it.

      for file in $archivelist
      do
	  echo $dir/output/$file
	  archive $dir/output/$file
      done

  elif [ $MACHINE == "BLUE_WATERS" ]; then

      # For Blue Waters, we're using Globus Online, which has a cap on the number 
      # of simultaneous transfers a user can have. Therefore our strategy is
      # to sync the entire output directory of this location rather than 
      # transferring the files independently.

      archive $dir/output/

  fi

}



# Copies all relevant files needed for a CASTRO run into the target directory,
# and updates the inputs and probin according to any shell variables we set.

function copy_files {

    if [ -z $1 ]
    then
	echo "No directory passed to copy_files; exiting."
    return
    else
	dir=$1
    fi

    if [ ! -e $dir/$CASTRO ]; then
        cp $compile_dir/$CASTRO $dir
    fi

    if [ ! -e "$dir/helm_table.dat" ]; then
	if [ -e "$compile_dir/helm_table.dat" ]; then
	    cp $compile_dir/helm_table.dat $dir
	fi
    fi

    if [ ! -e "$dir/inputs" ]; then
        if [ -e "$source_dir/inputs" ]; then
            cp $source_dir/inputs $dir
        else
            cp $WDMERGER_HOME/source/inputs $dir
	fi
    fi

    if [ ! -e "$dir/probin" ]; then
	if [ -e "$source_dir/probin" ]; then
	    cp $source_dir/probin $dir
	else
	    cp $WDMERGER_HOME/source/probin $dir
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

	replace_inputs_var $var $dir
	replace_probin_var $var $dir

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

      # We assume that the directory we submitted from is eligible to 
      # work in, so cd to that directory.

      echo "cd \$PBS_O_WORKDIR" >> $dir/$job_script

      # Number of OpenMP threads

      echo "export OMP_NUM_THREADS=$OMP_NUM_THREADS" >> $dir/$job_script

      # Set the aprun options.

      aprun_opts="-n $num_mpi_tasks -N $tasks_per_node -d $OMP_NUM_THREADS"

      restartString=$(get_restart_string $dir)

      # Main job execution.

      echo "aprun $aprun_opts $CASTRO inputs $restartString" >> $dir/$job_script

   elif [ $batch_system == "batch" ]; then

      echo "echo \"mpiexec -n $nprocs $CASTRO inputs > $job_name$run_ext\" | batch" > $dir/$job_script

   fi

   # Restore the number of processors per node in case we changed it.

   ppn=$old_ppn

}



# Main submission function. Checks which Linux variant we're on,
# and uses the relevant batch submission script. If you want to
# use a different machine, you'll need to include a run script
# for it in the job_scripts directory.
# The first argument is the name of the directory where we want to 
# submit this job from.
# The second argument is the number of processors you want to run the job on.
# The third argument is the walltime you want the job to run for.
# The last two are optional and default to one node running for one hour.

function run {

  if [ ! -z $1 ]; then
      dir=$1
  else
      echo "No directory given to run; exiting."
      return
  fi

  if [ ! -z $2 ]; then
      nprocs=$2
  else
      nprocs=$ppn
  fi

  if [ ! -z $3 ]; then
      walltime=$3
  else
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

	echo "Job currently in process in directory "$dir"."

    else

      # Remove the dump_and_stop file if it exists.

      rm -f $dir/dump_and_stop

      # Archive any files that have not yet been saved to the storage system.

      cwd=$(pwd)

      archive_all $cwd/$dir

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

      # Assume we're continuing, by default.

      time_flag=1
      step_flag=1

      if [ -e $dir/$checkpoint/Header ]; then

	  # Extract the checkpoint time. It is stored in row 3 of the Header file.

	  chk_time=$(awk 'NR==3' $dir/$checkpoint/Header)

	  # Extract the current timestep. We can get it from the 
	  # name of the checkpoint file. cut will do the trick;
	  # just capture everything after the 'k' of 'chk'.

	  chk_step=$(echo $checkpoint | cut -d"k" -f2)

	  time_flag=$(echo "$chk_time < $stop_time" | bc -l)
	  step_flag=$(echo "$chk_step < $max_step" | bc)

      elif [ -e $dir/$last_output ]; then

	  output_time=$(grep "STEP =" $dir/$last_output | tail -1 | awk '{print $6}')
	  output_step=$(grep "STEP =" $dir/$last_output | tail -1 | awk '{print $3}')

	  # bc can't handle numbers in scientific notation, so use printf to convert it to floating point.

	  output_time=$(printf "%f" $output_time)

	  time_flag=$(echo "$output_time < $stop_time" | bc -l)
	  step_flag=$(echo "$output_step < $max_step" | bc)

      fi

      if [ $time_flag -eq 1 ] && [ $step_flag -eq 1 ]; then

	  echo "Continuing job in directory "$dir"."

	  do_job=1

      else

	  # If we make it here, then we've already reached either stop_time
	  # or max_step, so we should conclude that the run is done.

	  echo "Job has already been completed in directory "$dir"."

      fi
    fi

  fi

  # If we are continuing or starting a job, change into the run directory, 
  # submit the job, then come back to the main directory.

  if [ $do_job -eq 1 ]; then

    copy_files $dir
    create_job_script $dir $nprocs $walltime
    cd $dir
    $exec $job_script
    cd - > /dev/null

    # Set up the function that periodically checks whether we should
    # terminate the run. Only do this for jobs with large enough
    # processor counts. Run the process in the background and disown it.

    nohup bash $WDMERGER_HOME/job_scripts/check_to_stop.sh $dir >&/dev/null &

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
    globus_dst_endpoint="ncsa#Nearline/projects/sciteam"

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
source_dir="source"

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
