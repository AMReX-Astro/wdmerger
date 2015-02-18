# run_utils.sh: helper functions for job submission scripts.

# This uses the functionality built into the CASTRO makefile setup,
# where make print-$VAR finds the variable VAR in the makefile
# variable list and prints it out to stdout. It is the last word
# on the last line of the make output.

function get_wdmerger_make_var {

    make print-$1 -C $compile_dir | tail -2 | head -1 | awk '{ print $NF }'

}

function get_castro_make_var {

    make print-$1 -f $local_makefile -C $compile_dir  | tail -2 | head -1 | awk '{ print $NF }'

}



# Return the name of the machine we're currently working on.

function get_machine {

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

  echo $MACHINE

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
    # with six digits, and then will fall back to ones with five digits.
    # On recent versions of GNU sort, one can simplify this with sort -V.

    checkpoint=$(find $dir -type d -name "*chk??????" | sort | tail -1)

    if [ -z $checkpoint ]; then

        checkpoint=$(find $dir -type d -name "*chk?????"  | sort | tail -1)

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



# Archive the file or directory given in the first argument, 
# to the same path on the archive machine relative to the machine's $workdir.

function archive {

  if [ ! -z $1 ]; then
      if [ -d $1 ]; then
	  echo "Archiving directory "$1"."
      else
	  echo "Archiving file " $1"."
      fi
  else
      echo "No file to archive; exiting."
      return
  fi

  # We may get a directory to archive, so call basename to make sure $file
  # doesn't appear with a trailing slash.

  file=$(basename $1)
  dir=$(dirname $1)

  # Remove the $workdir from the name. We test on both $workdir and $workdir/
  # so that the trailing slash doesn't matter here.

  storage_dir=$(echo ${dir#$workdir/})
  storage_dir=$(echo ${dir#$workdir})

  # Archive based on the method chosen for this machine.

  if   [ $archive_method == "htar" ]; then

      $HTAR ${storage_dir}/${file}.tar $dir/$file

  elif [ $archive_method == "globus" ]; then

      archive_dir=/projects/sciteam/$allocation/$USER/$storage_dir

      src=$globus_src_endpoint$dir/$file
      dst=$globus_dst_endpoint$archive_dir/$file

      if [ -d $dir/$file ]; then
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

  archivelist=""

  # Archive the plotfiles and checkpoint files.
  # Make sure that they have been completed by checking if
  # the Header file exists, which is the last thing created.
  # If there are no plotfiles or checkpoints to archive,
  # then assume we have completed the run and exit.

  pltlist=$(find $dir -maxdepth 1 -type d -name "*plt*" | sort)
  chklist=$(find $dir -maxdepth 1 -type d -name "*chk*" | sort)

  if [[ -z $pltlist ]] && [[ -z $chklist ]]; then
      return
  fi

  # Move all completed plotfiles and checkpoints to the output
  # directory, and add them to the list of things to archive.

  for file in $pltlist
  do
      if [ -e $file/Header ]; then
	  mv $file $dir/output/
	  f=$(basename $file)
	  archivelist=$archivelist" "$f
      fi
  done

  for file in $chklist
  do
      if [ -e $file/Header ]; then
	  mv $file $dir/output/
	  f=$(basename $file)
	  archivelist=$archivelist" "$f
      fi
  done

  diaglist=$(find $dir -maxdepth 1 -name "*diag*.out")

  # For the diagnostic files, we just want to make a copy and move it to the 
  # output directory; we can't move it, since the same file needs to be there
  # for the duration of the simulation if we want a continuous record.

  for file in $diaglist
  do
      cp $file $dir/output/
      f=$(basename $file)
      archivelist=$archivelist" "$f
  done

  # Same thing for the runtime stdout files.

  outlist=$(find $dir -maxdepth 1 -name *$job_name*)

  for file in $outlist
  do
      cp $file $dir/output/
      f=$(basename $file)
      archivelist=$archivelist" "$f
  done

  # Now we'll do the archiving for all files in $archivelist.
  # Determine the archiving method based on machine.

  if   [ $MACHINE == "TITAN"       ]; then

      # For Titan, just loop over every file we're archiving and htar it.

      for file in $archivelist
      do
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



# Copies all relevant files needed for a CASTRO run into the target directory.

function copy_files {

    cp $compile_dir/$CASTRO $1
    if [ -e "$compile_dir/helm_table.dat" ]; then
	cp $compile_dir/helm_table.dat $1
    fi
    cp $compile_dir/$inputs $1
    cp $compile_dir/$probin $1

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

  # If the number of processors is less than the number of processors per node,
  # there are scaling tests where this is necessary; we'll assume the user understands
  # what they are doing and set it up accordingly.

  old_ppn=$ppn

  if [ $nodes -eq 0 ]; then
      nodes="1"
      ppn=$nprocs
  fi

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

      # Number of threads for OpenMP

      echo "export OMP_NUM_THREADS=$OMP_NUM_THREADS" >> $dir/$job_script

      restartString=$(get_restart_string $dir)

      echo "aprun -n $nprocs -N $ppn $CASTRO $inputs $restartString" >> $dir/$job_script

   elif [ $batch_system == "batch" ]; then

      echo "echo \"mpiexec -n $nprocs $CASTRO $inputs > info.out\" | batch" > $dir/$job_script

   fi

   # Restore the number of processors per node in case we changed it.

   ppn=$old_ppn

}



# Main submission script. Checks which Linux variant we're on,
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

  if [ ! -d $dir ]; then

    echo "Submitting job in directory "$dir"."

    mkdir -p $dir

    # Change into the run directory, submit the job, then come back to the main directory.

    copy_files $dir
    create_job_script $dir $nprocs $walltime
    cd $dir
    $exec $job_script
    cd - > /dev/null

  else

    # First as a sanity check, make sure the desired job isn't already running.

    if [ -e $dir/*$run_ext ]; then

	echo "Job currently in process in directory "$dir"."

    else

      # Archive any files that have not yet been saved to the storage system.

      cwd=$(pwd)

      archive_all $cwd/$dir

      # If the directory already exists, check to see if we've reached the desired stopping point.

      checkpoint=$(get_last_checkpoint $dir)

      time_flag=1
      step_flag=1

      if [ -e $dir/$checkpoint/Header ]; then

	  # Extract the checkpoint time. It is stored in row 3 of the Header file.

	  time=$(awk 'NR==3' $dir/$checkpoint/Header)

	  # Extract the current timestep. It is stored in row 12 of the Header file.

	  step=$(awk 'NR==12' $dir/$checkpoint/Header)

	  # Now determine if we are both under max_step and stop_time. If so, re-submit the job.
	  # The job script already knows to start from the latest checkpoint file.

	  stop_time=$(grep "stop_time" $dir/$inputs | awk '{print $3}')
	  max_step=$(grep "max_step" $dir/$inputs | awk '{print $3}')

	  time_flag=$(echo "$time < $stop_time" | bc -l)
	  step_flag=$(echo "$step < $max_step" | bc -l)

      fi

      if [ $time_flag -eq 1 ] && [ $step_flag -eq 1 ]; then

	  echo "Continuing job in directory "$dir"."

	  create_job_script $dir $nprocs $walltime
	  cd $dir
	  $exec $job_script
	  cd - > /dev/null

	  # If we make it here, then we've already reached either stop_time
	  # or max_step, so we should conclude that the run is done.

      else

	  echo "Job has already been completed in directory "$dir"."

      fi
    fi

  fi

}



########################################################################

# Define variables

# Get current machine and set preferences accordingly.
# Note: workdir is the name of the directory you submit 
# jobs from (usually scratch directories).

MACHINE=$(get_machine)

job_name="wdmerger"
job_script="run_script"

OMP_NUM_THREADS="1"

archive_method="none"

if [ $MACHINE == "GENERICLINUX" ]; then

    exec="bash"
    ppn="16"
    batch_system="batch"

elif [ $MACHINE == "BLUE_WATERS" ]; then

    allocation="jni"
    exec="qsub"
    COMP="Cray"
    FCOMP="Cray"
    ppn="16"
    node_type="xe"
    run_ext=".OU"
    workdir="/scratch/sciteam/$USER/"
    batch_system="PBS"
    archive_method="globus"
    globus_src_endpoint="ncsa#BlueWaters"
    globus_dst_endpoint="ncsa#Nearline"

elif [ $MACHINE == "TITAN" ]; then

    allocation="ast106"
    exec="qsub"
    COMP="Cray"
    FCOMP="Cray"
    ppn="8"
    run_ext=".OU"
    workdir="/lustre/atlas/scratch/$USER/$allocation/"
    batch_system="PBS"
    archive_method="htar"

elif [ $MACHINE == "HOPPER" ]; then
    
    allocation="m1400"
    exec="qsub"
    COMP="Cray"
    FCOMP="Cray"
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

# Upon initialization, store some variables and create results directory.
# We only want to initialize these variables if we're currently in a root problem directory.

if [ -d $compile_dir ]; then

    if [ -e $compile_dir/makefile ]; then

	local_makefile=$(get_wdmerger_make_var local_makefile)
	inputs=$(get_wdmerger_make_var inputs)
	probin=$(get_wdmerger_make_var probin)
	CASTRO=$(get_castro_make_var executable)

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
