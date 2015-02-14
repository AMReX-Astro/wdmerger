# First, define some functions that are useful.

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

# Given a directory as the first argument, return the numerically last checkpoint file.

function get_last_checkpoint {

    # The -V option to sort tells it to sort checkpoints in numerically increasing order, which is 
    # necessary since some checkfiles will have five digits and others six. This may only work on
    # certain GNU versions of sort; if it doesn't work for you, perhaps try one of these suggestions:
    # http://stackoverflow.com/questions/7992689/bash-how-to-loop-all-files-in-sorted-order

    checkpoint=$(find $1 -name "chk*" | sort -V | tail -1)

    echo $checkpoint

}

# Copies all relevant files needed for a CASTRO run into the target directory.

function copy_files {

    cp $compile_dir/$CASTRO $1
    if [ -e "$compile_dir/helm_table.dat" ]; then
	cp $compile_dir/helm_table.dat $1
    fi
    cp $compile_dir/$inputs $1
    cp $compile_dir/$probin $1
    cp $compile_dir/$job_script $1

}

# Main submission script. Checks which Linux variant we're on,
# and uses the relevant batch submission script. If you want to
# use a different machine, you'll need to include a run script
# for it in the job_scripts directory.
# The first argument is the name of the directory where we want to 
# submit this job from.
# The second argument is the number of processors you want to run the job on.
# The third argument is the walltime you want the job to run for.
# The last two are optional and default to one processor running for one hour.

function run {

  dir=$1
  if [ ! -z $2 ]; then
      nprocs=$2
  else
      nprocs=1
  fi
  if [ ! -z $3 ]; then
      walltime=$3
  else
      walltime=1:00:00
  fi

  if [ ! -d $dir ]; then

    echo "Submitting job in directory "$dir"."

    mkdir -p $dir

    ntasks=$(expr $nprocs / $ppn)

    # If the number of processors is less than the number of processors per node,
    # there are scaling tests where this is necessary; we'll assume the user understands
    # what they are doing and set it up accordingly.

    old_ppn=$ppn

    if [ $ntasks -eq 0 ]; then
	ntasks="1"
	ppn=$nprocs
    fi

    if [ $MACHINE == "GENERICLINUX" ] ; then 
	echo "echo \"mpiexec -n $nprocs $CASTRO $inputs > info.out\" | batch" > $compile_dir/$job_script
    elif [ $MACHINE == "BLUE_WATERS" ]; then
	sed -i "/#PBS -l nodes/c #PBS -l nodes=$ntasks:ppn=$ppn:xe" $compile_dir/$job_script
	sed -i "/#PBS -l walltime/c #PBS -l walltime=$walltime" $compile_dir/$job_script
	sed -i "/aprun/c aprun -n $nprocs -N $ppn $CASTRO $inputs \$\{restartString\}" $compile_dir/$job_script
    fi

    # Change into the run directory, submit the job, then come back to the main directory.

    copy_files $dir
    cd $dir
    $exec $job_script
    cd - > /dev/null

    # Restore the number of processors per node in case we changed it.

    ppn=$old_ppn

  else

    # If the directory already exists, check to see if we've reached the desired stopping point.

    checkpoint=$(get_last_checkpoint)

    # Extract the checkpoint time. It is stored in row 3 of the Header file.

    time=$(awk 'NR==3' $checkpoint/Header)

    # Extract the current timestep. It is stored in row 12 of the Header file.

    step=$(awk 'NR==12' $checkpoint/Header)

    # Now determine if we are both under max_step and stop_time. If so, re-submit the job.
    # The job script already knows to start from the latest checkpoint file.

    stop_time=$(grep "stop_time" $dir/$inputs | awk '{print $3}')
    max_step=$(grep "max_step" $dir/$inputs | awk '{print $3}')

    time_flag=$(echo "$time < $stop_time" | bc -l)
    step_flag=$(echo "$step < $max_step" | bc -l)

    # First as a sanity check, make sure the desired job isn't already running.

    if [ -e $dir/*$run_ext ]; then

	echo "Job currently in process in directory "$dir"."

    else
 
	if [ $time_flag -eq 1 ] && [ $step_flag -eq 1 ]; then

	echo "Continuing job in directory "$dir"."

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

# Directory to compile the executable in

compile_dir="compile"

# Upon initialization, store some variables and create results directory.

local_makefile=$(get_wdmerger_make_var local_makefile)
inputs=$(get_wdmerger_make_var inputs)
probin=$(get_wdmerger_make_var probin)
CASTRO=$(get_castro_make_var executable)

# This returns the Linux variant the current machine uses, defined in
# $(BOXLIB_HOME)/Tools/C_mk/Make.defs.
# This is used to select which batch submission script we want to use.

MACHINE=$(get_castro_make_var WHICHLINUX)

if [ $MACHINE == "GENERICLINUX" ]; then
    exec="bash"
    job_script="linux.run"
    ppn="16"
elif [ $MACHINE == "BLUE_WATERS" ]; then
    exec="qsub"
    job_script="bluewaters.run"
    COMP="Cray"
    FCOMP="Cray"
    ppn="16"
    run_ext=".OU"
elif [ $MACHINE == "TITAN"]; then
    exec="qsub"
    job_script="titan.run"
    COMP="Cray"
    FCOMP="Cray"
    ppn="8"
fi

if [ ! -e $compile_dir/$job_script ]; then
    cp $WDMERGER_HOME/job_scripts/$job_script $compile_dir
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

