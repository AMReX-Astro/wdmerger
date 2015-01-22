# First, define some functions that are useful.

# This uses the functionality built into the CASTRO makefile setup,
# where make print-$VAR finds the variable VAR in the makefile
# variable list and prints it out to stdout. It is the last word
# on the last line of the make output.

function get_wdmerger_make_var {

    make print-$1 | tail -1 | awk '{ print $NF }'

}

function get_castro_make_var {

    make print-$1 -f $local_makefile  | tail -1 | awk '{ print $NF }'

}


function get_machine {
    
    MACHINE=$(get_make_var WHICHLINUX)

}

# Copies all relevant files needed for a CASTRO run into the target directory.

function copy_files {

    cp $CASTRO $1
    if [ -e "helm_table.dat" ]; then
	cp helm_table.dat $1
    fi
    cp $inputs $1
    cp $probin $1
    cp $job_script $1

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

    echo "Submitting job in directory" $dir

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
	echo "echo \"mpiexec -n $nprocs $CASTRO $inputs > info.out\" | batch" > $job_script
    elif [ $MACHINE == "BLUE_WATERS" ]; then
	sed -i "/#PBS -l nodes/c #PBS -l nodes=$ntasks:ppn=$ppn:xe" $job_script
	sed -i "/#PBS -l walltime/c #PBS -l walltime=$walltime" $job_script
	sed -i "/aprun/c aprun -n $nprocs -N $ppn $CASTRO $inputs \$\{restartString\}" $job_script
    fi

    # Change into the run directory, submit the job, then come back to the main directory.

    copy_files $dir
    cd $dir
    $exec $job_script
    cd - > /dev/null

    # Restore the number of processors per node in case we changed it.

    ppn=$old_ppn

  else

    echo "Directory" $1 "already exists; skipping it."

  fi

}

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
elif [ $MACHINE == "TITAN"]; then
    exec="qsub"
    job_script="titan.run"
    COMP="Cray"
    FCOMP="Cray"
    ppn="8"
fi

if [ ! -e $job_script ]; then
    cp $WDMERGER_HOME/job_scripts/$job_script .
fi

results_dir="results"

if [ ! -d $results_dir ]; then
  mkdir $results_dir
fi

# Create the plots directory, for saving output

plots_dir="plots"

if [ ! -d $plots_dir ]; then
    mkdir $plots_dir
fi
