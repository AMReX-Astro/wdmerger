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
      elif [[ $UNAMEN == *"lired"*  ]]; then
	MACHINE=LIRED
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


# Set variables that depend on which system we're using.

function set_machine_params {

    # The default choice.

    if [ $MACHINE == "GENERICLINUX" ]; then

	exec="bash"
	ppn="16"
	batch_system="batch"
	launcher="aprun"
	run_ext=".OU"

    # Blue Waters at NCSA

    elif [ $MACHINE == "BLUE_WATERS" ]; then

	allocation="jni"
	exec="qsub"
	cancel_job="qdel"
	ppn="32"
	threads_per_task="8"
	node_type="xe"
	run_ext=".OU"
	batch_system="PBS"
	launcher="aprun"
	archive_method="globus"
	globus_src_endpoint="ncsa#BlueWaters"
	globus_dst_endpoint="ncsa#Nearline/projects/sciteam/$allocation/$USER"
	time_remaining_column="5"

    # Titan at OLCF

    elif [ $MACHINE == "TITAN" ]; then

	allocation="ast106"
	exec="qsub"
	cancel_job="qdel"
	ppn="16"
	threads_per_task="8"
	run_ext=".OU"
	batch_system="PBS"
	launcher="aprun"
	archive_method="htar"
	time_remaining_column="5"

    # Hopper at NERSC

    elif [ $MACHINE == "HOPPER" ]; then
	
	allocation="m1400"
	exec="qsub"
	cancel_job="qdel"
	ppn="24"
	run_ext=".OU"
	batch_system="PBS"
	launcher="aprun"
	queue="regular"

    # LIRED at Stony Brook University

    elif [ $MACHINE == "LIRED" ]; then

	exec="qsub"
	cancel_job="qdel"
	ppn="24"
	threads_per_task="1"
	batch_system="PBS"
	launcher="mpirun"
	queue="medium"
	run_ext=".OU"
	time_remaining_column="5"

    fi

}