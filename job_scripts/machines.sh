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


# Set variables that depend on which system we're using.

function set_machine_params {

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

}