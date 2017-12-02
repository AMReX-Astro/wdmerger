# Return the name of the machine we're currently working on.

function get_machine {

  # The generic default.

  MACHINE=GENERICLINUX

  # Check to see if the machine name has been cached
  # in a local job_scripts run directory; if so, grab
  # it from there (this is useful on compute nodes at
  # various HPC sites that don't name the nodes with
  # a useful identifier).

  # Otherwise, assume we are on a login node and so
  # we can identify the machine based on the uname.

  if [ -d "job_scripts" ]; then

      if [ -e "job_scripts/machine" ]; then
          MACHINE=$(cat "job_scripts/machine")
      fi

  else

      UNAMEN=$(uname -n)

      if   [[ $UNAMEN == *"h2o"*    ]]; then
          MACHINE=BLUE_WATERS
      elif [[ $UNAMEN == *"titan"*  ]]; then
          MACHINE=TITAN
      elif [[ $UNAMEN == *"hopper"* ]]; then
          MACHINE=HOPPER
      elif [[ $UNAMEN == *"edison"* ]]; then
          MACHINE=EDISON
      elif [[ $UNAMEN == *"cori"* ]]; then
          MACHINE=CORI
      elif [[ $UNAMEN == *"lired"*  ]]; then
          MACHINE=LIRED
      elif [[ $UNAMEN == "login"  ]]; then
          MACHINE=SEAWULF
      elif [[ $UNAMEN == *"mira"*   ]]; then
          MACHINE=MIRA
      fi

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
	pause_job="qhold"
	resume_job="qrls"
	ppn="32"
	threads_per_task="8"
	node_type="xe"
	run_ext=".OU"
	batch_system="PBS"
	launcher="aprun"
	archive_method="globus"
	globus_src_endpoint="ncsa#BlueWaters"
	globus_dst_endpoint="ncsa#Nearline/projects/sciteam/$allocation/$USER"

    # Titan at OLCF

    elif [ $MACHINE == "TITAN" ]; then

	allocation="ast106sbu"
	exec="qsub"
	cancel_job="qdel"
	pause_job="qhold"
	resume_job="qrls"
	ppn="16"
	threads_per_task="8"
	run_ext=".OU"
	batch_system="PBS"
	queue="batch"
	launcher="aprun"
	archive_method="htar"
	archive_queue="dtn"
	archive_wclimit="24:00:00"

    # Hopper at NERSC

    elif [ $MACHINE == "HOPPER" ]; then
	
	allocation="m1400"
	exec="qsub"
	cancel_job="qdel"
	pause_job="qhold"
	resume_job="qrls"
	ppn="24"
	run_ext=".OU"
	batch_system="PBS"
	launcher="aprun"
	queue="regular"

    # Edison at NERSC

    elif [ $MACHINE == "EDISON" ]; then

	allocation="m1938"
	exec="sbatch"
	cancel_job="scancel"
	ppn="24"
	run_ext=".OU"
	batch_system="SLURM"
	launcher="srun"
	queue="regular"
        resource="SCRATCH"

    # Cori at NERSC (phase I)

    elif [ $MACHINE == "CORI" ]; then

	allocation="m1938"
	exec="sbatch"
	cancel_job="scancel"
	ppn="32"
        logical_ppn="64"
	run_ext=".OU"
	batch_system="SLURM"
	launcher="srun"
	queue="regular"
        constraint="haswell"
        resource="SCRATCH"

        job_prepend="export OMP_PROC_BIND=spread; export OMP_PLACES=threads"

    # LIRED at Stony Brook University

    elif [ $MACHINE == "LIRED" ]; then

	exec="qsub"
	cancel_job="qdel"
	pause_job="qhold"
	resume_job="qrls"
	ppn="24"
	threads_per_task="1"
	batch_system="PBS"
	launcher="mpirun"
	queue="medium"
	run_ext=".OU"
	job_prepend="module load shared; module load torque; module load maui; module load mvapich2; module load gcc"

    # SeaWulf at Stony Brook University

    elif [ $MACHINE == "SEAWULF" ]; then

	exec="qsub"
	cancel_job="qdel"
	pause_job="qhold"
	resume_job="qrls"
	ppn="28"
	threads_per_task="1"
	batch_system="PBS"
	launcher="mpirun"
	queue="medium"
	run_ext=".OU"
	job_prepend="module load shared; module load torque; module load maui; module load mvapich2; module load gcc"

    # Mira at ALCF

    elif [ $MACHINE == "MIRA" ]; then

	exec="qsub"
	cancel_job="qdel"
	pause_job="qhold"
	resume_job="qrls"
	ppn="16"
	queue="prod"
	run_ext=".OU"
	batch_system="COBALT"
	launcher="runjob"
	allocation="wdmerger"

    fi



    # Allow the user to request notification e-mails about jobs.
    # Requires the environment variable EMAIL_ADDRESS.

    if [ ! -z $EMAIL_ADDRESS ]; then

        if [ $exec == "qsub" ]; then
            opts_flag="-M $EMAIL_ADDRESS "

            mail_opts=""

            if [ ! -z $EMAIL_ON_ABORT ]; then
                mail_opts+="a"
            fi

            if [ ! -z $EMAIL_ON_START ]; then
                mail_opts+="b"
            fi

            if [ ! -z $EMAIL_ON_TERMINATE ]; then
                mail_opts+="e"
            fi

            if [ ! -z $EMAIL_ON_BAD_TERMINATE ]; then
                mail_opts+="f"
            fi

            if [ ! -z $NO_EMAIL ]; then
                mail_opts="p"
            fi

            if [ ! -z $mail_opts ]; then
                opts_flag+="-m $mail_opts "
            fi

            exec="qsub $opts_flag"

        fi

    fi
}
