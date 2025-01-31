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

      UNAMEN=$(uname -n)$(hostname -f)

      if [[ $UNAMEN == *"frontier"* ]]; then
          MACHINE=FRONTIER
      elif [[ $UNAMEN == *"lassen"* ]]; then
          MACHINE=LASSEN
      elif [[ $UNAMEN == *"lired"*  ]]; then
          MACHINE=LIRED
      elif [[ $UNAMEN == "login"  ]]; then
          MACHINE=SEAWULF
      elif [[ $UNAMEN == *"mira"*   ]]; then
          MACHINE=MIRA
      fi

      if [ ! -z $NERSC_HOST ]; then
          if   [ $NERSC_HOST == "perlmutter" ]; then
              MACHINE=PERLMUTTER
          fi
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

    # Frontier at OLCF

    elif [ $MACHINE == "FRONTIER" ]; then

        allocation="ast106"
        exec="sbatch"
        cancel_job="scancel"
        ppn="8"
        threads_per_task="7"
        run_ext=".OU"
        batch_system="SLURM"
        queue="batch"
        launcher="srun"
        #TODO
        #archive_exec="sbatch"
        #archive_method="htar"
        #archive_queue="batch"
        #archive_wclimit="24:00:00"
        #archive_cluster="dtn"
        #do_storage_in_job=0

    # Lassen at LLNL

    elif [ $MACHINE == "LASSEN" ]; then

        allocation="lcstaff"
        exec="bsub"
        cancel_job="bkill"
        pause_job="bstop"
        resume_job="bresume"
        ppn="4"
        run_ext=".OU"
        batch_system="LSF"
        queue="batch"
        launcher="jsrun"

    # Perlmutter at NERSC

    elif [ $MACHINE == "PERLMUTTER" ]; then

        allocation="m3018_g"
        exec="sbatch"
        cancel_job="scancel"
        ppn="4"
        logical_ppn="4"
        run_ext=".OU"
        batch_system="SLURM"
        launcher="srun"
        qos="regular"
        constraint="gpu"
        gpus_per_task="1"
        gpu_bind="map_gpu:0,1,2,3"

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
