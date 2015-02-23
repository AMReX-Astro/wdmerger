source $WDMERGER_HOME/job_scripts/run_utils.sh

# Problem-specific variables

max_step=10

amr_plot_files_output=0
amr_checkpoint_files_output=0

# For weak scaling, we increase the problem size in parallel with the number of processors.

for ncell in 64 128 256 512
do
  dir=$results_dir/$nprocs
  amr_n_cell="$ncell $ncell $ncell"

  if [ $ncell -eq "64" ]; then
      nprocs="16"
  elif [ $nprocs -eq "128" ]; then
      nprocs="128"
  elif [ $nprocs -eq "256" ]; then
      nprocs="1024"
  elif [ $nprocs -eq "512" ]; then
      nprocs="8192"
  fi

  run $dir $nprocs
done
