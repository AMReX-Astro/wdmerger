source $WDMERGER_HOME/job_scripts/run_utils.sh

max_step=10

amr_plot_files_output=0
amr_checkpoint_files_output=0

# For weak scaling, we increase the problem size in parallel with the number of processors.

for ncell in 128 256 512 1024 2048
do
  dir=$results_dir/n$ncell
  amr_n_cell="$ncell $ncell $ncell"

  if   [ $ncell -eq 128 ]; then
      amr_max_grid_size="64"
      nprocs=8
  elif [ $ncell -eq 256 ]; then
      amr_max_grid_size="64"
      nprocs=64
  elif [ $ncell -eq 512 ]; then
      amr_max_grid_size="64"
      nprocs=512
  elif [ $ncell -eq 1024 ]; then
      amr_max_grid_size="64"
      nprocs=4096
  elif [ $ncell -eq 2048 ]; then
      amr_max_grid_size="64"
      nprocs=32768
  fi

  walltime=1:00:00

  run $dir $nprocs $walltime
done
