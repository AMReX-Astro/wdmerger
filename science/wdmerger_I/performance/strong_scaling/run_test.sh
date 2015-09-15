source $WDMERGER_HOME/job_scripts/run_utils.sh

max_step=10

amr_plot_files_output=0
amr_checkpoint_files_output=0

# Run multiple problems that scale at different levels.

amr_max_level="1"
amr_ref_ratio="4"
amr_max_grid_size="32 48"

# For strong scaling, we fix the problem size and increase the number of processors.

for nprocs in 64 128 256 512 1024 2048
do
  dir=$results_dir/two_levels/n$nprocs
  walltime=4:00:00

  run $dir $nprocs $walltime
done

amr_max_level="2"
amr_ref_ratio="4 4"
amr_max_grid_size="32 48 64"

for nprocs in 2048 4096 8192 16384
do
  dir=$results_dir/three_levels/n$nprocs
  walltime=12:00:00

  run $dir $nprocs $walltime
done
