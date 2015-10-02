source $WDMERGER_HOME/job_scripts/run_utils.sh

max_step=10

amr_plot_files_output=0
amr_checkpoint_files_output=0

# Run multiple problems that scale at different levels.
# For strong scaling, we fix the problem size and increase the number of processors.

# Single level

amr_max_level="0"

for nprocs in 32 64 128 256
do
  dir=$results_dir/one_level/n$nprocs
  walltime=4:00:00

  run $dir $nprocs $walltime
done

# One level of refinement

amr_max_level="1"

for nprocs in 256 512 1024 2048
do
  dir=$results_dir/two_levels/n$nprocs
  walltime=4:00:00

  run $dir $nprocs $walltime
done

# Two levels of refinement

amr_max_level="2"

for nprocs in 2048 4096 8192 16384
do
  dir=$results_dir/three_levels/n$nprocs
  walltime=12:00:00

  run $dir $nprocs $walltime
done
