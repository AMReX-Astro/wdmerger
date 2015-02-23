source $WDMERGER_HOME/job_scripts/run_utils.sh

# Problem-specific variables

max_step=10

amr_plot_files_output=0
amr_checkpoint_files_output=0

# For strong scaling, we fix the problem size and increase the number of processors.

for nprocs in 2048 4096 8192 16384
do
  dir=$results_dir/$nprocs
  run $dir $nprocs
done
