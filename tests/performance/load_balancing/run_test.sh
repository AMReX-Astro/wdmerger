source $WDMERGER_HOME/job_scripts/run_utils.sh

# Set problem options and run time

max_step=10

amr_plot_files_output=0
amr_checkpoint_files_output=0

amr_max_level="1"
amr_ref_ratio="4"
amr_max_grid_size="32 48"

nprocs="512"

walltime=1:00:00

# Test various possibilities for the distribution mapping strategy.

strategy_list="ROUNDROBIN KNAPSACK SFC PFC RRSFC"

# Run this test multiple times to minimize risks from inter-run variability.

n_iters=3
iter_list=$(seq $n_iters)

for strat in $strategy_list
do

  for r in $iter_list
  do

    dir=$results_dir/strat_$strat/rep$r

    DistributionMapping_strategy=$strat

    run $dir $nprocs $walltime

  done
done
