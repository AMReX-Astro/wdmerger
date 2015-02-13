source $WDMERGER_HOME/job_scripts/run_utils.sh

# Loop over the resolutions in question

for problem in 1 2 3
do

  for ncell in 16 32 64 128 256
  do

    dir=$results_dir/problem$problem/$ncell
    sed -i "/problem/c problem = $problem" $compile_dir/$probin
    sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell" $compile_dir/$inputs
    run $dir

  done
done
