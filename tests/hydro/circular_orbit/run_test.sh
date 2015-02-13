source $WDMERGER_HOME/job_scripts/run_utils.sh

# Loop over the resolutions in question

for ncell in 128
do
  dir=$results_dir/$ncell
  sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell" $compile_dir/$inputs
  run $dir
done
