source $WDMERGER_HOME/job_scripts/run_utils.sh

# Loop over the resolutions in question

for ncell in 32 64
do
  dir=$results_dir/$ncell
  sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell" $inputs
  run $dir
done
