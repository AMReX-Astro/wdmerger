source $WDMERGER_HOME/job_scripts/run_utils.sh

# Loop over the resolutions in question

for ncell in 64 128 256
do
  dir=$results_dir/$ncell
  nprocs=16
  sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell" $inputs
  run $dir $nprocs
done
