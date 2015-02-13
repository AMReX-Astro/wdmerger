source $WDMERGER_HOME/job_scripts/run_utils.sh

# Loop over the resolutions in question

for gst in 1 2 3 4
do
  for ncell in 64 128 256
  do
    dir=$results_dir/gs$gst/ncell$ncell
    nprocs=16
    sed -i "/castro.grav_source_type/c castro.grav_source_type = $gst" $compile_dir/$inputs
    sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell" $compile_dir/$inputs
    run $dir $nprocs
  done
done
