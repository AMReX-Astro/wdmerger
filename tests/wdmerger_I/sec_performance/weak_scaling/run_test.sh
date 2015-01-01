source $WDMERGER_HOME/job_scripts/run_utils.sh

# For weak scaling, we increase the problem size in parallel with the number of processors.

for ncell in 64 128 256 512
do
  dir=$results_dir/$nprocs
  sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell" $inputs
  if [ $ncell -eq "64" ]; then
      nprocs="16"
  elif [ $nprocs -eq "128" ]; then
      nprocs="128"
  elif [ $nprocs -eq "256" ]; then
      nprocs="1024"
  elif [ $nprocs -eq "512" ]; then
      nprocs="8192"
  fi
  run $dir $nprocs
done
