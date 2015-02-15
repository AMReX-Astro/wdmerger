source $WDMERGER_HOME/job_scripts/run_utils.sh

if [ $MACHINE == "BLUE_WATERS" ]; then
    nprocs=16
    walltime=1:00:00
fi

# Loop over the multipole orders we want to examine

for l in {0..20}
do
  dir=$results_dir/$l
  sed -i "/gravity.max_multipole_order/c gravity.max_multipole_order = $l" $compile_dir/$inputs
  run $dir $nprocs $walltime
done

# Now do the 'exact' direct summation, for comparison purposes

dir=$results_dir/true
sed -i "/gravity.direct_sum_bcs/c gravity.direct_sum_bcs = 1" $compile_dir/$inputs
run $dir $nprocs $walltime

