source $WDMERGER_HOME/job_scripts/run_utils.sh

# Make sure the necessary gravity BC options are actually in the inputs file,
# since we'll be modifying them as we go.

if [ $(grep -F "gravity.direct_sum_bcs" $inputs | wc -l) -lt 1 ]; then
  echo "" >> $inputs
  echo "gravity.direct_sum_bcs = 0" >> $inputs
else
  sed -i "/gravity.direct_sum_bcs/c gravity.direct_sum_bcs = 0" $inputs
fi

if [ $(grep -F "gravity.max_multipole_order" $inputs | wc -l) -lt 1 ]; then
  echo "" >> $inputs
  echo "gravity.max_multipole_order = 0" >> $inputs
fi

# Loop over the multipole orders we want to examine

for l in {0..20}
do
  dir=$results_dir/$l
  sed -i "/gravity.max_multipole_order/c gravity.max_multipole_order = $l" $inputs
  run $dir
done

# Now do the 'exact' direct summation, for comparison purposes

dir=$results_dir/true
sed -i "/gravity.direct_sum_bcs/c gravity.direct_sum_bcs = 1" $inputs
run $dir

