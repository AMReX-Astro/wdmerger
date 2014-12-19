#!/bin/bash

# Some variable names

exec='mpiexec -n 8'
Castro='Castro3d.Linux.g++.gfortran.MPI.ex'
inputs='inputs'
probin='probin'

results_dir=results

function copy_files {

    cp $Castro $1
    cp helm_table.dat $1    
    cp $inputs $1
    cp $probin $1

}

function run {

    echo "$exec $Castro $inputs > info.out" | batch

}

if [ ! -d $results_dir ]; then
  mkdir $results_dir
fi

# Make sure the necessary gravity BC options are actually in the inputs file,
# since we'll be modifying them as we go.

if [ $(grep -F "gravity.direct_sum_bcs" $inputs | wc -l) -lt 1 ]; then
  echo "" >> $inputs
  echo "gravity.direct_sum_bcs = 0" >> $inputs
fi

if [ $(grep -F "gravity.max_multipole_order" $inputs | wc -l) -lt 1 ]; then
  echo "" >> $inputs
  echo "gravity.max_multipole_order = 0" >> $inputs
fi

# Loop over the multipole orders we want to examine

for l in {0..20}
do
  dir=$results_dir/$l
  if [ ! -d $dir ]; then
    echo "Submitting l =" $l
    mkdir $dir
    sed -i "/gravity.max_multipole_order/c gravity.max_multipole_order = $l" $inputs
    copy_files $dir
    cd $dir
    run
    cd -
  fi
done

# Now do the 'exact' direct summation, for comparison purposes

dir=$results_dir/true

if [ ! -d $dir ]; then
  echo "Submitting exact solution"
  mkdir $dir
  sed -i "/gravity.direct_sum_bcs/c gravity.direct_sum_bcs = 1" $inputs
  copy_files $dir
  sed -i "/gravity.direct_sum_bcs/c gravity.direct_sum_bcs = 0" $inputs
  cd $dir
  run
  cd -
fi
