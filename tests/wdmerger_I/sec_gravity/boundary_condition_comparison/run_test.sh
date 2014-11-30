#!/bin/bash

# Some variable names

exec='mpiexec -n 8'
Castro='Castro3d.Linux.g++.gfortran.MPI.ex'
inputs='inputs'
probin='probin'

# Define a function that moves all of the output
# data from a Castro run to the directory in the first argument.

function move_results {

  if [ -d "$1" ]; then
    rm -rf $1/
  fi
  mkdir $1
  mv plt* $1/
  mv *.out $1/
  cp $inputs $1/
  cp $probin $1/

}

results_dir=results

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

for l in 0 2 6 20
do
  dir=$results_dir/$l
  if [ ! -d $dir ]; then
    mkdir $dir
    echo "Now doing l =" $l
    sed -i "/gravity.max_multipole_order/c gravity.max_multipole_order = $l" $inputs
    $exec $Castro $inputs > info.out
    move_results $dir
  fi
done

# Now do the 'exact' direct summation, for comparison purposes

dir=$results_dir/true

if [ ! -d $dir ]; then
  mkdir $dir
  echo "Now computing exact solution"
  sed -i "/gravity.direct_sum_bcs/c gravity.direct_sum_bcs = 1" $inputs
  $exec $Castro $inputs > info.out
  move_results $dir
  sed -i "/gravity.direct_sum_bcs/c gravity.direct_sum_bcs = 0" $inputs
fi

sed -i "/gravity.direct_sum_bcs/c gravity.direct_sum_bcs = 0" $inputs
