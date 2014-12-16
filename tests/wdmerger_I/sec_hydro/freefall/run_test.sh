#!/bin/bash

# Some variable names

exec='mpiexec -n 8'
Castro='Castro3d.Linux.g++.gfortran.MPI.ex'
inputs='inputs'
probin='probin'

function copy_files {

    cp $Castro $1
    cp helm_table.dat $1    
    cp $inputs $1
    cp $probin $1
    cp sub* $1

}

function run {

    echo "$exec $Castro $inputs > info.out" | batch

}

# Check if results directory already exists, and if not then create it.

results_dir=results

if [ ! -d $results_dir ]; then
  mkdir $results_dir
fi

# Loop over the resolutions in question

for ncell in 32 64
do
  dir=$results_dir/$ncell
  if [ ! -d $dir ]; then
    echo "Now doing ncell =" $ncell
    mkdir $dir
    sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell" $inputs
    copy_files $dir
    cd $dir
    run
    cd -
  fi
done
