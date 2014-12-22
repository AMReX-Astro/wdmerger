#!/bin/bash

# Some variable names

exec='mpiexec -n 8'
Castro='Castro3d.Linux.g++.gfortran.MPI.ex'
inputs='inputs'
probin='probin'

function copy_files {

    cp $Castro $1
    cp $inputs $1
    cp $probin $1

}

function run {

    echo "$exec $Castro $inputs > info.out" | batch

}

# Check if results directory already exists, and if not then create it.

results_dir=results

if [ ! -d $results_dir ]; then
  mkdir $results_dir
fi

# Loop over problem choices

for problem in 1 2
do

  dir=$results_dir/problem$problem/
  if [ ! -d $dir ]; then
      echo "Submitting problem =" $problem
      mkdir $dir
      sed -i "/problem/c problem = $problem" $probin
  fi

  # Loop over possible bulk flow velocities

  for vel in 0 1 3 10 30 100
  do
      dir=$results_dir/problem$problem/velocity$vel
      if [ ! -d $dir ]; then
	  echo "$ubmitting velocity =" $vel
	  mkdir $dir
	  sed -i "/bulk_velocity/c bulk_velocity = $vel" $probin
      fi

      # Loop over the resolutions in question

      for ncell in 64
      do
	dir=$results_dir/problem$problem/velocity$vel/$ncell
	if [ ! -d $dir ]; then
	  echo "Submitting ncell =" $ncell
	  mkdir $dir
	  ncell_x=$ncell
	  ncell_y=$ncell
	  ncell_z=8

	  # Since CASTRO only supports dx = dy = dz,
	  # we must change the size of the problem in the z
	  # direction to ensure that dz stays equal to dx and dy.

	  problo_x=0.0
	  probhi_x=1.0
	  problo_y=0.0
	  probhi_y=1.0
	  problo_z=0.0
	  probhi_z=$(echo "$probhi_x / $ncell_x * $ncell_z" | bc -l)

	  sed -i "/geometry.prob_lo/c geometry.prob_lo = $problo_x $problo_y $problo_z" $inputs
	  sed -i "/geometry.prob_hi/c geometry.prob_hi = $probhi_x $probhi_y $probhi_z" $inputs
	  sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell_z" $inputs
	  copy_files $dir
	  cd $dir
	  run
	  cd -
	fi
      done
  done
done
