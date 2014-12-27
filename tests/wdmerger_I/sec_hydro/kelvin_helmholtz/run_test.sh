source $WDMERGER_HOME/job_scripts/run_utils.sh

# Loop over problem choices

for problem in 1 2
do

  # Loop over possible bulk flow velocities

  for vel in 0 1 3 10 30 100
  do

    # Loop over the resolutions in question

    for ncell in 64
    do
      dir=$results_dir/problem$problem/velocity$vel/$ncell
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

      sed -i "/problem/c problem = $problem" $probin
      sed -i "/geometry.prob_lo/c geometry.prob_lo = $problo_x $problo_y $problo_z" $inputs
      sed -i "/geometry.prob_hi/c geometry.prob_hi = $probhi_x $probhi_y $probhi_z" $inputs
      sed -i "/amr.n_cell/c amr.n_cell = $ncell $ncell $ncell_z" $inputs

      run $dir

    done
  done
done
