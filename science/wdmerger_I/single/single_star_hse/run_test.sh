source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

  if [ $MACHINE == "BLUE_WATERS" ]; then

      if   [ $ncell -eq 256 ]; then
	  nprocs=512
	  walltime=24:00:00
      elif [ $ncell -eq 512 ]; then
	  nprocs=512
	  walltime=24:00:00
      elif [ $ncell -eq 1024 ]; then
	  nprocs=512
	  walltime=24:00:00
      elif [ $ncell -eq 2048 ]; then
	  nprocs=512
	  walltime=24:00:00
      elif [ $ncell -eq 4096 ]; then
	  nprocs=512
	  walltime=24:00:00
      fi

  elif [ $MACHINE == "LIRED" ]; then

      if   [ $ncell -eq 256 ]; then
	  nprocs=144
	  walltime=12:00:00
      elif [ $ncell -eq 512 ]; then
	  nprocs=288
	  walltime=12:00:00
      elif [ $ncell -eq 1024 ]; then
	  nprocs=576
	  walltime=12:00:00
      elif [ $ncell -eq 2048 ]; then
	  nprocs=576
	  walltime=12:00:00
      elif [ $ncell -eq 4096 ]; then
	  nprocs=576
	  walltime=12:00:00
      fi

  fi

  if [ -z $ncell ]; then

      if   [ $ncell -eq 512 ]; then
	  amr_n_cell="256 256 256"
	  amr_ref_ratio="2"
	  amr_max_level="1"
      elif [ $ncell -eq 1024 ]; then
	  amr_n_cell="256 256 256"
	  amr_ref_ratio="4"
	  amr_max_level="1"
	  amr_max_grid_size="32 32"
      elif [ $ncell -eq 2048 ]; then
	  amr_n_cell="256 256 256"
	  amr_ref_ratio="4 2"
	  amr_max_level="2"
	  amr_max_grid_size="32 32 32"
      elif [ $ncell -eq 4096 ]; then
	  amr_n_cell="256 256 256"
	  amr_ref_ratio="4 4"
	  amr_max_level="2"
	  amr_max_grid_size="32 32 48"
      else
	  amr_n_cell="$ncell $ncell $ncell"
      fi

  fi



}



# Problem-specific variables

mass_P=" 0.90"
mass_S="-1.00"

stop_time=200.0
castro_do_react=0
castro_do_rotation=0

amr_check_per=20.0
amr_plot_per=20.0

# Loop over the resolutions for the stationary case

for ncell in 256 512 1024
do

  dir=$results_dir/static/n$ncell

  set_run_opts
  run

done




# Now we'll do the case where the star moves.
# We want the star to evolve for the same amount
# of time as the previous simulations. One way to
# do this is to start the star in the bottom left
# corner and have it move toward the top right.
# We widen the grid and use the AMR to zoom in
# so that it has the same effective resolution
# as the 256**3 run.

geometry_prob_lo="-2.048e10 -2.048e10 -2.048e10"
geometry_prob_hi=" 2.048e10  2.048e10  2.048e10"

center_fracx=0.125d0
center_fracy=0.125d0
center_fracz=0.125d0

bulk_velocity=256000000

bulk_velx=$(echo "$bulk_velocity / sqrt(3)" | bc -l)
bulk_vely=$(echo "$bulk_velocity / sqrt(3)" | bc -l)
bulk_velz=$(echo "$bulk_velocity / sqrt(3)" | bc -l)

for ncell in 1024 2048 4096
do

  dir=$results_dir/motion/n$ncell

  set_run_opts
  run

done
