source $WDMERGER_HOME/job_scripts/run_utils.sh

# Functions we'll need

function set_run_opts {

  if [ $MACHINE == "LIRED" ]; then
      if   [ $ncell == 128 ]; then
	  nprocs=24
	  walltime=12:00:00
      elif [ $ncell == 256 ]; then
	  nprocs=144
	  walltime=12:00:00
      elif [ $ncell == 512 ]; then
	  nprocs=288
	  walltime=12:00:00
      elif [ $ncell == 1024 ]; then
	  nprocs=576
	  walltime=12:00:00
      fi
  elif [ $MACHINE == "TITAN" ]; then
      if   [ $ncell == 128 ]; then
	  nprocs=16
	  walltime=2:00:00
      elif [ $ncell == 256 ]; then
	  nprocs=2048
	  walltime=6:00:00
      elif [ $ncell == 512 ]; then
	  nprocs=2048
	  walltime=6:00:00
      fi
  fi


  if [ ! -z $ncell ]; then

    if [ $ncell -eq 256 ]; then
	amr_n_cell="256 256 256"
	amr_max_level="0"
	amr_max_grid_size="32"
    elif [ $ncell -eq 512 ]; then
	amr_n_cell="256 256 256"
	amr_max_level="1"
	amr_ref_ratio="2"
	amr_max_grid_size="32 32"
    elif [ $ncell -eq 1024 ]; then
	amr_n_cell="256 256 256"
	amr_max_level="1"
	amr_ref_ratio="4"
	amr_max_grid_size="32 48"
    elif [ $ncell -eq 2048 ]; then
	amr_n_cell="256 256 256"
	amr_max_level="2"
	amr_ref_ratio="4 2"
	amr_max_grid_size="32 32 64"
    fi

  fi

}

amr_check_per=10.0

amr_small_plot_per=0.1
amr_small_plot_vars="density"

ncell_list="256 512 1024"

castro_do_react=0
castro_react_T_min=1.0e8

mass_P=0.90
mass_S=0.60

stop_time=10000.0

castro_riemann_solver=0

castro_sum_interval=1
castro_sum_per=-1.0

# Use piecewise-constant interpolation for regridding.

castro_state_interp_order=0





# First set up the problem using the approximate ICs.

problem=2
castro_do_sponge=0
update_stopping_criteria=1

for mass_P in 0.90
do
  for mass_S in 0.60 0.90
  do
    for ncell in $ncell_list
    do
      for roche_radius_factor in 0.90 1.00
      do
        for castro_do_rotation in 0 1
	do
          for castro_hybrid_hydro in 0 1
	  do

	      if [ $mass_S == "0.90" ]; then

		  if [ $roche_radius_factor != "0.90" ]; then
		      continue
		  fi

		  if [ $ncell != "512" ]; then
		      continue
		  fi

		  if [ $castro_hybrid_hydro != "1" ]; then
		      continue
		  fi

		  stop_time=100.0
		  amr_plot_per=5.0

	      elif [ $mass_S == "0.60" ]; then

		  if [ $roche_radius_factor == "1.00" ]; then

		      if [ $ncell == "256" ]; then

			  stop_time=200.0

		      elif [ $ncell == "512" ]; then

			  if [ $castro_hybrid_hydro == "1" ]; then
			      stop_time=200.0
			  else
			      stop_time=100.0
			  fi

		      elif [ $ncell == "1024" ]; then

			  stop_time=50.0

		      else

			  echo "Unknown value for ncell"
			  exit

		      fi

		      amr_plot_per=1.0

		  elif [ $roche_radius_factor == "0.90" ]; then

		      if [ $ncell != "512" ]; then
			  continue
		      fi

		      if [ $castro_hybrid_hydro != "1" ]; then
			  continue
		      fi

		      stop_time=175.0
		      amr_plot_per=5.0

		  fi

	      else

		  continue

	      fi

	      dir=$results_dir/approximate/mass_P_$mass_P/mass_S_$mass_S/roche$roche_radius_factor/rot$castro_do_rotation/hybrid$castro_hybrid_hydro/n$ncell

	      set_run_opts
	      run

	  done
        done
      done
    done
  done
done



# Now set up the problem using the accurate ICs.

problem=3
castro_do_rotation=1
castro_do_sponge=1

ncell_list="256"

roche_radius_factor=3.0

amr_max_level=0

# First we want to study equilibrium systems alone.
# So make it impossible to disable the relaxation,
# use a negative radial damping factor, and figure 
# out what the right damping timescale is.

relaxation_density_cutoff="1.0d10"
radial_damping_factor="-1.0d0"

timescale_list="50.0d0" #"10.0d0" # 20.0d0"
roche_list="3.0d0"

for roche_radius_factor in $roche_list
do
    for relaxation_damping_timescale in $timescale_list
    do
	for castro_hybrid_hydro in 0 1
	do
	    for castro_state_in_rotating_frame in 1 #0 1
	    do
		for ncell in $ncell_list
		do

		    dir=$results_dir/equilibrium/r$roche_radius_factor/t$relaxation_damping_timescale/h$castro_hybrid_hydro/f$castro_state_in_rotating_frame/n$ncell/

		    set_run_opts
		    run

		done
	    done
	done
    done
done
