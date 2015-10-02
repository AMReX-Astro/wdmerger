source $WDMERGER_HOME/job_scripts/run_utils.sh

castro_do_rotation=0

collision=T
collision_separation=4.0

amr_plot_per=0.1
amr_check_per=0.1

stop_time=10.0

castro_do_react=1

mass_P_list=(0.50 0.55 0.60 0.64 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00)
mass_S_list=(0.50 0.55 0.60 0.64 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00)

list_length=${#mass_P_list[@]}

for index in $(seq 0 $(($list_length-1)))
do
    
  mass_P=${mass_P_list[$index]}
  mass_S=${mass_S_list[$index]}

  for ncell in 256 512 1024
  do
      dir=$results_dir/m_P_$mass_P/m_S_$mass_S/n$ncell

      if [ $ncell -eq 256 ]; then
	  amr_n_cell="$ncell $ncell $ncell"
      elif [ $ncell -eq 512 ]; then
	  amr_n_cell="256 256 256"
	  amr_max_level="1"
	  amr_ref_ratio="2"
      elif [ $ncell -eq 1024 ]; then
	  amr_n_cell="256 256 256"
	  amr_max_level="1"
	  amr_ref_ratio="4"
      fi

      if [ $MACHINE == "BLUE_WATERS" ]; then

	  if [ $ncell -eq 256 ]; then
	      nprocs="2048"
	      walltime="12:00:00"
	  elif [ $ncell -eq 512 ]; then
	      nprocs="2048"
	      walltime="24:00:00"
	  elif [ $ncell -eq 1024 ]; then
	      nprocs="4096"
	      walltime="24:00:00"
	  fi
	  
      fi

      run
  done
done
