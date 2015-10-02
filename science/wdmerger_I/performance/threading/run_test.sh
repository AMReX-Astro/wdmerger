source $WDMERGER_HOME/job_scripts/run_utils.sh

# Problem-specific variables

mass_P=0.90
mass_S=0.90

castro_rotational_period=100.0

max_step=10

amr_plot_files_output=0
amr_checkpoint_files_output=0

# Loop over resolutions; the optimal 
# threading will depend on problem size.

for ncell in 256 1024 4096
do

    # Loop over the number of threads

    for OMP_NUM_THREADS in 1 2 4 8
    do

	dir=$results_dir/n$ncell/omp$OMP_NUM_THREADS/

	amr_n_cell="256 256 256"

	if [ $ncell -eq 1024 ]; then
	    amr_max_level="1"
	    amr_ref_ratio="4"
	    amr_max_grid_size="32 32"
	elif [ $ncell -eq 4096 ]; then
	    amr_max_level="2"
	    amr_ref_ratio="4 4"
	    amr_max_grid_size="32 48 64"
	fi

	if [ $MACHINE == "BLUE_WATERS" ]; then
	    if   [ $ncell -eq 256 ]; then
		nprocs=256
		walltime=1:00:00
	    elif [ $ncell -eq 1024 ]; then
		nprocs=2048
		walltime=1:00:00
	    elif [ $ncell -eq 4096 ]; then
		nprocs=16384
		walltime=6:00:00
	    fi
	fi

	run $dir $nprocs $walltime

    done

done
