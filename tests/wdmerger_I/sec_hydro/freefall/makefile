SHELL = /bin/bash

# Pull in useful directories and filenames, and include any others

include var_names
ca_make = castro_makefile

# Main make target

all:

	# Copy model files from main wdmerger model_files location,
	# and insert their names into the probin file.
	cp -f $(model_dir)/$(model_P) .
	cp -f $(model_dir)/$(model_S) .
	sed -i "/model_P_name/c model_P_name = \"$(model_P)\"" $(probin)
	sed -i "/model_S_name/c model_S_name = \"$(model_S)\"" $(probin)

	# Insert the orbital period and stopping time into the inputs file
	sed -i "/castro.rotational_period/c castro.rotational_period = $(rotational_period)" $(inputs)
	sed -i "/stop_time/c stop_time = $(stop_time)" $(inputs)

	# Copy necessary setup files from the main directory,
	# and then build the Castro executable.
	for file in Prob_3d.f90 probdata.f90 com.f90 Make.package sum_integrated_quantities.cpp Problem.H Problem_F.H ; \
	do \
		if [ ! -e $$file ]; then cp $(main_dir)/$$file .; fi \
	done

	if [ ! -e $(ca_make) ]; then cp $(main_dir)/GNUmakefile $(ca_make); fi

	# Make options: turn off OpenMP, use gcc, and turn off hydro and rotation.
	sed -i "/USE_MPI/c USE_MPI = TRUE" $(ca_make)
	sed -i "/USE_OMP/c USE_OMP = FALSE" $(ca_make)
	sed -i "/USE_GRAV/c USE_GRAV = TRUE" $(ca_make)
	sed -i "/USE_ROTATION/c USE_ROTATION = TRUE" $(ca_make)

	make -f $(ca_make)

# Clean everything up with make clean

clean:
	rm -f $(model_P) $(model_S) 
	rm -f Prob_3d.f90 probdata.f90 com.f90 
	rm -f Problem.H Problem_F.H
	rm -f sum_integrated_quantities.cpp
	rm -f helm_table.dat
	if [ -e  "$(ca_make)" ]; then \
		make realclean -f $(ca_make); \
		rm -f Make.package $(ca_make); \
	fi;
	rm -rf plt*
	rm -rf chk*
	rm -rf results
	rm -f *.out
	rm -f *.png
	git checkout $(probin)
	git checkout $(inputs)
