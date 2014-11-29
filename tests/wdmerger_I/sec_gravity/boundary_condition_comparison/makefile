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

	# Copy necessary setup files from the main directory,
	# and then build the Castro executable.
	if [ ! -e Prob_3d.f90 ]; then cp $(main_dir)/Prob_3d.f90 .; fi
	cp $(main_dir)/probdata.f90 .
	cp $(main_dir)/com.f90 .
	cp $(main_dir)/Make.package .
	cp $(main_dir)/GNUmakefile $(ca_make)

	# Make options: turn off OpenMP, use gcc, and turn off hydro and rotation.
	sed -i "/USE_MPI/c USE_MPI = TRUE" $(ca_make)
	sed -i "/USE_OMP/c USE_OMP = FALSE" $(ca_make)
	sed -i "/USE_HYDRO/c USE_HYDRO = FALSE" $(ca_make)
	sed -i "/USE_GRAV/c USE_GRAV = TRUE" $(ca_make)
	sed -i "/USE_ROTATION/c USE_ROTATION = TRUE" $(ca_make)

	make -f $(ca_make)

# Clean everything up with make clean

clean:
	rm -f $(model_P) $(model_S) 
	rm -f Prob_3d.f90 probdata.f90 com.f90 
	rm -f helm_table.dat
	if [ -e  "$(ca_make)" ]; then \
		make realclean -f $(ca_make); \
		rm -f Make.package $(ca_make); \
	fi;
	git checkout $(probin)
	git checkout $(inputs)
	rm -rf plt*
	rm -rf chk*
	rm -rf results
	rm -f *.out

