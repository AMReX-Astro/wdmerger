SHELL = /bin/bash

# Pull in useful directories and filenames, and include any others

include var_names
ca_make = castro_makefile

# Main make target

all:

	# Copy model files from main wdmerger model_files location,
	# and insert their names into the probin file.
	cp -f $(model_dir)/$(model) .
	sed -i "/model_name/c model_name = \"$(model)\"" $(probin)


	# Copy CASTRO makefile from main directory
	if [ ! -e $(ca_make) ]; then cp $(main_dir)/GNUmakefile $(ca_make); fi

	# Make options: turn off OpenMP, use gcc, and turn off hydro and rotation.
	sed -i "/USE_MPI/c USE_MPI = TRUE" $(ca_make)
	sed -i "/USE_OMP/c USE_OMP = FALSE" $(ca_make)
	sed -i "/USE_GRAV/c USE_GRAV = TRUE" $(ca_make)
	sed -i "/USE_ROTATION/c USE_ROTATION = TRUE" $(ca_make)

	make -f $(ca_make)

# Clean everything up with make clean

clean:
	rm -f $(model)
	rm -f helm_table.dat
	if [ -e  "$(ca_make)" ]; then \
		make realclean -f $(ca_make); \
	fi;
	rm -rf plt*
	rm -rf chk*
	rm -rf results
	rm -f *.out
	rm -f *.png
	git checkout $(probin)
	git checkout $(inputs)
