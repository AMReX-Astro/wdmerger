SHELL = /bin/bash

# Pull in useful directories and filenames, and include any others

include var_names

# Main make target

all:

	cp $(WDMERGER_HOME)/source/$(merger_make) $(makefile)

	# Make options: turn off OpenMP, use gcc, and turn off hydro and rotation.
	sed -i "/USE_MPI/c USE_MPI = TRUE" $(makefile)
	sed -i "/USE_OMP/c USE_OMP = FALSE" $(makefile)
	sed -i "/USE_HYDRO/c USE_HYDRO = FALSE" $(makefile)
	sed -i "/USE_GRAV/c USE_GRAV = FALSE" $(makefile)
	sed -i "/EOS_dir/c EOS_dir := gamma_law_general" $(makefile)
	sed -i "/USE_ROTATION/c USE_ROTATION = TRUE" $(makefile)

	# Override the source includes since this is a custom setup.
	sed -i "/Blocs/c Blocs := ." $(makefile)
	sed -i "/Bpack/c Bpack := ./Make.package" $(makefile)

	make -f $(makefile)

# Clean everything up with make clean

clean:
	make realclean -f $(makefile)
	rm -f $(makefile)
	rm -f helm_table.dat
	rm -rf plt*
	rm -rf chk*
	rm -rf results
	rm -f *.out
	rm -f *.png
	rm -f *.eps
	git checkout $(probin)
	git checkout $(inputs)
