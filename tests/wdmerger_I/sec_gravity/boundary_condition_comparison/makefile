SHELL = /bin/bash

# Pull in useful directories and filenames, and include any others

include var_names

# Main make target

all:
	cp $(merger_makefile) $(makefile)

	# Insert masses into probin file.
	sed -i "/mass_P/c mass_P = $(mass_P)" $(probin)
	sed -i "/mass_S/c mass_S = $(mass_S)" $(probin)

	# Make options: turn off OpenMP, use gcc, and turn off hydro and rotation.
	sed -i "/USE_MPI/c USE_MPI = TRUE" $(makefile)
	sed -i "/USE_OMP/c USE_OMP = FALSE" $(makefile)
	sed -i "/USE_GRAV/c USE_GRAV = TRUE" $(makefile)
	sed -i "/USE_ROTATION/c USE_ROTATION = TRUE" $(makefile)

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
	git checkout $(probin)
	git checkout $(inputs)

