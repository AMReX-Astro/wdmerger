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
	sed -i "/model_P_name/c model_P_name = $(model_P)" $(probin)
	sed -i "/model_S_name/c model_S_name = $(model_S)" $(probin)

	# Copy necessary setup files from the main directory,
	# and then build the Castro executable.
	cp $(main_dir)/Prob_3d.f90 .
	cp $(main_dir)/probdata.f90 .
	cp $(main_dir)/com.f90 .
	cp $(main_dir)/Make.package .
	cp $(main_dir)/GNUmakefile castro_makefile

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
