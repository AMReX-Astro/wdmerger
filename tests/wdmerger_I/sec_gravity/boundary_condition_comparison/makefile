# Define some useful directories and filenames

main_dir  = ../../../..
model_dir = $(main_dir)/model_files
model_P   = sub_chandra.M_WD-0.90.hse.CO.100.00km.nocoulomb.hot
model_S   = sub_chandra.M_WD-0.81.hse.CO.100.00km.nocoulomb.hot

inputs  = inputs_3d
probin  = probin
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

clean:
	rm -f $(model_P) $(model_S) 
	rm -f Prob_3d.f90 probdata.f90 com.f90 
	rm -f helm_table.dat
	if [ -e  "$(ca_make)" ]; then \
		make realclean -f $(ca_make); \
		rm -f Make.package $(ca_make); \
	fi;

