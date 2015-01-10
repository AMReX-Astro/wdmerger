SHELL=/bin/bash

merger_makefile=$(WDMERGER_HOME)/source/GNUmakefile
local_makefile=castro_makefile
inputs=inputs
probin=probin
run_utils=$(WDMERGER_HOME)/job_scripts/run_utils.sh

# Prints a variable VAR when makefile is called as
# make print-VAR

print-%: ; @echo $* is $($*)

# Makefile options. We'll use the ?= method so that
# the including makefile can override these for custom tests
# (e.g. changing the number of dimensions or the EOS).

USE_MPI ?= TRUE
USE_OMP ?= FALSE
COMP ?= g++
FCOMP ?= gfortran
DIM ?= 3
EOS_dir ?= helmeos
USE_GRAV ?= TRUE

# Sets up the CASTRO makefile in the local directory
# and inserts our preferred make options like compiler
# and whether to use MPI and/or OpenMP. Should be called
# by the main make target in any run directory.

make_castro:

	# Check if makefile exists; if not, copy it from the source directory,
	# then edit the makefile variables that we want to use, including
	# compiler options. If we want to use a custom setup for a particular
	# problem, instead of using the main source files, check for a file named
	# full_problem in the source directory.
	# The < and > characters in sed mark off word boundaries, and are useful
	# so that COMP does not overwrite FCOMP, which also matches.

	if [ ! -e $(local_makefile) ] ; \
	then \
		cp $(merger_makefile) $(local_makefile) ; \
		sed -i "/\<COMP\>/c COMP = $(COMP)" $(local_makefile) ; \
		sed -i "/\<FCOMP\>/c FCOMP = $(FCOMP)" $(local_makefile) ; \
		sed -i "/DIM/c DIM = $(DIM)" $(local_makefile) ; \
		sed -i "/EOS_dir/c EOS_dir := $(EOS_dir)" $(local_makefile) ; \
		sed -i "/USE_GRAV/c USE_GRAV = $(USE_GRAV)" $(local_makefile) ; \
		sed -i "/USE_MPI/c USE_MPI = $(USE_MPI)" $(local_makefile) ; \
		sed -i "/USE_OMP/c USE_OMP = $(USE_OMP)" $(local_makefile) ; \
		if [ -e "source/full_source" ] ; \
		then \
			sed -i "/Blocs/c Blocs := source/" $(local_makefile) ; \
			sed -i "/Bpack/c Bpack := source/Make.package" $(local_makefile) ; \
		elif [ -e "source/partial_source" ] ; \
		then \
			sed -i "/Blocs/c Blocs := source/ $(WDMERGER_HOME)/source/" $(local_makefile) ; \
		fi \
	fi

	make -f $(local_makefile)

# Clean everything up with make clean

clean_targets = *.mod *.run plt* chk* *.out *.png helm_table.dat results

clean:
	if [ -e $(local_makefile) ]; then \
		make realclean -f $(local_makefile); \
		rm -f $(local_makefile); \
	fi
	rm -rf $(clean_targets)
	git checkout $(probin)
	git checkout $(inputs)
