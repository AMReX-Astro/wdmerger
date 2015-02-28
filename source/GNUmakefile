SHELL=/bin/bash

# Makefile options. We'll use the ?= method so that
# any including makefile can override these for custom tests
# (e.g. changing the number of dimensions or the EOS).

PRECISION ?= DOUBLE
PROFILE ?= FALSE

DEBUG ?= FALSE

USE_MPI ?= TRUE
USE_OMP ?= FALSE

DIM ?= 3

USE_GRAV ?= TRUE
USE_REACT ?= FALSE
USE_MODELPARSER ?= TRUE
USE_ROTATION ?= TRUE

# Override default compiler for certain machines.

UNAMEN := $(shell uname -n)

ifeq ($(findstring titan, $(UNAMEN)), titan)
	# OLCF Titan
	COMP ?= Cray
	FCOMP ?= Cray
else ifeq ($(findstring h2o, $(UNAMEN)), h2o)
	# NCSA Blue Waters
	COMP ?= Cray
	FCOMP ?= Cray
else
	COMP ?= g++
	FCOMP ?= gfortran
endif

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     ?= helmeos

# This sets the network directory in $(NETWORK_HOME)
Network_dir ?= general_null
GENERAL_NET_INPUTS := $(CASTRO_DIR)/Networks/general_null/triple_alpha_plus_o.net 

# Determine where our source files are coming from
# based on the relevant file in the source directory.
# For building the full program, this will be irrelevant
# and we will simply compile everything in the 
# WDMERGER_HOME source directory.

diag_locs = $(WDMERGER_HOME)/source/diagnostics/
diag_pack = $(WDMERGER_HOME)/source/diagnostics/Make.package

root_locs = $(WDMERGER_HOME)/source/
root_pack = $(WDMERGER_HOME)/source/Make.package

prob_locs = ../source/
prob_pack = ../source/Make.package

ifneq ($(wildcard ../source/full_source),)
	Blocs := $(prob_locs) $(diag_locs)
	Bpack := $(prob_pack) $(diag_pack)
else
	Blocs := $(root_locs) $(diag_locs)
	Bpack := $(root_pack) $(diag_pack)
endif

BUILD_GIT_NAME := wdmerger
BUILD_GIT_DIR  := $(WDMERGER_HOME)

# This is a hack that we need to tell BoxLib
# not to build the dependency list before we execute the rule
# if we're doing a cleanall.
# MAKECMDGOALS is an automatically filled variable created by make
# that stores the name of the rule we gave to make
# at the command line. We probably shouldn't be messing 
# with it, but this works.

ifeq ($(MAKECMDGOALS),cleanall)
	MAKECMDGOALS = realclean
endif

include $(CASTRO_DIR)/Exec/Make.Castro

# Clean executable up with make clean;
# clean up all build files including objects with make realclean;
# clean up everything, including all results and analysis, with make cleanall.

clean_targets = *.mod helm_table.dat $(inputs) $(probin) 
cleanall_targets = $(clean_targets) ../results/ ../plots/

.PHONY cleanall: realclean	
	rm -rf $(cleanall_targets)

clean::
	rm -rf $(clean_targets)

