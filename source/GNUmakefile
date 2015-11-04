SHELL=/bin/bash

# Makefile options. We'll use the ?= method so that
# any including makefile can override these for custom tests
# (e.g. changing the number of dimensions or the EOS).

PRECISION ?= DOUBLE
PROFILE ?= FALSE

DEBUG ?= FALSE
DEBUG = TRUE

USE_MPI ?= TRUE
USE_OMP ?= TRUE

DIM ?= 3

ifeq ($(DIM), 1)

  $(error Error: wdmerger is not configured to run in only one dimension. Please use DIM = 2 or DIM = 3)

endif

USE_GRAV ?= TRUE
USE_REACT ?= TRUE
USE_MODELPARSER ?= TRUE
USE_ROTATION ?= TRUE

DIMENSION_AGNOSTIC ?= TRUE

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

# This sets the EOS directory in $(MICROPHYSICS_DIR)/eos
EOS_dir ?= helmholtz

# This sets the network directory in $(MICROPHYSICS_DIR)/networks
Network_dir ?= aprox13

ifeq ($(USE_REACT), FALSE)
VODE_dir    := $(CASTRO_DIR)/Util/VODE
include $(VODE_dir)/Make.package
EXTERN_CORE += $(CASTRO_DIR)/Util/VODE
INCLUDE_LOCATIONS += $(VODE_dir) 
VPATH_LOCATIONS   += $(VODE_dir) 

BLAS_dir    := $(CASTRO_DIR)/Util/BLAS
include $(BLAS_dir)/Make.package
EXTERN_CORE += $(CASTRO_DIR)/Util/BLAS
INCLUDE_LOCATIONS += $(BLAS_dir) 
VPATH_LOCATIONS   += $(BLAS_dir) 
endif

# Location where our source files are coming from.

Blocs := $(WDMERGER_HOME)/source/
Bpack := $(WDMERGER_HOME)/source/Make.package

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

