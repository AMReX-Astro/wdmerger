SHELL=/bin/bash

UNAMEN := $(shell uname -n)

# Makefile options. We'll use the ?= method so that
# any including makefile can override these for custom tests
# (e.g. changing the number of dimensions or the EOS).

PRECISION ?= DOUBLE
PROFILE ?= FALSE

DEBUG ?= FALSE
BACKTRACE ?= FALSE
TEST ?= FALSE

USE_MPI ?= TRUE

# Some systems don't play as nicely with
# OpenMP, so disable it by default.

ifeq ($(findstring mira, $(UNAMEN)), mira)
	# ALCF Mira
	USE_OMP ?= FALSE
else
	USE_OMP ?= TRUE
endif

DIM ?= 3

USE_GRAV ?= TRUE
USE_REACT ?= TRUE
USE_MODELPARSER ?= TRUE
USE_ROTATION ?= TRUE

DIMENSION_AGNOSTIC ?= TRUE

USE_HYBRID_MOMENTUM ?= TRUE

# Set the default integrator for the burner.

INTEGRATOR_DIR ?= BS

# Override default compiler for certain machines.

ifeq ($(findstring titan, $(UNAMEN)), titan)
	# OLCF Titan
	COMP ?= Cray
	FCOMP ?= Cray
else ifeq ($(findstring h2o, $(UNAMEN)), h2o)
	# NCSA Blue Waters
	COMP ?= Cray
	FCOMP ?= Cray
else ifeq ($(findstring mira, $(UNAMEN)), mira)
	COMP ?= IBM
	FCOMP ?= IBM
else
	COMP ?= g++
	FCOMP ?= gfortran
endif

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_dir ?= helmholtz

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
Network_dir ?= aprox13

# If using general_null, make sure we have a reasonable set of isotopes.
ifeq ($(Network_dir), general_null)
	GENERAL_NET_INPUTS ?= $(MICROPHYSICS_HOME)/networks/$(Network_dir)/aprox13.net
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

include $(CASTRO_HOME)/Exec/Make.Castro

# Clean executable up with make clean;
# clean up all build files including objects with make realclean;
# clean up everything, including all results and analysis, with make cleanall.

clean_targets = *.mod helm_table.dat $(inputs) $(probin) 
cleanall_targets = $(clean_targets) ../results/ ../plots/

.PHONY cleanall: realclean	
	rm -rf $(cleanall_targets)

clean::
	rm -rf $(clean_targets)

