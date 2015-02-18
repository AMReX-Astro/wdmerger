CASTRO_DIR ?= /path/to/Castro

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

USE_MPI = TRUE
USE_OMP = FALSE

DIM        = 3

COMP	   = g++
FCOMP	   = gfortran

USE_GRAV = TRUE
USE_REACT = FALSE
USE_MODELPARSER = TRUE
USE_ROTATION = TRUE
#USE_TRACING = TRUE

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := helmeos

# This sets the network directory in $(NETWORK_HOME)
Network_dir := general_null
GENERAL_NET_INPUTS := $(CASTRO_DIR)/Networks/general_null/triple_alpha_plus_o.net 

diagpack := $(WDMERGER_HOME)/source/diagnostics/Make.package
diaglocs := $(WDMERGER_HOME)/source/diagnostics/

modspack := $(WDMERGER_HOME)/source/modules/Make.package
modslocs := $(WDMERGER_HOME)/source/modules/

srcepack := $(WDMERGER_HOME)/source/Make.package
srcelocs := $(WDMERGER_HOME)/source/

Bpack   := $(srcepack) $(modspack) $(diagpack)
Blocs   := $(srcelocs) $(modslocs) $(diaglocs)

BUILD_GIT_NAME := wdmerger
BUILD_GIT_DIR := $(WDMERGER_HOME)

include $(CASTRO_DIR)/Exec/Make.Castro

