CASTRO_DIR ?= /path/to/Castro

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

USE_MPI    = TRUE
USE_OMP    = TRUE

DIM        = 3

COMP	   = g++
FCOMP	   = gfortran

USE_GRAV   = TRUE
USE_REACT  = FALSE
USE_MODELPARSER  = FALSE
USE_ROTATION = FALSE

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := helmeos

# This sets the network directory in $(NETWORK_HOME)
Network_dir := general_null
GENERAL_NET_INPUTS := $(CASTRO_DIR)/Networks/general_null/triple_alpha_plus_o.net 

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_DIR)/Exec/Make.Castro
