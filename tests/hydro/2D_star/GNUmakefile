
PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

#DIM        = 1
DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE
USE_OMP    = FALSE

# This sets the EOS directory in Castro/EOS
EOS_dir     := helmeos

# This sets the network directory in Castro/Networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(CASTRO_DIR)/Networks/$(Network_dir)/ignition.net

Bpack   := ./Make.package
Blocs   := 

include $(CASTRO_DIR)/Exec/Make.Castro
