# Pull in the main make targets

DIM=2
TINY_PROFILE=TRUE
NETWORK_DIR=aprox13
USE_MPI=TRUE

USE_GPU=FALSE
USE_CUDA=FALSE
USE_HIP=FALSE

ifneq ($(ROCM_PATH),)
  USE_HIP=TRUE
else ifneq ($(CUDA_PATH),)
  USE_CUDA=TRUE
endif

USE_OMP=FALSE

PROBLEM_DIR=$(CASTRO_HOME)/Exec/science/wdmerger

include $(PROBLEM_DIR)/GNUmakefile
