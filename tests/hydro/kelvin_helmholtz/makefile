# Problem-specific makefile variables
EOS_dir=gamma_law_general
DIM=2
USE_GRAV=FALSE

# Problem-specific variables

geometry.is_periodic = 1 1

castro.lo_bc       =  0 0
castro.hi_bc       =  0 0

castro.do_grav = 0

# Main make target

all: make_castro

# Pull in main makefile

include $(WDMERGER_HOME)/job_scripts/makefile
