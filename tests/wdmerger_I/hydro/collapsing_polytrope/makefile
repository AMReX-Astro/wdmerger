# Problem-specific makefile variables
EOS_dir=polytrope

# Problem-specific variables

geometry.is_periodic = 1 1 1

castro.lo_bc       =  0 0 0
castro.hi_bc       =  0 0 0

geometry.prob_lo = 0.00e9 0.00e9 0.00e9
geometry.prob_hi = 5.12e9 5.12e9 5.12e9

castro.do_rotation = 0

# Main make target

all: make_castro

# Pull in main makefile

include $(WDMERGER_HOME)/job_scripts/makefile
