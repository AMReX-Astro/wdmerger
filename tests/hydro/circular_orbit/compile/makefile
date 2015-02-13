# Problem-specific variables

mass_P = 0.90
mass_S = 0.60

castro.rotational_period = 100.0

num_periods = 50

stop_time = $$(echo "$(num_periods) * $(castro.rotational_period)" | bc -l)

# Main make target

all: make_castro

# Pull in main makefile

include $(WDMERGER_HOME)/job_scripts/makefile
