# Problem-specific variables

mass_P:=0.90
mass_S:=0.60

rotational_period:=100.0

max_step:=10

# Main make target

all: make_castro

	# Insert masses into probin file.
	sed -i "/mass_P/c mass_P = $(mass_P)" $(probin)
	sed -i "/mass_S/c mass_S = $(mass_S)" $(probin)

	# Insert the orbital period and stopping time into the inputs file
	sed -i "/castro.rotational_period/c castro.rotational_period = $(rotational_period)" $(inputs)
	sed -i "/max_step/c max_step = $(max_step)" $(inputs)

# Pull in main makefile

include $(WDMERGER_HOME)/job_scripts/makefile
