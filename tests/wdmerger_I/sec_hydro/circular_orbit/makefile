# Problem-specific variables

mass_P:=0.90
mass_S:=0.60

rotational_period:=100.0

num_periods:=50

stop_time:=$$(echo "$(num_periods) * $(rotational_period)" | bc -l)

# Main make target

all: make_castro

	# Insert masses into probin file.
	sed -i "/mass_P/c mass_P = $(mass_P)" $(probin)
	sed -i "/mass_S/c mass_S = $(mass_S)" $(probin)

	# Insert the orbital period and stopping time into the inputs file
	sed -i "/castro.rotational_period/c castro.rotational_period = $(rotational_period)" $(inputs)
	sed -i "/stop_time/c stop_time = $(stop_time)" $(inputs)

# Pull in main makefile

include $(WDMERGER_HOME)/job_scripts/makefile
