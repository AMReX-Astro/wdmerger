# Problem-specific variables

mass_P:=0.90
mass_S:=-1.00

# Main make target

all: make_castro

	# Insert mass into probin file
	sed -i "/mass_P/c mass_P = $(mass_P)" $(probin)
	sed -i "/mass_S/c mass_S = $(mass_S)" $(probin)

# Pull in main makefile

include $(WDMERGER_HOME)/job_scripts/makefile
