# Problem-specific variables

mass_P:=0.90
mass_S:=0.81

rotational_period:=100.0

# We can work out the stopping time using the formula
# t_freefall = rotational_period / (4 * sqrt(2))
# bc doesn't have a native way to evaluate fractional exponents like square roots,
# so we perform the trick described here:
# http://www.linuxquestions.org/questions/programming-9/bc-and-exponents-containing-decimals-and-fractions-755260/
# The argument is that a^b = exp(b*ln(a)).
# By calling bc with the -l option we include the math library, 
# which knows how to evaluate exponents with e() and logarithms with l().
# The $$() tells the makefile that we want to expand the following in the shell.
# Of course, I could have just used sqrt(2) = 1.414, but where's the fun in that?

stop_time:=$$(echo "0.90 * $(rotational_period) / 4.0 / e(0.5*l(2))" | bc -l)

# Main make target

all: make_castro

	# Insert masses into probin file.
	sed -i "/mass_P/c mass_P = $(mass_P)" $(probin)
	sed -i "/mass_S/c mass_S = $(mass_S)" $(probin)

	# Insert the orbital period and stopping time into the inputs file
	sed -i "/castro.rotational_period/c castro.rotational_period = $(rotational_period)" $(inputs)
	sed -i "/stop_time/c stop_time = $(stop_time)" $(inputs)

# Pull in the main make targets

include $(WDMERGER_HOME)/job_scripts/makefile
