# Problem-specific variables

mass_P = 0.90
mass_S = 0.81

castro.do_rotation = 0
castro.rotational_period = 100.0

amr.plot_per = 0.1
amr.check_per = 0.1

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

stop_time = $$(echo "0.90 * $(castro.rotational_period) / 4.0 / e(0.5*l(2))" | bc -l)

# Main make target

all: make_castro

# Pull in the main make targets

include $(WDMERGER_HOME)/job_scripts/makefile
