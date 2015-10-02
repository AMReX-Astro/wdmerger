import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger

# List of possible load-balancing strategies
strat_arr = [ "SFC", "RRSFC", "PFC", "KNAPSACK", "ROUNDROBIN" ]

# Storage array for run-times
time_arr  = np.zeros( (len(strat_arr)) )

# Number of test repetitions
n_iters = 3
iter_arr = range(1, n_iters)
t_rep = np.zeros( (len(iter_arr)) )

# Output data file
dat_file_name = 'plots/timing.out'
dat_file = open(dat_file_name, 'w')

dat_file.write('Timing data for various load-balancing strategies:\n\n')

idx1 = 0

for strat in strat_arr:

    idx2 = 0

    for iter in iter_arr:

        # Open up the standard output for analysis.

        dir = "results/strat_" + strat_arr[idx1] + "/rep" + str(iter)

        output_filename = wdmerger.get_last_output(dir)

        # Get the timing information

        [ t_rep[idx2], _ ] = wdmerger.timing(output_filename)

        idx2 += 1

    # Now take the median of the various repetitions

    time_arr[idx1] = np.median(t_rep)

    dat_file.write(strat_arr[idx1] + '   ' + str(time_arr[idx1]) + '\n')

    idx1 += 1

minloc = time_arr.argmin()

print "The best distribution mapping strategy is: " + strat_arr[minloc]
print "See the " + dat_file_name + " file for the timing data."
