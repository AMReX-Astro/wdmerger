import yt
import numpy as np
from matplotlib import pyplot as plt
import os
import wdmerger

one_level_nproc_arr = np.array([32, 64, 128, 256])
one_level_time_arr = np.zeros(len(one_level_nproc_arr))

two_level_nproc_arr = np.array([256, 512, 1024, 2048])
two_level_time_arr = np.zeros(len(two_level_nproc_arr))

three_level_nproc_arr = np.array([2048, 4096, 8192, 16384])
three_level_time_arr = np.zeros(len(three_level_nproc_arr))

idx = 0

for nprocs in one_level_nproc_arr:

    # Open up the standard output for analysis.

    dir = "results/one_level/n" + str(nprocs)

    output_filename = wdmerger.get_last_output(dir)

    # Get the timing information

    [ one_level_time_arr[idx], _ ] = wdmerger.timing(output_filename)

    idx += 1

idx = 0

for nprocs in two_level_nproc_arr:

    # Open up the standard output for analysis.

    dir = "results/two_levels/n" + str(nprocs)

    output_filename = wdmerger.get_last_output(dir)

    # Get the timing information

    [ two_level_time_arr[idx], _ ] = wdmerger.timing(output_filename)

    idx += 1

idx = 0

for nprocs in three_level_nproc_arr:

    # Open up the standard output for analysis.

    dir = "results/three_levels/n" + str(nprocs)

    output_filename = wdmerger.get_last_output(dir)

    # Get the timing information

    [ three_level_time_arr[idx], _ ] = wdmerger.timing(output_filename)

    idx += 1

eps_filename = 'plots/strong_scaling.eps'

# Divide time array by initial value to get speedup

one_level_speedup_arr = (one_level_time_arr / one_level_time_arr[0])**(-1)
two_level_speedup_arr = (two_level_time_arr / two_level_time_arr[0])**(-1)
three_level_speedup_arr = (three_level_time_arr / three_level_time_arr[0])**(-1)

# Analytical result for perfect speedup

one_level_perfect_arr = np.ones(len(one_level_nproc_arr))
one_level_perfect_arr = one_level_perfect_arr * one_level_nproc_arr / one_level_nproc_arr[0]

two_level_perfect_arr = np.ones(len(two_level_nproc_arr))
two_level_perfect_arr = two_level_perfect_arr * two_level_nproc_arr / two_level_nproc_arr[0]

three_level_perfect_arr = np.ones(len(three_level_nproc_arr))
three_level_perfect_arr = three_level_perfect_arr * three_level_nproc_arr / three_level_nproc_arr[0]

# Move the figure up and to the right so that the larger 
# axis labels fit inside the plot.

ax = plt.gca()
box=[0.125, 0.125, 0.85, 0.85]
ax.set_position(box)

plt.plot(one_level_nproc_arr, one_level_perfect_arr, color='0.25', linestyle='-', lw=3.0, label="Perfect strong scaling")
plt.plot(one_level_nproc_arr, one_level_speedup_arr, 'go', ms=12.0, label="One level")
plt.plot(two_level_nproc_arr, two_level_perfect_arr, color='0.25', linestyle='-', lw=3.0)
plt.plot(two_level_nproc_arr, two_level_speedup_arr, 'bD', ms=12.0, label="Two levels")
plt.plot(three_level_nproc_arr, three_level_perfect_arr, color='0.25', linestyle='-', lw=3.0)
plt.plot(three_level_nproc_arr, three_level_speedup_arr, 'rs', ms=12.0, label="Three levels")
plt.xlabel("Number of processors", fontsize=20.0)
plt.ylabel("Relative speedup", fontsize=20.0)
plt.xscale("log", basex=2)
plt.yscale("log")
plt.tick_params(labelsize=16,pad=10)
plt.axis([0.8 * one_level_nproc_arr[0], 1.2 * three_level_nproc_arr[-1], 0.8 * one_level_perfect_arr[0], 2.0 * two_level_perfect_arr[-1]])
plt.legend(loc="upper right", numpoints=1)
plt.savefig(eps_filename)
wdmerger.insert_commits_into_eps(eps_filename, output_filename, 'info')
