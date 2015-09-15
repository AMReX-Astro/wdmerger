import yt
import numpy as np
from matplotlib import pyplot as plt
import os
import wdmerger

small_nproc_arr = np.array([64, 128, 256, 512, 1024, 2048])
small_time_arr = np.zeros(len(small_nproc_arr))

large_nproc_arr = np.array([2048, 4096, 8192, 16384])
large_time_arr = np.zeros(len(large_nproc_arr))

idx = 0

for nprocs in small_nproc_arr:

    # Open up the standard output for analysis.

    dir = "results/two_levels/n" + str(nprocs)

    output_filename = wdmerger.get_last_output(dir)

    # Get the timing information

    [ small_time_arr[idx], _ ] = wdmerger.timing(output_filename)

    idx += 1

idx = 0

for nprocs in large_nproc_arr:

    # Open up the standard output for analysis.

    dir = "results/three_levels/n" + str(nprocs)

    output_filename = wdmerger.get_last_output(dir)

    # Get the timing information

    [ large_time_arr[idx], _ ] = wdmerger.timing(output_filename)

    idx += 1

eps_filename = 'plots/strong_scaling.eps'

# Divide time array by initial value to get speedup

small_speedup_arr = (small_time_arr / small_time_arr[0])**(-1)
large_speedup_arr = (large_time_arr / large_time_arr[0])**(-1)

# Analytical result for perfect speedup

small_perfect_arr = np.ones(len(small_nproc_arr))
small_perfect_arr = small_perfect_arr * small_nproc_arr / small_nproc_arr[0]

large_perfect_arr = np.ones(len(large_nproc_arr))
large_perfect_arr = large_perfect_arr * large_nproc_arr / large_nproc_arr[0]

# Move the figure up and to the right so that the larger 
# axis labels fit inside the plot.

ax = plt.gca()
box=[0.125, 0.125, 0.85, 0.85]
ax.set_position(box)

plt.plot(small_nproc_arr, small_perfect_arr, color='0.25', linestyle='-', lw=3.0, label="Perfect strong scaling")
plt.plot(small_nproc_arr, small_speedup_arr, 'bo', ms=12.0, label="Two levels")
plt.plot(large_nproc_arr, large_perfect_arr, color='0.25', linestyle='-', lw=3.0)
plt.plot(large_nproc_arr, large_speedup_arr, 'rs', ms=12.0, label="Three levels")
plt.xlabel("Number of processors", fontsize=20.0)
plt.ylabel("Relative speedup", fontsize=20.0)
plt.xscale("log", basex=2)
plt.yscale("log")
plt.tick_params(labelsize=16)
plt.axis([0.8 * small_nproc_arr[0], 1.2 * large_nproc_arr[-1], 0.8 * small_perfect_arr[0], 2.0 * small_perfect_arr[-1]])
plt.legend(loc="upper right", numpoints=1)
plt.savefig(eps_filename)
wdmerger.insert_commits_into_eps(eps_filename, output_filename, 'info')
