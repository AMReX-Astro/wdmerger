import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger

ncell_arr = [ 128, 256,  512,  1024,  2048]
nproc_arr = [   8,  64,  512,  4096, 32768]
time_arr  = np.zeros(len(nproc_arr))

idx = 0

for ncell in ncell_arr:

    # Open up the standard output for analysis.

    dir = "results/n" + str(ncell)

    output_filename = wdmerger.get_last_output(dir)

    # Get the timing information

    [ time_arr[idx], _ ] = wdmerger.timing(output_filename)

    idx += 1

eps_filename = 'plots/weak_scaling.eps'

# Divide time array by value for smallest number of processors

scaling_arr = time_arr / time_arr[0]

# Analytical result for perfect weak scaling

perfect_arr = np.ones(len(nproc_arr))

# Move the figure up and to the right so that the larger 
# axis labels fit inside the plot.

ax = plt.gca()
box=[0.125, 0.125, 0.85, 0.75]
ax.set_position(box)

ax.plot(ncell_arr, perfect_arr, color='0.25', linestyle='-', lw=3.0)
ax.plot(ncell_arr, scaling_arr, 'bo', ms=12.0)
plt.xlabel("Number of zones per dimension", fontsize=20.0)
plt.ylabel("Relative Timestep Length", fontsize=20.0)
plt.xscale("log", basex=2)
ax.tick_params(width=2,length=6,labelsize=16,pad=10)
ax.set_xticks(ncell_arr)
ax.set_yticks(np.arange(0.0,2.25,0.25))
ax.set_xlim(0.8 * ncell_arr[0], 1.2 * ncell_arr[-1])
ax.set_ylim(-0.025, 2.025)
ax.set_xscale("log", basex=2)

# Add a second x-axis so that we can plot the number of processors
# http://matplotlib.org/faq/howto_faq.html#multiple-y-axis-scales
# http://stackoverflow.com/questions/10514315/how-to-add-a-second-x-axis-in-matplotlib

ax2 = ax.twiny()
ax2.set_position(box)
ax2.set_xlabel("Number of processors", fontsize=20.0)
ax2.tick_params(width=2,length=6,labelsize=16,pad=10)
ax2.set_xscale("log", basex=2)
ax2.set_xticks(ncell_arr)
ax2.set_xbound(ax.get_xbound())
ax2.set_xticklabels(nproc_arr)

plt.savefig(eps_filename)
wdmerger.insert_commits_into_eps(eps_filename, output_filename, 'info')
