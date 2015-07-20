import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger
import wdmerger_star_diag_analysis.py as star_diag

# Set the name of the diagnostic output file

ncell = 256

diag_filename = "results/n" + str(ncell) + "/output/star_diag.out"

# Get the desired columns for the time and distance

time = star_diag.get_column('TIME',        diag_filename)
dist = star_diag.get_column('WD DISTANCE', diag_filename)

# Let's divide the time array by free-fall timescale

rotational_period = 100.0
t_ff = rotational_period / (4.0 * np.sqrt(2.0))

time = time / t_ff

# Divide distance by initial distance

dist = dist / dist[0]

# Generate the analytical result

d_exact = np.arange(0.0, 1.0, 0.001)

t_exact = np.arccos(np.sqrt(d_exact)) + np.sqrt(d_exact * (1.0 - d_exact))
t_exact *= 2.0 / np.pi

# Let's plot with reversed axes for aesthetic purposes

eps_filename = 'plots/freefall.eps'

# Stride for our arrays; the markers will be too dense if we plot every data point

stride = 25

plt.plot(dist[::stride], time[::stride], 'o', ms=12.0)
plt.plot(d_exact, t_exact, '-', lw=4.0)
plt.axis([0.0, 1.0, 0.0, 1.0])
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.xlabel("Distance / Initial Distance", fontsize=20)
plt.ylabel("Time / Free Fall Timescale", fontsize=20)
plt.tick_params(labelsize=16)
plt.savefig(eps_filename)
wdmerger.insert_commits_into_eps(eps_filename, diag_filename, 'diag')
