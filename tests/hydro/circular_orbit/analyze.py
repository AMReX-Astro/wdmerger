import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger

# Set the name of the diagnostic output file

ncell = 128

diag_filename = "results/" + str(ncell) + "/wdmerger_diag.out"

# Obtain the time column, and the distance between the stars

time = wdmerger.get_column('TIME',        diag_filename)
dist = wdmerger.get_column('WD DISTANCE',  diag_filename)

# Normalize disatnce by initial distance

dist = dist / dist[0]

eps_filename = 'plots/circular_orbit.eps'

plt.plot(time, dist)
plt.xlabel("Time (s)")
plt.ylabel("Distance / Initial Distance")
plt.savefig(eps_filename)
wdmerger.insert_commits_into_eps(eps_filename, diag_filename, 'diag')
