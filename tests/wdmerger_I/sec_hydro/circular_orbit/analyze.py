import yt
import numpy as np
from matplotlib import pyplot as plt

# Open up the diagnostic output for analysis

ncell = 128

diag_filename = "results/" + str(ncell) + "/wdmerger_diag.out"

diag_out = np.loadtxt(diag_filename)

# Now we need to obtain the column locations for the left and right centers of mass

diag_file = open(diag_filename,'r')
col_names = diag_file.readline().split('  ')
diag_file.close()

# Let's do some cleanup

col_names.pop(0)                                        # Get rid of the # at the beginning
col_names = [string.strip() for string in col_names]    # Remove any leading or trailing whitespace
col_names = filter(None, col_names)                     # Remove any remaining blank entries
col_names.pop(len(col_names)-1)                         # Remove the ending newline

# Obtain the time column, and the locations of the center of mass of both stars

col_t = col_names.index('TIME')
col_l = col_names.index('LEFT X COM')
col_r = col_names.index('RIGHT X COM')

time = diag_out[:,col_t]
dist = abs(diag_out[:,col_r] - diag_out[:,col_l])

# Normalize disatnce by initial distance

dist = dist / dist[0]

plt.plot(time, dist)
plt.xlabel("Time (s)")
plt.ylabel("Distance / Initial Distance")
plt.savefig('circular_orbit.png')
