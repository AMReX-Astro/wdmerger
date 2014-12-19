import yt
import numpy as np
from matplotlib import pyplot as plt

# Open up the diagnostic output for analysis

ncell = 64

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

# Obtain the time column, and the locations of the center of mass of both stars

col_t = col_names.index('TIME')
col_l = col_names.index('LEFT X COM')
col_r = col_names.index('RIGHT X COM')

time = diag_out[:,col_t]
dist = abs(diag_out[:,col_r] - diag_out[:,col_l])

# Let's divide the time array by free-fall timescale

rotational_period = 100.0
t_ff = rotational_period / (4.0 * np.sqrt(2.0))

time = time / t_ff

# Divide distance by initial distance

dist = dist / dist[0]

# Generate the analytical result.

d_exact = np.arange(0.0, 1.0, 0.001)

t_exact = np.arccos(np.sqrt(d_exact)) + np.sqrt(d_exact * (1.0 - d_exact))
t_exact *= 2.0 / np.pi

# Let's plot with reversed axis for aesthetic purposes

plt.plot(dist[::5], time[::5], 'o', d_exact, t_exact, '-')
plt.axis([0.0, 1.0, 0.0, 1.0])
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.xlabel("Distance / Initial Distance")
plt.ylabel("Time / Free Fall Timescale")
plt.savefig('freefall.png')
