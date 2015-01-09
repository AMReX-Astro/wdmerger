import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger

# Set the name of the diagnostic output file

ncell = 256

diag_filename = "results/" + str(ncell) + "/wdmerger_diag.out"

# Get the desired columns for the time and the left and right centers of mass

time      = wdmerger.get_column('TIME',        diag_filename)
left_com  = wdmerger.get_column('LEFT X COM',  diag_filename)
right_com = wdmerger.get_column('RIGHT X COM', diag_filename)

dist = abs(right_com - left_com)

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

plt.plot(dist[::5], time[::5], 'o', d_exact, t_exact, '-')
plt.axis([0.0, 1.0, 0.0, 1.0])
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.xlabel("Distance / Initial Distance")
plt.ylabel("Time / Free Fall Timescale")
plt.savefig('freefall.png')
