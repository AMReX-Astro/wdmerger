import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger

# Set the name of the diagnostic output file

ncell = 128

diag_filename = "results/" + str(ncell) + "/wdmerger_diag.out"

# Obtain the time column, and the locations of the center of mass of both stars

time      = wdmerger.get_column('TIME',        diag_filename)
left_com  = wdmerger.get_column('LEFT X COM',  diag_filename)
right_com = wdmerger.get_column('RIGHT X COM', diag_filename)

dist = abs(right_com - left_com)

# Normalize disatnce by initial distance

dist = dist / dist[0]

plt.plot(time, dist)
plt.xlabel("Time (s)")
plt.ylabel("Distance / Initial Distance")
plt.savefig('circular_orbit.png')
