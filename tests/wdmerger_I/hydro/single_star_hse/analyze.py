import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger

# Set the name of the diagnostic output file

ncell = 64

diag_filename = "results/" + str(ncell) + "/single_star.out"

# Obtain the time column, and the spherical radii denoting different density cutoffs

time = wdmerger.get_column('TIME', diag_filename)

nCols = 7

rad = np.zeros((len(time), nCols))

for n in range(nCols):
    rad[:,n] = wdmerger.get_column('1E' + str(n) + ' BOUNDARY RAD.')

    # Divide radii by initial radii

    rad[:,n] = rad[:,n] / rad[0,n]

eps_filename = "plots/single_star_hse.eps"

# Let's plot with reversed axis for aesthetic purposes

plt.plot(time, rad[:,2], label=r'$\rho=10^2$')
plt.plot(time, rad[:,4], label=r'$\rho=10^4$')
plt.plot(time, rad[:,6], label=r'$\rho=10^6$')
plt.legend(loc='lower right')
plt.xlabel("Time (s)")
plt.ylabel("Radius / Initial Radius")
plt.savefig(eps_filename)
wdmerger.insert_commits_into_eps(eps_filename, diag_filename, 'diag')
