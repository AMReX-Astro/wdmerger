import yt
import numpy as np
from matplotlib import pyplot as plt

# Open up the diagnostic output for analysis

ncell = 64

diag_filename = "results/" + str(ncell) + "/single_star.out"

diag_out = np.loadtxt(diag_filename)

# Now we need to obtain the column names

diag_file = open(diag_filename,'r')
col_names = diag_file.readline().split('  ')
diag_file.close()

# Let's do some cleanup

col_names.pop(0)                                        # Get rid of the # at the beginning
col_names = [string.strip() for string in col_names]    # Remove any leading or trailing whitespace and newlines
col_names = filter(None, col_names)                     # Remove any remaining blank entries

# Obtain the time column, and the spherical radii denoting different density cutoffs

time = diag_out[:,col_names.index('TIME')]

nCols = 7

rad = np.zeros((len(time), nCols))

rad[:,0] = diag_out[:,col_names.index('1E0 BOUNDARY RAD.')]
rad[:,1] = diag_out[:,col_names.index('1E1 BOUNDARY RAD.')]
rad[:,2] = diag_out[:,col_names.index('1E2 BOUNDARY RAD.')]
rad[:,3] = diag_out[:,col_names.index('1E3 BOUNDARY RAD.')]
rad[:,4] = diag_out[:,col_names.index('1E4 BOUNDARY RAD.')]
rad[:,5] = diag_out[:,col_names.index('1E5 BOUNDARY RAD.')]
rad[:,6] = diag_out[:,col_names.index('1E6 BOUNDARY RAD.')]

# Divide radii by initial radii

for j in range(nCols):
    rad[:,j] = rad[:,j] / rad[0,j]

# Let's plot with reversed axis for aesthetic purposes

plt.plot(time, rad[:,2], label=r'$\rho=10^2$')
plt.plot(time, rad[:,4], label=r'$\rho=10^4$')
plt.plot(time, rad[:,6], label=r'$\rho=10^6$')
plt.legend(loc='lower right')
plt.xlabel("Time (s)")
plt.ylabel("Radius / Initial Radius")
plt.savefig('single_star_hse.png')
