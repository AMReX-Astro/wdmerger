import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger

ncell_arr = [256, 1024, 4096]
omp_arr = [1, 2, 4, 8]

time_arr = np.zeros( (len(ncell_arr),len(omp_arr)) ) # Average timestep
hydro_arr = np.zeros( (len(ncell_arr),len(omp_arr)) ) # Average timestep without counting gravity

index = 0

idx1 = 0

for n in ncell_arr:

    idx2 = 0

    minloc = 0

    for nthreads in omp_arr:

        dir = 'results/n' + str(n) + '/omp' + str(nthreads)

        # Open up the standard output for analysis.

        output = wdmerger.get_last_output(dir)

        [time_arr[idx1,idx2], hydro_arr[idx1,idx2]] = wdmerger.timing(output)

        if (time_arr[idx1,idx2] < time_arr[idx1,minloc]):
            minloc = idx2

        idx2 += 1

    idx1 += 1

    # Now print out our best combination of settings

    print 'For ' + str(n) + ' zones per spatial dimension, the optimal number of threads per MPI task is ' + str(omp_arr[minloc])
