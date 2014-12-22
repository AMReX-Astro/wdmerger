import yt
import numpy as np
from matplotlib import pyplot as plt
import os

num_problems = 2

ncell = 64

vel_arr = [0, 1, 3, 10, 30, 100]

for p in range(num_problems):

    problem = p + 1

    for v in vel_arr:

        # First we load in the numerical data. We want to use the last plotfile
        # since that is at the time of interest.

        dir = "results/problem" + str(problem) + "/velocity" + str(v) + "/" + str(ncell)

        dir_contents = os.listdir(dir)

        # Strip out non-plotfiles, then select the last one in the list.

        plotfiles = sorted(filter(lambda s: s[0:3] == 'plt', dir_contents))

        pf_name = dir + "/" + plotfiles[len(plotfiles)-1]

        pf = yt.load(pf_name)

        plot = yt.SlicePlot(pf, 'z', "density", width=(1.0, 'cm'))

        file_name = 'density_t2_p_' + str(problem) + '_v_' + str(v) + '_n_' + str(ncell) + '.eps'

        plot.save(file_name)
