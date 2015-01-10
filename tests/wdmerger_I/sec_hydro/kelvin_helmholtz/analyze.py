import yt
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import os
import wdmerger

problem_arr = [1, 2, 3]

ncell_arr = [64, 128, 256, 1024, 2048, 4096]

vel_arr = [0, 1, 3, 10, 30, 100]

plots_dir = 'plots/'

fig = plt.figure()

for problem in problem_arr:

    p = problem - 1 # For array indexing

    for v in vel_arr:

        for ncell in ncell_arr:

            # Set up the simulation times for generating plots

            if (problem == 3):
                plot_per = 0.1
                time_arr = [ 0.0, 1.5, 2.5, 4.7, 9.2 ]
            else:
                plot_per = 0.05
                time_arr = [ 0.0, 2.0 ]

            for t in time_arr:

                # Generate the plot name, and skip this iteration if it already exists.

                eps_filename = plots_dir + 'density_t' + str(t) + '_p' + str(problem) + '_v' + str(v) + '_n' + str(ncell) + '.eps'

                if (os.path.isfile(eps_filename)):
                    print "Plot with filename " + eps_filename + " already exists; skipping."
                    continue

                # First we load in the numerical data.

                dir = "results/problem" + str(problem) + "/velocity" + str(v) + "/" + str(ncell)

                # Make sure there's actually data for this combination of settings.
                # For example, we don't do every flow velocity at every resolution.

                if (not os.path.isdir(dir)):
                    continue

                print "Generating plot with filename " + eps_filename

                # Get the list of plotfiles in the directory.

                plotfiles = wdmerger.get_plotfiles(dir)

                # Figure out which one in the list we want by dividing the simulation time
                # by the interval between plotfile outputs.

                index = int(round(t / plot_per))

                # As a redundancy, guard against the possibility that this plotfile doesn't exist.
                # Could happen if the full run hasn't been completed yet.

                if (index >= len(plotfiles)):
                    print "Error: the plotfile does not exist. Skipping this iteration."
                    continue

                pf_name = dir + "/" + plotfiles[index]

                pf = yt.load(pf_name)

                # Load the data and create the plot

                p = yt.SlicePlot(pf, 'z', "density", width=(1.0, 'cm'))

                p.set_log("density", False)

                p.set_cmap(field="density", cmap='bwr')

                plot = p.plots['density']

                colorbar = plot.cb

                colorbar.solids.set_edgecolor("face") # Prevents lines from appearing in the colorbar at color boundaries, which can happen for the EPS output

                colorbar.set_clim(0.8,2.2)
                colorbar.set_ticks([0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2])

                p.save(eps_filename)

                wdmerger.insert_commits_into_eps(eps_filename, pf_name, 'plot')
                
