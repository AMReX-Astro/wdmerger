import yt
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import os
import wdmerger

problem_arr = [1, 2, 3, 4]

ncell_arr = [64, 128, 256, 512, 1024, 2048, 4096]

vel_arr = [0, 1, 3, 10, 30, 100]

plots_dir = 'plots/'

fig = plt.figure()

# Generate an individual image for every problem at every resolution

for problem in problem_arr:

    p = problem - 1 # For array indexing

    for v in vel_arr:

        for ncell in ncell_arr:

            # Set up the simulation times for generating plots

            if (problem == 3):
                plot_per = 0.1
                time_arr = [ 0.0, 1.5, 2.0, 2.5, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 4.7, 9.2 ]
            else:
                plot_per = 0.05
                time_arr = [ 0.0, 2.0 ]

            for t in time_arr:

                dir = "results/problem" + str(problem) + "/velocity" + str(v) + "/n" + str(ncell) + "/output"

                # Make sure there's actually data for this combination of settings.
                # For example, we don't do every flow velocity at every resolution.

                if (not os.path.isdir(dir)):
                    continue

                # Generate the plot name, and skip this iteration if it already exists.

                eps_filename = plots_dir + 'kh_t' + str(t) + '_p' + str(problem) + '_v' + str(v) + '_n' + str(ncell) + '.eps'

                if (os.path.isfile(eps_filename)):
                    print "Plot with filename " + eps_filename + " already exists; skipping."
                    continue

                print "Generating plot with filename " + eps_filename

                # First we load in the numerical data.

                # Get the list of plotfiles in the directory.

                plotfiles = wdmerger.get_plotfiles(dir)

                # Figure out which one in the list we want by searching through
                # and checking the Header file. We'll take the first plotfile 
                # after the desired time.

                plotfile_times = np.zeros(len(plotfiles))

                for n in range(len(plotfiles)):                
                    plotfile_times[n] = wdmerger.get_time_from_plotfile(dir + '/' + plotfiles[n])

                for n in range(len(plotfile_times)):

                    index = n

                    time_curr = plotfile_times[n]

                    if (n == 0):
                        time_old = time_curr
                    else:
                        time_old = plotfile_times[n-1]

                    if (time_old <= t and time_curr >= t):
                        break
                
                pf_name = dir + "/" + plotfiles[index]

                pf = yt.load(pf_name)

                # Load the data and create the plot

                p = yt.SlicePlot(pf, 'z', "density", width=(1.0, 'cm'))

                p.set_log("density", False)

                p.set_cmap(field="density", cmap='bwr')

                p.set_zlim('density',0.75,2.25)

                plot = p.plots['density']

                cb = plot.cb

                # Prevents lines from appearing in the colorbar at color boundaries, which can happen for the EPS output
                cb.solids.set_rasterized(True)

                p.save(eps_filename)

                wdmerger.insert_commits_into_eps(eps_filename, pf_name, 'plot')
                


# Now create the collated frames for each problem, for all the moderate resolutions.

ncell_arr = [ 64, 128, 256, 512 ]

pROW = len(vel_arr)
pCOL = len(ncell_arr)

grid_padding = 0.01 # padding between axes in inches
label_mode = 'L' # Label_mode='L" will only apply labels to the outer sides of the grid of figures


font = {'family' : 'serif',
        'fontname' : 'times',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16}

for problem in problem_arr:

    p = problem - 1 # For array indexing

    # Set up the simulation times for generating plots

    if (problem == 3):
        plot_per = 0.1
        time_arr = [ 0.0, 1.5, 2.0, 2.5, 4.7, 9.2 ]
    else:
        plot_per = 0.05
        time_arr = [ 0.0, 2.0 ]

    for t in time_arr:

        # Generate the plot name, and skip this iteration if it already exists.

        eps_filename = plots_dir + 'kh_t' + str(t) + '_p' + str(problem) + '_low_res_collage.eps'

        if (os.path.isfile(eps_filename)):
            print "Plot with filename " + eps_filename + " already exists; skipping."
            continue

        print "Generating plot with filename " + eps_filename

        p_idx = 0

        fig = plt.figure()

        grid = AxesGrid(fig,(0.025,0.10,0.975,0.8),
                        nrows_ncols=(pROW,pCOL),
                        axes_pad=grid_padding,
                        label_mode=label_mode,
                        share_all="true",
                        cbar_location="right",
                        cbar_mode="single")

        plt.figtext(0.475,0.050,'x (cm)', fontdict=font)
        plt.figtext(0.200,0.525,'y (cm)', fontdict=font, rotation='vertical')

        plt.figtext(0.290,0.925,r'n = 64',  fontdict=font)        
        plt.figtext(0.395,0.925,r'n = 128', fontdict=font)
        plt.figtext(0.510,0.925,r'n = 256', fontdict=font)
        plt.figtext(0.625,0.925,r'n = 512', fontdict=font)

        plt.figtext(0.125,0.825,r'v = 0',   fontdict=font)
        plt.figtext(0.125,0.690,r'v = 1',   fontdict=font)
        plt.figtext(0.125,0.555,r'v = 3',   fontdict=font)
        plt.figtext(0.125,0.420,r'v = 10',  fontdict=font)
        plt.figtext(0.125,0.295,r'v = 30',  fontdict=font)
        plt.figtext(0.125,0.160,r'v = 100', fontdict=font)

        for v in vel_arr:

            for ncell in ncell_arr:

                # First we load in the numerical data.

                dir = "results/problem" + str(problem) + "/velocity" + str(v) + "/n" + str(ncell) + "/output"

                # Make sure there's actually data for this combination of settings.
                # For example, we don't do every flow velocity at every resolution.

                if (not os.path.isdir(dir)):
                    print "Needed data at velocity " + str(v) + " and resolution " + str(ncell) + " for problem " + str(problem) + " does not exist; quitting."
                    exit()

                # Get the list of plotfiles in the directory.

                plotfiles = wdmerger.get_plotfiles(dir)

                # Figure out which one in the list we want by searching through
                # and checking the Header file. We'll take the first plotfile 
                # after the desired time.

                plotfile_times = np.zeros(len(plotfiles))

                for n in range(len(plotfiles)):                
                    plotfile_times[n] = wdmerger.get_time_from_plotfile(dir + '/' + plotfiles[n])

                for n in range(len(plotfile_times)):

                    index = n

                    time_curr = plotfile_times[n]

                    if (n == 0):
                        time_old = time_curr
                    else:
                        time_old = plotfile_times[n-1]

                    if (time_old <= t and time_curr >= t):
                        break
                
                pf_name = dir + "/" + plotfiles[index]

                pf = yt.load(pf_name)

                # Load the data and create the plot

                p = yt.SlicePlot(pf, 'z', "density", width=(1.0, 'cm'), origin='left-domain')

                p.set_zlim('density', 0.75, 2.25)

                p.set_log("density", False)

                p.set_cmap(field="density", cmap='bwr')

                plot = p.plots['density']

                plot.figure = fig

                plot.axes = grid[p_idx].axes
                plot.cax = grid.cbar_axes[p_idx]

                # Redraw plot onto the AxesGrid

                p._setup_plots()

                plot = p.plots['density']

                cb = plot.cb

                # Prevents lines from appearing in the colorbar at color boundaries, which can happen for the EPS output
                cb.solids.set_rasterized(True)

                grid[p_idx].axes.set_xlabel('')
                grid[p_idx].axes.set_ylabel('')

                grid[p_idx].xaxis.set_ticks([0.25, 0.75])
                grid[p_idx].yaxis.set_ticks([0.25, 0.75])

                p_idx += 1

        # Save output

        plt.savefig(eps_filename)
        wdmerger.insert_commits_into_eps(eps_filename, pf_name, 'plot')



# Now generate a high resolution time series for problem 3

time_arr = [3.0, 3.2, 3.4, 3.6]

problem = 3
p = 2
v = 0

ncell_arr = [1024, 2048, 4096]

pROW = len(ncell_arr)
pCOL = len(time_arr)

eps_filename = plots_dir + 'kh_p3_high_res_collage.eps'

if (os.path.isfile(eps_filename)):
    print "Plot with filename " + eps_filename + " already exists; skipping."
    exit()

p_idx = 0

fig = plt.figure()

grid = AxesGrid(fig,(0.20,0.10,0.7,0.8),
                nrows_ncols=(pROW,pCOL),
                axes_pad=grid_padding,
                label_mode=label_mode,
                share_all="true",
                cbar_location="right",
                cbar_mode="single")

plt.figtext(0.515,0.150,'x (cm)', fontdict=font)
plt.figtext(0.125,0.525,'y (cm)', fontdict=font, rotation='vertical')

plt.figtext(0.250,0.825,r't = 3.0', fontdict=font)
plt.figtext(0.415,0.825,r't = 3.2', fontdict=font)
plt.figtext(0.585,0.825,r't = 3.4', fontdict=font)
plt.figtext(0.755,0.825,r't = 3.6', fontdict=font)

plt.figtext(0.01,0.7,r'n = 1024', fontdict=font)
plt.figtext(0.01,0.5,r'n = 2048', fontdict=font)
plt.figtext(0.01,0.3,r'n = 4096', fontdict=font)

for ncell in ncell_arr:

    for t in time_arr:

        # First we load in the numerical data.

        dir = "results/problem" + str(problem) + "/velocity" + str(v) + "/n" + str(ncell) + "/output"

        # Make sure there's actually data for this combination of settings.
        # For example, we don't do every flow velocity at every resolution.

        if (not os.path.isdir(dir)):
            print "Needed data at velocity " + str(v) + " and resolution " + str(ncell) + " for problem " + str(problem) + " does not exist; quitting."
            exit()

        # Get the list of plotfiles in the directory.

        plotfiles = wdmerger.get_plotfiles(dir)

        # Figure out which one in the list we want by searching through
        # and checking the Header file. We'll take the first plotfile 
        # after the desired time.

        plotfile_times = np.zeros(len(plotfiles))

        for n in range(len(plotfiles)):                
            plotfile_times[n] = wdmerger.get_time_from_plotfile(dir + '/' + plotfiles[n])

        for n in range(len(plotfile_times)):

            index = n

            time_curr = plotfile_times[n]

            if (n == 0):
                time_old = time_curr
            else:
                time_old = plotfile_times[n-1]

            if (time_old <= t and time_curr >= t):
                break
                
        pf_name = dir + "/" + plotfiles[index]

        pf = yt.load(pf_name)

        # Load the data and create the plot

        p = yt.SlicePlot(pf, 'z', "density", width=(1.0, 'cm'), origin='left-domain')

        p.set_zlim('density', 0.75, 2.25)

        p.set_log("density", False)

        p.set_cmap(field="density", cmap='bwr')

        plot = p.plots['density']

        plot.figure = fig

        plot.axes = grid[p_idx].axes
        plot.cax = grid.cbar_axes[p_idx]

        # Redraw plot onto the AxesGrid

        p._setup_plots()

        plot = p.plots['density']

        cb = plot.cb

        # Prevents lines from appearing in the colorbar at color boundaries, which can happen for the EPS output
        cb.solids.set_rasterized(True)

        grid[p_idx].axes.set_xlabel('')
        grid[p_idx].axes.set_ylabel('')

        grid[p_idx].xaxis.set_ticks([0.25, 0.75])
        grid[p_idx].yaxis.set_ticks([0.25, 0.75])

        p_idx += 1

# Save output

plt.savefig(eps_filename)
wdmerger.insert_commits_into_eps(eps_filename, pf_name, 'plot')
