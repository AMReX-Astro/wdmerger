import os
import wdmerger
import numpy as np
import matplotlib.pyplot as plt



#
# Get the data column from a wdmerger diagnostic output file whose header name is col_name.
#

def get_column(col_name, diag_filename):

    # Open up the file for reading. Get the names of the columns, as well as a 2D list
    # with all the data.

    diag_file = open(diag_filename,'r')

    vc_line = 'git'

    # Skip the first few lines, they store the version control information

    line = diag_file.readline()

    while (line.split()[2] == vc_line):
        line = diag_file.readline()

    # The very next line will be the column headers

    col_names = line.split('  ')

    # Now read in the data

    diag_list = diag_file.readlines()
    data = []
    for line in diag_list:
        data.append( map(float, line[0:-1].split()) )

    diag_file.close()

    # Convert the data list into a 2D numpy array.

    data = np.array(data)

    # Let's do some cleanup.

    col_names.pop(0)                                        # Get rid of the # at the beginning
    col_names = [string.strip() for string in col_names]    # Remove any leading or trailing whitespace
    col_names = filter(None, col_names)                     # Remove any remaining blank entries

    # Obtain the column index and then return the column with that index.

    col_index = col_names.index(col_name)

    return data[:,col_index]



#
# Plot the change in distance of the white dwarf binary over time.
# This can accept multiple runs, expressed as a list.
#

def plot_wd_distance(diag_filenames, output_filename, labels=''):

    # If we got only a single diagnostic file,
    # create a list out of it for the coming loop.

    if (type(diag_filenames) not in [list, tuple]):
        diag_filenames = [ diag_filenames ]
        labels = [ labels ]

    # Cycle through markers for multiple curves on the same plot.

    markers = ['o', '+', '.', ',', '*']

    j = 0

    for diag_filename, label in zip(diag_filenames, labels):

        diag_file = open(diag_filename, 'r')

        time = get_column("TIME", diag_filename)
        dist = get_column("WD DISTANCE",diag_filename)

        # Normalize distance by initial distance.

        dist = dist / dist[0]

        # Normalize time by rotational period.

        rot_period = 0.0

        dir = os.path.dirname(diag_filename)

        inputs_filename = dir + '/' + wdmerger.get_inputs_filename(dir)
        rot_period = wdmerger.get_inputs_var(inputs_filename, "castro.rotational_period")

        if (rot_period > 0.0):
            time = time / rot_period
            xlabel = "Time / Rotational Period"
        else:
            xlabel = "Time (s)"
        ylabel = "WD Distance / Initial Distance"

        if (label is not ''):
            plot_label = label

        if (len(diag_filenames) > 1):
            marker = markers[j]
        else:
            marker = ''

        # Use a markevery stride because of how densely packed the star_diag.out data usually is.

        plt.plot(time, dist, marker = marker, ms = 12.0, markevery = 500, lw = 4.0, label=plot_label)

        j += 1

    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)

    # Use the 'best' location for the legend, since for a generic function like this
    # it is hard to know ahead of time where the legend ought to go.
    # The alpha value controls the transparency, since we may end up covering some data.

    if (label is not ''):
        plt.legend(loc='best', fancybox=True, framealpha=0.5)

    # The padding ensures that the lower-left ticks on the x- and y-axes don't overlap.

    plt.tick_params(labelsize=16, pad=10)

    # We have now increased the size of both the ticks and the axis labels,
    # which may have caused the latter to fall off the plot. Use tight_layout
    # to automatically adjust the plot to fix this.

    plt.tight_layout()

    # Save it into our designated file, which is usually EPS format.

    plt.savefig(output_filename)

    # Insert git commit hashes into this file from the various code sources.

    wdmerger.insert_commits_into_eps(output_filename, diag_filename, 'diag')

    plt.close()



#
# Plot the 2D location of both WDs.
# This can accept multiple runs, expressed as a list, and will 
# use subplots to separate them.
# This currently assumes you are using rot_axis = 3 and so 
# the stars are moving on the xy-plane.
#

def plot_wd_location(diag_filenames, output_filename):

    # If we got only a single diagnostic file,
    # create a list out of it for the coming loop.

    if (type(diag_filenames) not in [list, tuple]):
        diag_filenames = [ diag_filenames ]

    # Cycle through markers for multiple curves on the same plot.

    N = len(diag_filenames)

    # Determine the number of subplots and their layout.

    nRows = 1
    nCols = 1

    if (N == 2):
        nCols = 2
    elif (N == 3):
        nCols = 3
    elif (N == 4):
        nRows = 2
        nCols = 2

    if (N > 1):  
        fig, ax = plt.subplots(nRows, nCols)
    else:
        fig = plt.gcf()
        ax  = plt.gca()

    row = 0
    col = 0
    j   = 0

    labels = ['(a)', '(b)', '(c)', '(d)']

    for diag_filename in diag_filenames:

        diag_file = open(diag_filename, 'r')

        # x and y locations of each star

        xp   = get_column("PRIMARY X COM",diag_filename)
        yp   = get_column("PRIMARY Y COM",diag_filename)
        xs   = get_column("SECONDARY X COM",diag_filename)
        ys   = get_column("SECONDARY Y COM",diag_filename)

        # Width of domain, from inputs file.

        dir = os.path.dirname(diag_filename)

        inputs_filename = dir + '/' + wdmerger.get_inputs_filename(dir)

        problo = wdmerger.get_inputs_var(inputs_filename, "geometry.prob_lo")
        probhi = wdmerger.get_inputs_var(inputs_filename, "geometry.prob_hi")

        half_width = (probhi - problo) / 2.0

        # If we are in the rotating frame, we want to transform back to the inertial frame.

        do_rotation = wdmerger.get_inputs_var(inputs_filename, "castro.do_rotation")

        if (do_rotation == 1):
            time = get_column("TIME",diag_filename)
            rot_period = wdmerger.get_inputs_var(inputs_filename, "castro.rotational_period")

            theta = (2.0 * np.pi / rot_period) * time        

            # Rotate from rotating frame to inertial frame using rotation matrix

            xp_old = xp
            yp_old = yp

            xp = xp_old * np.cos(theta) - yp_old * np.sin(theta)
            yp = xp_old * np.sin(theta) + yp_old * np.cos(theta)

            xs_old = xs
            ys_old = ys

            xs = xs_old * np.cos(theta) - ys_old * np.sin(theta)
            ys = xs_old * np.sin(theta) + ys_old * np.cos(theta)

        # Normalize the locations by this domain size, so that they are in the range [-0.5, 0.5]

        xp = xp / half_width[0]
        yp = yp / half_width[1]
        xs = xs / half_width[0]
        ys = ys / half_width[1]

        xlabel = "x"
        ylabel = "y"

        if (N > 1):
            curr_ax = ax[row,col]
        else:
            curr_ax = ax

        # Determine the label and position of this subplot based on the number of inputs.

        curr_ax.plot(xp, yp, 'b--')
        curr_ax.plot(xs, ys, 'r')

        curr_ax.set_xlabel(xlabel, fontsize=20)
        curr_ax.set_ylabel(ylabel, fontsize=20)

        curr_ax.set_xlim([-0.55, 0.55])
        curr_ax.set_ylim([-0.55, 0.55])

        curr_ax.xaxis.set_ticks([-0.50, 0.0, 0.50])
        curr_ax.yaxis.set_ticks([-0.50, 0.0, 0.50])

        curr_ax.set_aspect('auto')

        # The padding ensures that the lower-left ticks on the x- and y-axes don't overlap.

        curr_ax.tick_params(labelsize=16, pad=10)

        # Give it a letter label, for multipanel figures.

        if (N > 1):

            curr_ax.annotate(labels[j], xy=(0.15, 0.85), xycoords='axes fraction', fontsize=16,
                             horizontalalignment='right', verticalalignment='bottom')

        col += 1
        
        if (col > nCols-1):
            col = 0
            row += 1

        j += 1

    # Tighten up figure layout.

    fig.tight_layout()

    # Save it into our designated file, which is usually EPS format.

    plt.savefig(output_filename)

    # Insert git commit hashes into this file from the various code sources.

    wdmerger.insert_commits_into_eps(output_filename, diag_filename, 'diag')

    plt.close()


