import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger
import os

# Problem parameters

ncell_arr = [ 256, 512, 1024 ]
vel_arr   = [ 'motion', 'static' ]

linestyles = ["-", "--","-"]
markers = ["", "", "o"]

nCols = 7

# Effective radii of the star

for vel in vel_arr:

    for c in range(nCols):

        # Set the name of the diagnostic output file

        eps_filename = "plots/single_star_" + str(vel) + "_1e" + str(c) + "_radius.eps"

        if (os.path.isfile(eps_filename)):
            print "Plot with filename " + eps_filename + " already exists; skipping."            
            continue        

        print "Generating plot with filename " + eps_filename

        idx = 0

        for n in ncell_arr:

            # Account for the fact that we run at higher resolution when we're in the moving frame.

            if ( vel == 'motion' ):
                nc = n * 4
            elif ( vel == 'static' ):
                nc = n

            diag_filename = "results/" + str(vel) + "/n" + str(nc) + "/output/star_diag.out"

            timeCol = wdmerger.get_column('TIME', diag_filename)

            # Obtain the spherical radii denoting different density cutoffs

            radCol = wdmerger.get_column('PRIMARY 1E' + str(c) + ' RADIUS', diag_filename)

            # Divide radii by initial radii, then calculate difference with respect to
            # the initial radius.

            radCol /= radCol[0]

            plt.plot(timeCol, radCol, linestyle=linestyles[idx], linewidth=4.0, marker = markers[idx], markevery=500, markersize=8.0, label=r'n = ' + str(nc))

            idx += 1

        plt.legend(loc='upper center',ncol=len(ncell_arr))
        plt.xlabel("Time (s)", fontsize=20)
        plt.ylabel("Radius / Initial Radius", fontsize=20)

        plt.xkcd()

        if ( vel == 'motion' ):
            ylim_l = 0.95
            ylim_r = 2.05
        elif ( vel == 'static' ):
            ylim_l = 0.995
            ylim_r = 1.1

        plt.axis([timeCol[0], timeCol[-1], ylim_l, ylim_r])
        plt.tick_params(labelsize=16)

        plt.savefig(eps_filename)
        wdmerger.insert_commits_into_eps(eps_filename, diag_filename, 'diag')

        plt.clf()



# Kinetic energy of the star (inertial frame only)

vel = 'static'

eps_filename = 'plots/single_star_' + vel + '_ke.eps'

if (os.path.isfile(eps_filename)):
    print "Plot with filename " + eps_filename + " already exists; skipping."            
else:
    print "Generating plot with filename " + eps_filename

    idx = 0

    for n in ncell_arr:

        # Account for the fact that we run at higher resolution when we're in the moving frame.

        if ( vel == 'motion' ):
            nc = n * 4
        elif ( vel == 'static' ):
            nc = n

        diag_filename = "results/" + str(vel) + "/n" + str(nc) + "/output/grid_diag.out"

        timeCol = wdmerger.get_column('TIME', diag_filename)

        # Obtain the spherical radii denoting different density cutoffs

        keCol = wdmerger.get_column('KIN. ENERGY', diag_filename)

        plt.plot(timeCol, keCol, linestyle=linestyles[idx], linewidth=4.0, marker = markers[idx], markevery=500, markersize=8.0, label=r'n = ' + str(nc))

        idx += 1

    plt.legend(loc='upper center',ncol=len(ncell_arr))
    plt.xlabel("Time (s)", fontsize=20)
    plt.ylabel("Kinetic Energy (erg)", fontsize=20)
    plt.yscale("log")

    plt.axis([timeCol[0], timeCol[-1], 1.e43, 1.e47])
    plt.tick_params(labelsize=16)

    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, diag_filename, 'diag')

    plt.clf()


# Now generate plots that compare the reference frames

for c in range(nCols):

    eps_filename = "plots/single_star_compare_1e" + str(c) + "_radius.eps"

    if (os.path.isfile(eps_filename)):
        print "Plot with filename " + eps_filename + " already exists; skipping."
        continue        

    print "Generating plot with filename " + eps_filename

    idx = 0

    for n in ncell_arr:

        for vel in vel_arr:

            # Account for the fact that we run at higher resolution when we're in the moving frame.

            if ( vel == 'motion' ):
                nc = n * 4
            elif ( vel == 'static' ):
                nc = n

            # First do the static case

            diag_filename = "results/" + str(vel) + "/n" + str(nc) + "/output/star_diag.out"

            timeCol = wdmerger.get_column('TIME', diag_filename)

            # Obtain the spherical radii denoting different density cutoffs

            radCol = wdmerger.get_column('PRIMARY 1E' + str(c) + ' RADIUS', diag_filename)

            # Divide radii by initial radii, then calculate difference with respect to
            # the initial radius.

            radCol /= radCol[0]

            plt.plot(timeCol, radCol, linestyle=linestyles[idx], linewidth=4.0, marker=markers[idx], markevery=500, markersize=8.0, label=str(vel) + ", " + str(nc))

        idx += 1

    plt.legend(loc='upper center',ncol=len(ncell_arr),columnspacing=0.5)
    plt.xlabel("Time (s)", fontsize=20)
    plt.ylabel("Radius / Initial Radius", fontsize=20)
    plt.axis([timeCol[0], timeCol[-1], 0.95, 2.05])
    plt.tick_params(labelsize=16)

    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, diag_filename, 'diag')

    plt.clf()
