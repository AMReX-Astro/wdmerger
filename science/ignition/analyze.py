import numpy as np
import wdmerger
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import os
import math
import yt

linestyles = ['-', '--', ':', '-.']
markers = ['', '', '', 'o']

# Colors are chosen with the assistance of colorbrewer2.org,
# using 4 qualitative classes.

colors = ['#1f78b4', '#33a02c', '#a6cee3', '#b2df8a']

def amr_ignition(file_base, results_base, do_amr = True):
    """Plot the effect of the resolution on the detonation conditions."""

    if not os.path.isdir(results_base):
        return

    size = 1.6384e9

    ncell_list = sorted([int(n[1:]) for n in os.listdir(results_base)])

    res_lists = [[]] * len(ncell_list)
    dist_lists = [[]] * len(ncell_list)
    time_lists = [[]] * len(ncell_list)

    base_res_list = []

    for i, ncell in enumerate(ncell_list):

        base_res_list.append(size / ncell / 1.e5)

        dir = results_base + '/n' + str(ncell)

        if not os.path.isdir(dir):
            continue

        # Generate a plot of the distance from the center of the ignition point,
        # and also the time of the ignition.

        # Get the list of refinement values we have tried.

        r_list = wdmerger.get_parameter_list(dir)

        if (r_list == []):
            continue

        # Sort them numerically.

        r_list = ['r' + str(r) for r in sorted([int(r[1:]) for r in r_list])]

        # Cycle through the plot files.
        # Get the distance from the origin of the point where T > 4e9 K, if it exists.

        zonesPerDim = []
        res_list = []
        dist_list = []
        time_list = []

        for r in r_list:
            results_dir = results_base + '/n' + str(ncell) + '/' + r

            if not os.path.isdir(results_dir):
                continue

            # Only do AMR plots for certain resolutions, to avoid the
            # plot being too crowded.

            if (int(r[1:]) > 1) and (ncell not in [128, 4096, 16384]):
                continue

            print('Searching in directory ' + results_dir)

            prefix = 'det_x_plt'
            plot_list = [results_dir + '/output/' + plot for plot in wdmerger.get_plotfiles(results_dir, prefix = prefix)]

            # Since we'll have stopped after the ignition has occurred,
            # we'll search through the plotfiles backward.

            T_max = -1.0
            T_last = -1.0

            for plot in reversed(plot_list):

                [T_max, x, y, z] = wdmerger.get_maxloc(plot, 'Temp')

                if (T_max < 4.0e9 and T_last >= 4.0e9) or (T_max >= 4.0e9 and T_last < 0.0):
                    dist = abs(x.v) / 1.0e5
                    time = wdmerger.get_time_from_plotfile(plot)

                    zonesPerDim.append(ncell * int(r[1:]))

                    ds = yt.load(plot)
                    problo = ds.domain_left_edge.v
                    probhi = ds.domain_right_edge.v
                    size = probhi[0] - problo[0]

                    res_list.append(size / (ncell * int(r[1:])) / 1.0e5)
                    dist_list.append(dist)
                    time_list.append(time)

                    break

                else:
                    T_last = T_max

            res_lists[i] = res_list
            dist_lists[i] = dist_list
            time_lists[i] = time_list

    has_plot = False

    # First add all plots with refinement.

    if do_amr:
        plt_count = 0
        for i, (res_list, dist_list) in enumerate(zip(res_lists, dist_lists)):
            if len(res_list) > 1:
                has_plot = True
                lbl = '{} km base + AMR'.format(int(base_res_list[i]))
                plt.plot(res_list, dist_list, linestyle=linestyles[plt_count], marker=markers[plt_count], lw=2.0, label=lbl)
                plt_count += 1

    # Now add the base refinement cases.

    coarse_list = []
    r_list = []
    for i, ncell in enumerate(ncell_list):
        if len(dist_lists[i]) > 0:
            has_plot = True
            r_list.append(res_lists[i][0])
            coarse_list.append(dist_lists[i][0])

    if len(coarse_list) > 0:
        lbl = 'Uniform grid'
        plt.plot(r_list, coarse_list, color='k', lw=2.0, label=lbl, marker='o', markersize=12)

    if has_plot:
        plt.tick_params(axis='both', which='major', pad=10, labelsize=16)
        plt.xscale('log', basex=10)
        plt.xlabel(r"Finest resolution (km)", fontsize=24)
        plt.ylabel(r"Ignition location (km)", fontsize=24)
        plt.legend(loc='best', prop={'size':12}, fontsize=12)
        plt.rcParams["figure.figsize"] = (11, 8.5)
        plt.tight_layout()
        plt.savefig(file_base + '.eps')
        plt.savefig(file_base + '.png')

        plt.close()



# Main execution: do all of the analysis routines.

if __name__ == "__main__":
    """Generate the plots and tables for the collision tests."""

    import os
        
    plots_dir = 'plots/'

    file_base = 'amr_ignition'
    results_base = 'results/'

    for burning_mode in ['self-heat', 'suppressed']:

        run_dir = '/'.join([results_base, burning_mode])
        run_str = '_'.join([file_base, burning_mode])

        if burning_mode == 'self-heat':
            do_amr = True
        else:
            do_amr = False

        amr_ignition(plots_dir + run_str, run_dir, do_amr)
