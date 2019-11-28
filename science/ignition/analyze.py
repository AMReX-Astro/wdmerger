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

    # Skip this if we've already generated the plots.

    if os.path.isfile(file_base + '.eps'):
        return

    size = 3.2768e9

    ncell_list = sorted([int(n[1:]) for n in os.listdir(results_base)])

    res_lists = [[]] * len(ncell_list)
    dist_lists = [[]] * len(ncell_list)
    time_lists = [[]] * len(ncell_list)
    ts_te_lists = [[]] * len(ncell_list)

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
        ts_te_list = []

        for r in r_list:
            results_dir = results_base + '/n' + str(ncell) + '/' + r

            if not os.path.isdir(results_dir):
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

                    [ts_te_max, _, _, _] = wdmerger.get_maxloc(plot, 't_sound_t_enuc')

                    ts_te_list.append(ts_te_max)

                    break

                else:
                    T_last = T_max

            res_lists[i] = res_list
            dist_lists[i] = dist_list
            time_lists[i] = time_list
            ts_te_lists[i] = ts_te_list

    # Plot location of the ignition and time on the same axis.
    # We want to combine data at low resolution using their coarse
    # grid runs, with the AMR data at the highest resolution.

    res_list = []
    dist_list = []
    time_list = []
    ts_te_list = []

    for i in range(len(ncell_list)-1):
        res_list.append(res_lists[i][0])
        dist_list.append(dist_lists[i][0])
        time_list.append(time_lists[i][0])
        ts_te_list.append(ts_te_lists[i][0])

    for res in res_lists[-1]:
        res_list.append(res)

    for dist in dist_lists[-1]:
        dist_list.append(dist)

    for time in time_lists[-1]:
        time_list.append(time)

    for ts_te in ts_te_lists[-1]:
        ts_te_list.append(ts_te)

    fig, ax1 = plt.subplots()
    ax1.plot(res_list, dist_list, linestyle=linestyles[0], marker=markers[0], color=colors[0], lw=2.0, label='Ignition location')

    ax1.set_xscale('log', basex=10)
    ax1.set_xlabel(r"Finest resolution (km)", fontsize=24)
    ax1.set_ylabel(r"Ignition location (km)", fontsize=24, color=colors[0])
    ax1.tick_params(axis='y', labelcolor=colors[0], labelsize=20)
    ax1.tick_params(axis='x', labelsize=20)
    ax1.set_ylim(bottom=0.0)
    if 'co' in file_base:
        ax1.set_ylim(top=6000.0)
    else:
        ax1.set_ylim(top=1100.0)

    ax2 = ax1.twinx()
    ax2.plot(res_list, time_list, linestyle=linestyles[1], marker=markers[1], color=colors[1], lw=2.0, label='Ignition time')
    ax2.set_ylabel(r"Ignition time (s)", fontsize=24, color=colors[1])
    ax2.tick_params(axis='y', labelcolor=colors[1], labelsize=20)
    ax2.set_ylim(bottom=0.0)
    if 'co' in file_base:
        ax2.set_ylim(top=4.0)
    else:
        ax2.set_ylim(top=1.1)

    fig.set_size_inches(11, 8.5)
    fig.tight_layout()
    plt.savefig(file_base + '.eps')
    plt.savefig(file_base + '.png')

    plt.close()

    fig, ax1 = plt.subplots()

    ax1.plot(res_list, ts_te_list, linestyle=linestyles[0], marker=markers[0], color=colors[0], lw=2.0)

    ax1.set_xscale('log', basex=10)
    ax1.set_xlabel(r"Finest resolution (km)", fontsize=24)
    ax1.set_ylabel(r"$\tau_{\rm s}\, /\, \tau_{\rm e}$", fontsize=24)
    ax1.tick_params(axis='y', labelsize=20)
    ax1.tick_params(axis='x', labelsize=20)

    fig.set_size_inches(11, 8.5)
    fig.tight_layout()
    plt.savefig(file_base + '_ts_te.eps')
    plt.savefig(file_base + '_ts_te.png')

    plt.close()



# Main execution: do all of the analysis routines.

if __name__ == "__main__":
    """Generate the plots and tables for the collision tests."""

    import os
        
    plots_dir = 'plots/'

    # First do the tests with helium.

    results_base = 'results/he_co'
    file_base = 'amr_ignition_he_co'

    run_dir = results_base
    run_str = file_base

    amr_ignition(plots_dir + run_str, run_dir)

    # Now do the pure C/O plots.

    results_base = 'results/co/'
    file_base = 'amr_ignition_co'

    run_dir = results_base
    run_str = file_base

    amr_ignition(plots_dir + run_str, run_dir)
