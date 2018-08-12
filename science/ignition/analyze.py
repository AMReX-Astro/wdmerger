import numpy as np
import wdmerger
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import os

linestyles = ['-', '--', ':', '-.']
markers = ['', '', '', 'o']

# Colors are chosen with the assistance of colorbrewer2.org,
# using 4 qualitative classes.

colors = ['#1f78b4', '#33a02c', '#a6cee3', '#b2df8a']

def get_diagfile(results_dir):
    """Return the name of a diagnostics file in results_dir that is completed.

       Useful for the functions that insert git commits into files and need a 
       representative diagnostic file to get the git commits from. We arbitrarily
       pick the last complete directory."""

    dir_list = wdmerger.get_parameter_list(results_dir)

    if (dir_list == []):
        return

    dir_list = [results_dir + '/' + directory for directory in dir_list]

    for dir in dir_list:
        if (wdmerger.is_dir_done(dir)):
            return dir + '/output/species_diag.out'



def amr_ignition(filename_base, results_base):
    """Plot the effect of the refinement based on the nuclear burning rate on the detonation conditions."""

    import os

    import math
    import numpy as np
    from matplotlib import pyplot as plt
    import yt

    if not os.path.isdir(results_base):
        return

    size = 8.192e8

    tempgrad_list = []
    dxnuc_list = []

    tempgrad_r_list = []
    dxnuc_r_list = []

    # First we need to loop through the directories and figure out
    # what values we have for dxnuc and tempgrad, and their refinements.
    # Then we can reverse the order and loop through those, and make plots
    # as a function of ncell.

    ncell_list = sorted([int(n[1:]) for n in os.listdir(results_base)])

    for i, ncell in enumerate(ncell_list):

        tempgrad_list_loc = [tempgrad[8:] for tempgrad in os.listdir(results_base + '/n' + str(ncell))]

        for tempgrad in tempgrad_list_loc:

            if tempgrad not in tempgrad_list:
                tempgrad_list.append(tempgrad)

            dxnuc_list_loc = [dxnuc[5:] for dxnuc in os.listdir(results_base + '/n' + str(ncell) + '/tempgrad' + str(tempgrad))]

            for dxnuc in dxnuc_list_loc:

                if dxnuc not in dxnuc_list:
                    dxnuc_list.append(dxnuc)

                tempgrad_r_list_loc = [tempgrad_r[10:] for tempgrad_r in os.listdir(results_base + '/n' + str(ncell) + '/tempgrad' + str(tempgrad) + '/dxnuc' + dxnuc)]

                for tempgrad_r in tempgrad_r_list_loc:

                    if tempgrad_r not in tempgrad_r_list:
                        tempgrad_r_list.append(tempgrad_r)

                    dxnuc_r_list_loc = [dxnuc_r[7:] for dxnuc_r in os.listdir(results_base + '/n' + str(ncell) + '/tempgrad' + str(tempgrad) + '/dxnuc' + dxnuc + '/tempgrad_r' + tempgrad_r)]

    # Now we can do the actual plot generation.

    for tempgrad in tempgrad_list:
        for dxnuc in dxnuc_list:
            for tempgrad_r in tempgrad_r_list:

                filename = filename_base + '_dxnuc_' + dxnuc + '_tempgrad_' + tempgrad + '_tempgrad_r' + tempgrad_r + '.eps'

                if os.path.isfile(filename):
                    continue

                res_lists = [[]] * len(ncell_list)
                dist_lists = [[]] * len(ncell_list)
                time_lists = [[]] * len(ncell_list)

                base_res_list = []

                for i, ncell in enumerate(ncell_list):

                    base_res_list.append(size / ncell / 1.e5)

                    dir = results_base + '/n' + str(ncell) + '/tempgrad' + tempgrad + '/dxnuc' + dxnuc + '/tempgrad_r' + tempgrad_r

                    if not os.path.isdir(dir):
                        continue

                    # Generate a plot of the distance from the center of the ignition point,
                    # and also the time of the ignition.

                    # Get the list of refinement values we have tried

                    r_list = wdmerger.get_parameter_list(dir)

                    # Strip out the non-refinement directories

                    r_list = [r for r in r_list if r[0:7] == 'dxnuc_r']

                    if (r_list == []):
                        continue

                    # Sort them numerically.

                    r_list = ['r' + str(r) for r in sorted([int(r[7:]) for r in r_list])]

                    # Cycle through the plot files.
                    # Get the distance from the origin of the point where T > 4e9 K, if it exists.

                    zonesPerDim = []
                    res_list = []
                    dist_list = []
                    time_list = []

                    for r in r_list:
                        results_dir = results_base + '/n' + str(ncell) + '/tempgrad' + tempgrad + '/dxnuc' + dxnuc + '/tempgrad_r' + tempgrad_r + '/dxnuc_' + r

                        if not os.path.isdir(results_dir):
                            continue

                        print('Searching in directory ' + results_dir)

                        prefix = 'det_x_plt'
                        plot_list = [results_dir + '/output/' + plot for plot in wdmerger.get_plotfiles(results_dir, prefix = prefix)]

                        # Since most of the time we'll have stopped after the ignition has occurred,
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

                                res_list.append(size / (ncell * max(int(r[1:]), int(tempgrad_r))) / 1.0e5)
                                dist_list.append(dist)
                                time_list.append(time)

                                break

                            else:
                                T_last = T_max

                    res_lists[i] = res_list
                    dist_lists[i] = dist_list
                    time_lists[i] = time_list

                    print('res', res_list)
                    print('dist', dist_list)
                    print('time', time_list)

                has_plot = False

                # First generate a plot using the base refinement.

                coarse_list = []
                r_list = []
                for i, ncell in enumerate(ncell_list):
                    if len(dist_lists[i]) > 0:
                        has_plot = True
                        r_list.append(res_lists[i][0])
                        coarse_list.append(dist_lists[i][0])

                if len(coarse_list) > 0:
                    if int(tempgrad_r) == 1:
                        lbl = 'Uniform grid'
                    else:
                        lbl = tempgrad_r + 'x temperature gradient AMR'
                    plt.plot(r_list, coarse_list, color='k', label=lbl, marker='o', markersize=12)

                # Now add all plots with dxnuc AMR.

                for i, (res_list, dist_list) in enumerate(zip(res_lists, dist_lists)):
                    if len(res_list) > 1:
                        has_plot = True
                        lbl = '{} km base + auto-ignition AMR'.format(int(base_res_list[i]))
                        plt.plot(res_list, dist_list, lw=2.0, label=lbl)

                if has_plot:
                    plt.tick_params(axis='both', which='major', pad=10, labelsize=16)
                    plt.xscale('log', basex=10)
                    plt.xlabel(r"Finest resolution (km)", fontsize=24)
                    plt.ylabel(r"Ignition location (km)", fontsize=24)
                    plt.legend(loc='best', prop={'size':12}, fontsize=12)
                    plt.rcParams["figure.figsize"] = (11, 8.5)
                    plt.tight_layout()
                    plt.savefig(filename)
                    plt.savefig(filename.replace('.eps', '.png'))

                    plt.close()



# Main execution: do all of the analysis routines.

if __name__ == "__main__":
    """Generate the plots and tables for the collision tests."""

    import os
        
    plots_dir = 'plots/'

    cfrac = '0.5d0'
    ofrac = '0.0d0'

    file_base = 'amr_ignition'
    results_base = 'results/'

    dens_list = os.listdir(results_base)

    # Only do the ignition plot for the
    # self-heating burn, because that's
    # the only one we actually let
    # get to ignition.
    
    burning_mode = 'self-heat'

    g_list = os.listdir(results_base + '/' + burning_mode)

    for g in g_list:

        v_list = os.listdir(results_base + '/' + burning_mode + '/' + g)

        for v in v_list:

            cfrac_list = os.listdir(results_base + '/' + burning_mode + '/' + g + '/' + v)

            for cfrac in cfrac_list:

                ofrac_list = os.listdir(results_base + '/' + burning_mode + '/' + g + '/' + v + '/' + cfrac)

                for ofrac in ofrac_list:

                    run_dir = results_base + '/' + burning_mode + '/' + g + '/' + v + '/' + cfrac + '/' + ofrac
                    run_str = file_base + '_' + burning_mode + '_' + g + '_' + v + '_' + cfrac + '_' + ofrac

                    amr_ignition(plots_dir + run_str, run_dir)
