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



def amr_ignition(eps_filename, results_base):
    """Plot the effect of the refinement based on the nuclear burning rate on the detonation conditions."""

    import os

    if os.path.isfile(eps_filename):
        return
    else:
        print("Generating file %s" % eps_filename)

    import math
    import numpy as np
    from matplotlib import pyplot as plt
    import yt

    if not os.path.isdir(results_base):
        return

    # Loop over resolutions of the temperature refinement region.

    temperr_r_list = [1, 4, 16, 64, 256, 1024, 4096, 16384]

    size = 8.192e8

    for temperr_r in temperr_r_list:

        filename = eps_filename.replace('.eps', '_temperr_r_' + str(temperr_r) + '_location.eps')

        if os.path.isfile(filename):
            continue

        ncell_list = sorted([int(n[1:]) for n in os.listdir(results_base)])

        # Generate a plot of the distance from the center of the ignition point,
        # and also the time of the ignition.

        res_lists = [[]] * len(ncell_list)
        dist_lists = [[]] * len(ncell_list)
        time_lists = [[]] * len(ncell_list)

        base_res_list = []

        for i, ncell in enumerate(ncell_list):

            base_res_list.append(size / ncell / 1.e5)

            results_dir = results_base + '/n' + str(ncell)

            if not os.path.isdir(results_dir):
                continue

            # Get the list of parameter values we have tried

            r_list = wdmerger.get_parameter_list(results_dir, include_incomplete = True)

            # Strip out the non-refinement directories

            r_list = [r for r in r_list if r[0] == 'r']

            if (r_list == []):
                continue

            # Sort them numerically.

            r_list = ['r' + str(r) for r in sorted([int(r[1:]) for r in r_list])]

            # Cycle through the plot files.
            # Get the maximum value of t_sound_t_enuc, and its distance from the origin.
            # Stop if we hit an ignition (T > 5e9).

            zonesPerDim = []
            res_list = []
            dist_list = []
            time_list = []

            for r in r_list:
                results_dir = results_base + '/n' + str(ncell) + '/' + r + '/temperr1.0d8/temperr_r' + str(temperr_r) + '/nbuf4'

                if not os.path.isdir(results_dir):
                    continue

                print('Searching in directory ' + results_dir)

                prefix = 'det_x_plt'
                plot_list = [results_dir + '/output/' + plot for plot in wdmerger.get_plotfiles(results_dir, prefix = prefix)]

                # Assume the last output is the one where we cross the temperature threshold.

                if len(plot_list) > 1:
                    plot = plot_list[-1]
                elif len(plot_list) == 1:
                    plot = plot_list[0]
                else:
                    continue


                [T_max, x, y, z] = wdmerger.get_maxloc(plot, 'Temp')

                if T_max > 5.e9:
                    dist = abs(x.v) / 1.e5
                    time = wdmerger.get_time_from_plotfile(plot)

                    zonesPerDim.append(ncell * int(r[1:]))

                    ds = yt.load(plot)
                    problo = ds.domain_left_edge.v
                    probhi = ds.domain_right_edge.v
                    size = probhi[0] - problo[0]

                    res_list.append(size / (ncell * int(r[1:])) / 1.e5)
                    dist_list.append(dist)
                    time_list.append(time)

            res_lists[i] = res_list
            dist_lists[i] = dist_list
            time_lists[i] = time_list

        for i, (res_list, dist_list) in enumerate(zip(res_lists, dist_lists)):
            if len(res_list) > 1:
                lbl = '%3.2e' % base_res_list[i]
                plt.plot(res_list, dist_list, markersize=12, lw=1.0, label=lbl)

        plt.tick_params(axis='both', which='major', pad=10, labelsize=16)
        plt.xscale('log', basex=10)
        plt.xlabel(r"Effective finest resolution (km)", fontsize=24)
        plt.ylabel(r"Ignition location (km)", fontsize=24)
        plt.legend(loc='best', prop={'size':12}, ncol=2, title='Effective base resolution (km)', fontsize=12)
        plt.rcParams["figure.figsize"] = (11, 8.5)
        plt.tight_layout()
        plt.savefig(filename)
        png_filename = filename[0:-3] + "png"
        plt.savefig(png_filename)

        plt.close()

        filename = filename.replace('location', 'time')

        for i, (res_list, time_list) in enumerate(zip(res_lists, time_lists)):
            print(i, res_list, time_list)
            if len(res_list) > 1:
                lbl = '%3.2e' % base_res_list[i]
                plt.plot(res_list, time_list, markersize=12, lw=1.0, label=lbl)

        plt.tick_params(axis='both', which='major', pad=10, labelsize=16)
        plt.xscale('log', basex=10)
        plt.xlabel(r"Effective finest resolution (km)", fontsize=24)
        plt.ylabel(r"Ignition time (s)", fontsize=24)
        plt.legend(loc='best', prop={'size':12}, ncol=2, title='Effective base resolution (km)', fontsize=12)
        plt.rcParams["figure.figsize"] = (11, 8.5)
        plt.tight_layout()
        plt.savefig(filename)
        png_filename = filename[0:-3] + "png"
        plt.savefig(png_filename)

        plt.close()



# Main execution: do all of the analysis routines.

if __name__ == "__main__":
    """Generate the plots and tables for the collision tests."""

    import os
        
    plots_dir = 'plots/'

    size = '8.192e8'
    burning_mode = 'self-heat'

    file_base = burning_mode
    results_base = 'results/' + 'size' + size + '/' + burning_mode

    dtnuc_list = os.listdir(results_base)

    for dtnuc in dtnuc_list:

        tempgrad_list = os.listdir(results_base + '/' + dtnuc)

        for tempgrad in tempgrad_list:

            dens_list = os.listdir(results_base + '/' + dtnuc + '/' + tempgrad)

            for dens in dens_list:

                g_list = os.listdir(results_base + '/' + dtnuc + '/' + tempgrad + '/' + dens)

                for g in g_list:

                    v_list = os.listdir(results_base + '/' + dtnuc + '/' + tempgrad + '/' + dens + '/' + g + '/')

                    for v in v_list:

                        file_string = file_base + '_dtnuc_' + dtnuc[5:] + '_tempgrad_' + tempgrad[8:] + '_dens_' + dens[4:] + '_g_' + g[1:] + '_v_' + v[1:]
                        results_dir = results_base + '/' + dtnuc + '/' + tempgrad + '/' + dens + '/' + g + '/' + v

                        ncell_list = os.listdir(results_dir)

                        amr_ignition(plots_dir + 'amr_ignition_' + file_string + '.eps', results_dir)
