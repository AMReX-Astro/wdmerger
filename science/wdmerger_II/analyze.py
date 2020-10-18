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

def get_ni56(results_dir, prefix = ''):
    """Return a list of maximum 56Ni production from all completed sub-directories in results_dir."""

    import numpy as np

    ni56_arr = []

    dir_list = wdmerger.get_parameter_list(results_dir)

    # Strip out those directories that don't match the prefix.

    if (prefix != ''):
        dir_list = [dir for dir in dir_list if dir[0:len(prefix)] == prefix]

    if (dir_list == []):
        return ni56_arr

    diag_filename_list = [results_dir + '/' + directory + '/output/species_diag.out' for directory in dir_list]

    # Ensure all of the output files are there.

    for diag_file in diag_filename_list:
        if not os.path.isfile(diag_file):
            return ni56_arr

    ni56_arr = [np.amax(wdmerger.get_column('Mass ni56', diag_filename)) - np.amin(wdmerger.get_column('Mass ni56', diag_filename))
                for diag_filename in diag_filename_list]

    return ni56_arr



def get_abar(results_dir, prefix = ''):
    """Return a list of maximum 56Ni production from all completed sub-directories in results_dir."""

    import numpy as np

    abar_arr = []

    dir_list = wdmerger.get_parameter_list(results_dir)

    # Strip out those directories that don't match the prefix.

    if (prefix != ''):
        dir_list = [dir for dir in dir_list if dir[0:len(prefix)] == prefix]

    if (dir_list == []):
        return abar_arr

    diag_filename_list = [results_dir + '/' + directory + '/output/species_diag.out' for directory in dir_list]

    # Ensure all of the output files are there.

    for diag_file in diag_filename_list:
        if not os.path.isfile(diag_file):
            return abar_arr

    abar_arr = []

    for diag_filename in diag_filename_list:

        col_names, col_data = wdmerger.get_column_data(diag_filename)

        # The only thing we need is the last row, and we can
        # skip the first two columns (timestep, simulation time).

        spec_names = col_names[2:]
        spec_data = col_data[-1, 2:]

        # Normalize the species by the total mass.

        spec_data /= sum(spec_data)

        # Get the atomic masses of each element by stripping out
        # only the digits from each column header.

        spec_masses = [int(''.join(c for c in s if c.isdigit())) for s in spec_names]

        # Now compute abar.

        abar = sum(spec_masses * spec_data)

        abar_arr.append(abar)

    return abar_arr



def get_e_release(results_dir):
    """Return a list of energy generation from all completed sub-directories in results_dir."""

    import wdmerger_grid_diag_analysis as grid_diag
    import numpy as np

    dir_list = wdmerger.get_parameter_list(results_dir)

    if (dir_list == []):
        return

    diag_filename_list = [results_dir + '/' + directory + '/output/grid_diag.out' for directory in dir_list]

    e_arr = [np.amax(wdmerger.get_column('TOTAL ENERGY', diag_filename)) - np.amin(wdmerger.get_column('TOTAL ENERGY', diag_filename))
             for diag_filename in diag_filename_list]

    # Divide by 10**51

    e_arr = [e / 1.0e51 for e in e_arr]

    return e_arr



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



def write_ni56_table(results_dir, out_filename, v_list, n_list, comment, col1title, label, e_list = [], precision = 3, title = '', caption = ''):
    """Write out a LaTeX table of 56Ni mass as a function of a given variable.

    Optionally we can also include the energy release from the explosion.
    """

    if e_list:
        have_energy = True
    else:
        have_energy = False

    wdmerger.insert_commits_into_txt(out_filename, get_diagfile(results_dir), 'diag')

    out_file = open(out_filename, 'a')

    out_file.write(comment)

    if have_energy:
        col_just = 'ccc'
    else:
        col_just = 'cc'

    out_file.write('\\begin{deluxetable}{' + col_just + '}' + '\n')

    out_file.write('  \\tablecaption{ ' + title + ' \label{table:' + label + '} }' + '\n')

    # Pad the table columns so that they fix the text column. This is an experimentally derived solution.

    padlabel = '\\tbl' + label + 'pad'
    if have_energy:
        padwidth = '0.03in'
    else:
        padwidth = '0.325in'

    out_file.write('  \\newlength{' + padlabel + '}' + '\n')
    out_file.write('  \\setlength{' + padlabel + '}{' + padwidth + '}' + '\n')

    padstr = '\\hspace{\\tbl' + label + 'pad}'

    out_file.write('  \\tablehead{ \\colhead{ ' + padstr + col1title + padstr + ' } & \\colhead{ ' + padstr + 'Max. $^{56}$Ni (${M_{\odot}}$)' + padstr + ' } ')

    if have_energy:
        out_file.write('& \\colhead{ ' + padstr + 'Energy Release (' + '$10^{51}$' + ' erg)' + padstr + ' } ')

    out_file.write('}' + '\n')

    out_file.write('  \\startdata' + '\n')

    for i, v in enumerate(v_list):

        n = n_list[i]

        line_str = '  ' + padstr + str(v) + padstr + ' & ' + padstr + "%.*f" % (precision, float(n)) + padstr

        if have_energy:
            e = e_list[i]
            line_str += ' & ' + padstr + "%.3f" % float(e) + padstr

        if (v != v_list[-1]):
            line_str += " \\\\"  + "\n"
        else:
            line_str += "\n"

        out_file.write(line_str)

    out_file.write('  \\enddata' + '\n')
    out_file.write('  \\tablecomments{ ' + caption + ' }' + '\n')
    out_file.write('\end{deluxetable}' + '\n')
    out_file.close()

    # Write a version of the table using the standard LaTeX tabular environment.

    out_filename += ".nodeluxe"

    wdmerger.insert_commits_into_txt(out_filename, get_diagfile(results_dir), 'diag')

    out_file = open(out_filename, 'a')

    out_file.write(comment)

    if have_energy:
        col_just = '|c|c|c|'
    else:
        col_just = '|c|c|'

    out_file.write('\\begin{table}[h]' + '\n')
    out_file.write('  \\centering' + '\n')
    out_file.write('  \\begin{tabular}{' + col_just + '}' + '\n')
    out_file.write('    \\hline' + '\n')
    out_file.write('    ' + col1title + ' & ' + 'Max. $^{56}$Ni (${M_{\odot}}$)')

    if have_energy:
        out_file.write(' & ' + 'Energy Release ($10^{51}$ erg)')

    out_file.write(' \\\\' + '\n')
    out_file.write('    \\hline' + '\n')

    for i, v in enumerate(v_list):

        n = n_list[i]

        line_str = '    ' + str(v) + ' & ' + "%.3f" % float(n)

        if have_energy:
            e = e_list[i]
            line_str += ' & ' + "%.3f" % float(e)

        line_str += " \\\\"  + "\n"

        out_file.write(line_str)
        out_file.write('    \\hline' + '\n')

    out_file.write('  \\end{tabular}' + '\n')
    out_file.write('  \\caption[' + title + ']{' + caption + ' \label{table:' + label + '}}' + '\n')
    out_file.write('\end{table}' + '\n')
    out_file.close()



def amr_nickel(eps_filename, results_base):
    """Plot the effect of the refinement based on the nuclear burning rate on the nickel generation."""

    import os

    if os.path.isfile(eps_filename):
        return
    else:
        print("Generating file %s" % eps_filename)

    import math
    import numpy as np
    from matplotlib import pyplot as plt

    if not os.path.isdir(results_base):
        return

    diag_file = None

    # First, generate a plot of nickel production.

    ncell_list = sorted([int(n[1:]) for n in os.listdir(results_base)])

    i = 0

    for ncell in ncell_list:

        results_dir = results_base + 'n' + str(ncell) + '/self-heat/dxnuc/'

        if not os.path.isdir(results_dir):
            continue

        # Get the list of parameter values we have tried

        r_list = wdmerger.get_parameter_list(results_dir)

        # Strip out the non-refinement directories

        r_list = [r for r in r_list if r[0] == 'r']

        if (r_list == []):
            continue

        r_list = [float(r[1:]) for r in r_list]

        ni56_arr = get_ni56(results_dir, 'r')

        if ni56_arr == []:
            continue

        if diag_file == None and get_diagfile(results_dir) != None:
            diag_file = get_diagfile(results_dir)

        # Sort the lists

        ni56_arr = np.array([x for _,x in sorted(zip(r_list,ni56_arr))])
        r_list = np.array(sorted(r_list))

        zonesPerDim = [ncell * int(r) for r in r_list]

        plt.plot(zonesPerDim, ni56_arr, markers[i], markersize=12, linestyle=linestyles[i], lw = 4.0, label = 'n = ' + str(ncell))

        i += 1

    plt.tick_params(axis='both', which='major', pad=10)
    plt.xscale('log', basex=2)
    plt.ylim([0.0, 0.6])
    plt.xlabel(r"Effective zones/dimension", fontsize=24)
    plt.ylabel(r"$^{56}$Ni generated (M$_{\odot}$)", fontsize=24)
    plt.tick_params(labelsize=16)
    plt.legend(loc='best', prop={'size':20})
    plt.tight_layout()
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, diag_file, 'diag')
    png_filename = eps_filename[0:-3] + "png"
    plt.savefig(png_filename)

    plt.close()



def amr_detonation(eps_filename, results_base):
    """Plot the effect of the refinement based on the nuclear burning rate on the detonation conditions."""

    import os

    if os.path.isfile(eps_filename):
        return
    else:
        print("Generating file %s" % eps_filename)

    import math
    import numpy as np
    from matplotlib import pyplot as plt

    if not os.path.isdir(results_base):
        return

    diag_file = None

    # First, generate a plot of nickel production.

    ncell_list = sorted([int(n[1:]) for n in os.listdir(results_base)])

    # Now, generate a plot of the distance from the center of the ignition point.
    
    i = 0

    for ncell in ncell_list:

        results_dir = results_base + 'n' + str(ncell) + '/self-heat/dxnuc/'

        if not os.path.isdir(results_dir):
            continue

        # Get the list of parameter values we have tried

        r_list = wdmerger.get_parameter_list(results_dir)

        # Strip out the non-refinement directories

        r_list = [r for r in r_list if r[0] == 'r']

        if (r_list == []):
            continue

        # Sort them numerically.

        r_list = ['r' + str(r) for r in sorted([int(r[1:]) for r in r_list])]

        # Cycle through the plot files.
        # Get the maximum value of t_sound_t_enuc, and its distance from the origin.
        # Stop if we hit a detonation (t_sound_t_enuc > 1).

        for r in r_list:
            results_dir = results_base + 'n' + str(ncell) + '/self-heat/dxnuc/' + r

            prefix = 'smallplt'
            plt_list = [results_dir + '/output/' + plt for plt in wdmerger.get_plotfiles(results_dir, prefix = prefix)]

            for plt in plt_list:
                [ts_te_max, x, y, z] = wdmerger.get_maxloc(plt, 't_sound_t_enuc')

                if ts_te_max > 1.0:
                    [rho_max, x, y, z] = wdmerger.get_maxloc(plt, 'density')
                    dist = np.sqrt(x**2 + y**2)
                    time = wdmerger.get_time_from_plotfile(plt)
                    print(ncell, r, time, rho_max, ts_te_max, dist.v / 1.e5)

                    break

        i += 1


def rho_T_sliceplots_eps_rename(pltfile, output_dir, results_dir):

    eps_file = output_dir + "/rho_T_slice" + "_t_" + str("%.2f" % wdmerger.get_time_from_plotfile(results_dir + '/output/' + pltfile)) + '.eps'

    return eps_file



def rho_T_sliceplots_doit(inputs): #output_dir, results_dir, domain_frac, x_ticks, y_ticks, scale_exp):

    pltfile     = inputs[0]
    output_dir  = inputs[1]
    results_dir = inputs[2]
    domain_frac = inputs[3]
    x_ticks     = inputs[4]
    y_ticks     = inputs[5]
    scale_exp   = inputs[6]

    eps_file = rho_T_sliceplots_eps_rename(pltfile, output_dir, results_dir)

    if os.path.isfile(eps_file) and os.path.isfile(eps_file.replace('eps','png')):
        return

    print("Generating plot with filename " + eps_file)

    wdmerger.rho_T_sliceplot(eps_file, results_dir + '/output/' + pltfile,
                             domain_frac = domain_frac,
                             x_ticks = x_ticks,
                             y_ticks = y_ticks,
                             scale_exp = scale_exp)



def rho_T_sliceplots(output_dir, results_dir, smallplt = True, domain_frac = 1.0,
                     x_ticks = [2.0e9, 4.0e9], y_ticks = [2.0e9, 4.0e9], scale_exp = 9):
    """Create a rho/T sliceplot for every plotfile in a given directory."""

    import os
    from multiprocessing import Pool

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    from PIL import Image

    if (smallplt):
        plt_prefix = 'smallplt'

    plt_list = wdmerger.get_plotfiles(results_dir, plt_prefix)
    inputs_list = [[pltfile, output_dir, results_dir, domain_frac, x_ticks, y_ticks, scale_exp] for pltfile in plt_list]

    with Pool() as p:
        p.map(rho_T_sliceplots_doit, inputs_list)

    eps_list = [rho_T_sliceplots_eps_rename(pltfile, output_dir, results_dir) for pltfile in plt_list]
    jpg_list = [eps.replace('eps', 'jpg') for eps in eps_list]

    mpg_filename = output_dir + "/rho_T_slice.mpg"

    if not os.path.isfile(mpg_filename):
        print("Generating file %s" % mpg_filename)
        wdmerger.make_movie(output_dir, jpg_list, mpg_filename)



#
# Calculate the collision time.
#

def collision_time(results_dir):
    """Calculate the time of the collision."""

    out_file = results_dir + '/collision_time.txt'

    if os.path.isfile(out_file):
        return

    print("Calculating collision time in directory " + results_dir)

    import yt

    # Load all the output files in the output directory.
    # Search them one by one (in reverse) until we find
    # the time where the density threshold is crossed
    # at the center of the simulation.

    density_threshold = 1.0e2

    plots = wdmerger.get_plotfiles(results_dir, prefix = 'smallplt')
    plots.reverse()

    if (plots == []) or (plots == None):
        return

    for i, plot in enumerate(plots):

        ds = yt.load(results_dir + '/output/' + plot)

        # Get a region that is only one coarse zone-width wide

        dx = ds.domain_width / ds.domain_dimensions

        sph = ds.sphere([0.0, 0.0, 0.0], max(dx))

        # Find the maximum density in that region

        rho_max = sph.max("density")

        if rho_max < density_threshold:
            idx = i - 1
            break

    ds = yt.load(results_dir + '/output/' + plots[idx])

    file = open(out_file, 'w')

    file.write(str(ds.current_time.v))

    file.close()



# Main execution: do all of the analysis routines.

if __name__ == "__main__":
    """Generate the plots and tables for the collision tests."""

    import os

    mass = '0.64'
    he = '0.16'

    results_base = 'results/m' + mass + '/he' + he + '/'
    plots_dir = 'plots/'

    burning_r_list = os.listdir(results_base)

    for burning_r in sorted(burning_r_list):

        hydro_r_list = os.listdir(results_base + burning_r)

        for hydro_r in sorted(hydro_r_list):

            results_dir = results_base + burning_r + '/' + hydro_r + '/'
            plot_dir = plots_dir + 'slices/' + burning_r + '/' + hydro_r + '/'

            rho_T_sliceplots(plot_dir, results_dir)
