import numpy as np
import wdmerger
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import os

linestyles = ['-', '--', ':', '-', '--', ':']
markers = ['', '', '', 'o', 's', 'D']
colors = ['b', 'g', 'r', 'c', 'm', 'k']

def get_ni56(results_dir):
    """Return a list of maximum 56Ni production from all completed sub-directories in results_dir."""

    import wdmerger_spec_diag_analysis as spec_diag
    import numpy as np

    dir_list = wdmerger.get_parameter_list(results_dir)

    if (dir_list == []):
        return

    diag_filename_list = [results_dir + '/' + directory + '/output/species_diag.out' for directory in dir_list]

    ni56_arr = [np.amax(wdmerger.get_column('Mass ni56', diag_filename)) - np.amin(wdmerger.get_column('Mass ni56', diag_filename))
                for diag_filename in diag_filename_list]

    return ni56_arr



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



def burning_limiter_e(eps_filename, results_base):
    """Plot the effect of the burning timestep factor castro.dtnuc_e."""

    if (os.path.isfile(eps_filename)):
        return

    results_dir = results_base + 'burning_limiter_e/'

    # Get the list of parameter values we have tried

    dtnuc_list = wdmerger.get_parameter_list(results_dir)

    if (dtnuc_list == []):
        return

    ni56_arr = np.array(get_ni56(results_dir))

    dtnuc_list = np.array([float(dtnuc[len('dt'):]) for dtnuc in dtnuc_list])

    # Sort the lists

    ni56_arr = np.array([x for _,x in sorted(zip(dtnuc_list,ni56_arr))])
    dtnuc_list = sorted(dtnuc_list)

    plt.plot(dtnuc_list, ni56_arr, linestyle=linestyles[0], lw = 4.0)

    plt.xscale('log')
    plt.xlabel(r"Nuclear burning timestep factor $\Delta t_{be}$", fontsize=20)
    plt.ylabel(r"$^{56}$Ni generated (M$_{\odot}$)", fontsize=20)
    plt.tick_params(labelsize=20)
    plt.tight_layout()
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, get_diagfile(results_dir), 'diag')

    plt.close()



def burning_limiter_X(eps_filename, results_base):
    """Plot the effect of the burning timestep limiter castro.dtnuc_X."""

    if (os.path.isfile(eps_filename)):
        return

    results_dir = results_base + 'burning_limiter_X/'

    # Get the list of parameter values we have tried

    dtnuc_list = wdmerger.get_parameter_list(results_dir)

    if (dtnuc_list == []):
        return

    ni56_arr = np.array(get_ni56(results_dir))

    dtnuc_list = np.array([float(dtnuc[len('dt'):]) for dtnuc in dtnuc_list])

    # Sort the lists

    ni56_arr = np.array([x for _,x in sorted(zip(dtnuc_list,ni56_arr))])
    dtnuc_list = sorted(dtnuc_list)

    plt.plot(dtnuc_list, ni56_arr, linestyle=linestyles[0], lw = 4.0)

    plt.xscale('log')
    plt.xlabel(r"Nuclear burning timestep factor $\Delta t_{bX}$", fontsize=20)
    plt.ylabel(r"$^{56}$Ni generated (M$_{\odot}$)", fontsize=20)
    plt.tick_params(labelsize=20)
    plt.tight_layout()
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, get_diagfile(results_dir), 'diag')

    plt.close()



# Combined plot of both burning timestep factors

def burning_limiter(eps_filename, results_base):
    """Plot the effects of both burning timestep limiters (castro.dtnuc_e and castro.dtnuc_X)."""

    if (os.path.isfile(eps_filename)):
        return

    results_dir = results_base + 'burning_limiter_e/'

    # Get the list of parameter values we have tried

    dtnuc_list = wdmerger.get_parameter_list(results_dir)

    if (dtnuc_list == []):
        return

    ni56_arr = np.array(get_ni56(results_dir))

    dtnuc_list = np.array([float(dtnuc[len('dt'):]) for dtnuc in dtnuc_list])

    # Sort the lists

    ni56_arr = np.array([x for _,x in sorted(zip(dtnuc_list,ni56_arr))])
    dtnuc_list = sorted(dtnuc_list)

    plt.plot(dtnuc_list, ni56_arr, linestyle=linestyles[0], lw = 4.0, label=r'$\Delta t_{be}$')



    results_dir = results_base + 'burning_limiter_X/'

    # Get the list of parameter values we have tried

    dtnuc_list = wdmerger.get_parameter_list(results_dir)

    if (dtnuc_list == []):
        return

    ni56_arr = np.array(get_ni56(results_dir))

    dtnuc_list = np.array([float(dtnuc[len('dt'):]) for dtnuc in dtnuc_list])

    # Sort the lists

    ni56_arr = np.array([x for _,x in sorted(zip(dtnuc_list,ni56_arr))])
    dtnuc_list = sorted(dtnuc_list)

    plt.plot(dtnuc_list, ni56_arr, linestyle=linestyles[1], lw = 4.0, label=r'$\Delta t_{bX}$')

    plt.xscale('log')
    plt.xlabel(r"Nuclear burning timestep factor", fontsize=20)
    plt.ylabel(r"$^{56}$Ni generated (M$_{\odot}$)", fontsize=20)
    plt.tick_params(labelsize=20)
    plt.legend(loc='best', prop={'size':20})
    plt.tight_layout()
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, get_diagfile(results_dir), 'diag')

    plt.close()



# Effect of burning timestep limiter mode

def burning_limiter_mode(out_filename, results_base):
    """Create a table with the nickel generation as a function of the mode used for timestep limiting."""

    import os

    if os.path.isfile(out_filename):
        return

    # Get the list of parameter values we have tried

    results_dir = results_base + '/burning_limiter_mode'

    param_list = wdmerger.get_parameter_list(results_dir)

    if (param_list == []):
        return

    mode_list = [param[4:] for param in param_list]

    ni56_list = get_ni56(results_dir)

    comment = '%\n' + \
              '% Summary of the effect of the timestep limiting mode \n' + \
              '% on 56Ni yield in a white dwarf collision.\n' + \
              '% The first column is the limiter mode \n' + \
              '% and the second is the final nickel \n' + \
              '% mass produced, in solar masses.\n' + \
              '%\n'

    caption = 'The burning limiter modes are explained in \\autoref{sec:hydrocoupling}.'

    col1title = 'Mode'

    title = 'Burning Limiter Mode'

    label = 'burninglimitermode'

    write_ni56_table(results_dir, out_filename, mode_list, ni56_list, comment, col1title, label, title = title, caption = caption)



# Combined plot of effect of ODE tolerances

def ode_tolerances(eps_filename, results_base):
    """Plot the effects of all ODE tolerance options."""

    if (os.path.isfile(eps_filename)):
        return

    results_dir = results_base + 'spec_tol/'

    # Get the list of parameter values we have tried

    tol_list = wdmerger.get_parameter_list(results_dir)

    if (tol_list == []):
        return

    ni56_arr = np.array(get_ni56(results_dir))

    tol_list = [tol[len('tol'):] for tol in tol_list]
    tol_list = np.array([float(tol.replace('d','e')) for tol in tol_list])

    # Sort the lists

    ni56_arr = np.array([x for _,x in sorted(zip(tol_list,ni56_arr))])
    tol_list = sorted(tol_list)

    plt.plot(tol_list, ni56_arr, linestyle=linestyles[0], lw = 4.0, label=r'Species Tolerance')



    results_dir = results_base + 'temp_tol/'

    # Get the list of parameter values we have tried

    tol_list = wdmerger.get_parameter_list(results_dir)

    if (tol_list == []):
        return

    ni56_arr = np.array(get_ni56(results_dir))

    tol_list = [tol[len('tol'):] for tol in tol_list]
    tol_list = np.array([float(tol.replace('d','e')) for tol in tol_list])

    # Sort the lists

    ni56_arr = np.array([x for _,x in sorted(zip(tol_list,ni56_arr))])
    tol_list = sorted(tol_list)

    plt.plot(tol_list, ni56_arr, linestyle=linestyles[1], lw = 4.0, label=r'Temperature Tolerance')



    results_dir = results_base + 'enuc_tol/'

    # Get the list of parameter values we have tried

    tol_list = wdmerger.get_parameter_list(results_dir)

    if (tol_list == []):
        return

    ni56_arr = np.array(get_ni56(results_dir))

    tol_list = [tol[len('tol'):] for tol in tol_list]
    tol_list = np.array([float(tol.replace('d','e')) for tol in tol_list])

    # Sort the lists

    ni56_arr = np.array([x for _,x in sorted(zip(tol_list,ni56_arr))])
    tol_list = sorted(tol_list)

    plt.plot(tol_list, ni56_arr, linestyle=linestyles[2], lw = 4.0, label=r'Energy Tolerance')



    plt.xscale('log')
    plt.xlabel(r"ODE Error Tolerance", fontsize=20)
    plt.ylabel(r"$^{56}$Ni generated (M$_{\odot}$)", fontsize=20)
    plt.tick_params(labelsize=16)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, get_diagfile(results_dir), 'diag')

    plt.close()



# Effect of the temperature floor

def small_temp(out_filename, results_base):
    """Create a table with the nickel generation as a function of the temperature floor."""

    if (os.path.isfile(out_filename)):
        return

    # Get the list of parameter values we have tried

    results_dir = results_base + '/small_temp'

    param_list = wdmerger.get_parameter_list(results_dir)

    if (param_list == []):
        return

    temp_list = [float(param[1:]) for param in param_list]

    ni56_list = get_ni56(results_dir)

    # Sort the lists

    ni56_list = np.array([x for _,x in sorted(zip(temp_list,ni56_list))])
    temp_list = sorted(temp_list)

    # Return temp_list to scientific notation

    temp_list = np.array(["{:.1e}".format(t) for t in temp_list])

    comment = '%\n' + \
              '% Summary of the effect of the temperature floor \n' + \
              '% castro.small_temp on 56Ni yield in a white \n' + \
              '% dwarf collision.\n' + \
              '% The first column is the temperature floor and the \n' + \
              '% second is the final nickel mass produced, in solar masses.\n' + \
              '%\n'

    title = 'Effect of the temperature floor'

    col1title = '$T_{\\text{floor}}$'

    label = 'smalltemp'

    write_ni56_table(results_dir, out_filename, temp_list, ni56_list, comment, col1title, label, title = title)



# Effect of the minimum temperature for reactions

def t_min(out_filename, results_base):
    """Create a table with the nickel generation as a function of the minimum temperature for reactions."""

    if (os.path.isfile(out_filename)):
        return

    # Get the list of parameter values we have tried

    results_dir = results_base + '/T_min'

    param_list = wdmerger.get_parameter_list(results_dir)

    if (param_list == []):
        return

    temp_list = [float(param[1:]) for param in param_list]

    ni56_list = get_ni56(results_dir)

    # Sort the lists

    ni56_list = np.array([x for _,x in sorted(zip(temp_list,ni56_list))])
    temp_list = sorted(temp_list)

    # Return temp_list to scientific notation

    temp_list = np.array(["{:.1e}".format(t) for t in temp_list])

    comment = '%\n' + \
              '% Summary of the effect of the minimum temperature\n' + \
              '% for reactions on 56Ni yield in a white \n' + \
              '% dwarf collision.\n' + \
              '% The first column is the temperature minimum and the \n' + \
              '% second is the final nickel mass produced, in solar masses.\n' + \
              '%\n'

    title = 'Effect of the temperature floor for reactions'

    col1title = '$T_{\\text{min}}$'

    label = 'tmin'

    write_ni56_table(results_dir, out_filename, temp_list, ni56_list, comment, col1title, label, title = title)



# Effect of the minimum temperature for reactions

def rho_min(out_filename, results_base):
    """Create a table with the nickel generation as a function of the minimum density for reactions."""

    if (os.path.isfile(out_filename)):
        return

    # Get the list of parameter values we have tried

    results_dir = results_base + '/rho_min'

    param_list = wdmerger.get_parameter_list(results_dir)

    if (param_list == []):
        return

    temp_list = [float(param[3:]) for param in param_list]

    ni56_list = get_ni56(results_dir)

    # Sort the lists

    ni56_list = np.array([x for _,x in sorted(zip(temp_list,ni56_list))])
    temp_list = sorted(temp_list)

    # Return temp_list to scientific notation

    temp_list = np.array(["{:.1e}".format(t) for t in temp_list])


    comment = '%\n' + \
              '% Summary of the effect of the minimum density\n' + \
              '% for reactions on 56Ni yield in a white \n' + \
              '% dwarf collision.\n' + \
              '% The first column is the density minimum and the \n' + \
              '% second is the final nickel mass produced, in solar masses.\n' + \
              '%\n'

    title = 'Effect of the density floor for reactions'

    col1title = '$\\rho_{\\text{min}}$'

    label = 'tmin'

    write_ni56_table(results_dir, out_filename, temp_list, ni56_list, comment, col1title, label, title = title)



# Effect of the method for integrating the energy/temperature equations

def burning_mode(out_filename, results_base):
    """Create a table with the nickel generation as a function of the burning_mode parameter."""

    if (os.path.isfile(out_filename)):
        return

    # Get the list of parameter values we have tried

    results_dir = results_base + '/burning_mode'

    param_list = wdmerger.get_parameter_list(results_dir)

    if (param_list == []):
        return

    mode_list = [param for param in param_list]

    ni56_list = get_ni56(results_dir)

    # Sort the lists

    ni56_list = np.array([x for _,x in sorted(zip(mode_list,ni56_list))])
    mode_list = sorted(mode_list)

    # Replace the burning mode numbers with names for output purposes

    try:
        mode_list[mode_list.index('0')] = 'Hydrostatic'
    except:
        pass

    try:
        mode_list[mode_list.index('1')] = 'Self-heating'
    except:
        pass

    try:
        mode_list[mode_list.index('2')] = 'Hybrid'
    except:
        pass

    try:
        mode_list[mode_list.index('3')] = 'Suppressed'
    except:
        pass

    comment = '%\n' + \
              '% Summary of the effect of the burning mode \n' + \
              '% for reactions on 56Ni yield in a white \n' + \
              '% dwarf collision.\n' + \
              '% The first column is the density minimum and the \n' + \
              '% second is the final nickel mass produced, in solar masses.\n' + \
              '%\n'

    caption = 'The burning modes are explained in \\autoref{sec:burner}.'

    col1title = 'Burning Mode'

    title = 'Burning Mode'

    label = 'burningmode'

    precision = 4 # Because burning mode 3 creates so little nickel

    write_ni56_table(results_dir, out_filename, mode_list, ni56_list, comment, col1title, label, precision = precision, title = title, caption = caption)



def amr(eps_filename, results_base):
    """Plot the effect of the refinement based on the nuclear burning rate."""

    import os

    if os.path.isfile(eps_filename):
        return

    import math
    import numpy as np
    from matplotlib import pyplot as plt

    results_dir = results_base + 'amr'

    # Loop over base number of zones per dimension.

    ncell_list = sorted([int(n[1:]) for n in os.listdir(results_dir)])

    coarse_list = []

    i = 0

    for ncell in ncell_list:

        results_dir = results_base + 'amr/n' + str(ncell)

        # Get the list of parameter values we have tried

        r_list = wdmerger.get_parameter_list(results_dir)

        if (r_list == []):
            return

        r_list = [float(r[1:]) for r in r_list]

        ni56_arr = get_ni56(results_dir)

        # Sort the lists

        ni56_arr = np.array([x for _,x in sorted(zip(r_list,ni56_arr))])
        r_list = np.array(sorted(r_list))

        zonesPerDim = [ncell * int(r) for r in r_list]

        plt.plot(zonesPerDim, ni56_arr, markers[i], markersize=12, linestyle=linestyles[i], lw = 4.0, label = 'n = ' + str(ncell))

        i += 1

        coarse_list.append(ni56_arr[0])

    # Add a plot connecting the uniform grid cases

    plt.plot(ncell_list, coarse_list, lw = 4.0, color = 'black', marker = 'o', markersize = 12, label = 'Uniform grid')

    plt.tick_params(axis='both', which='major', pad=10)
    plt.xscale('log', basex=2)
    plt.ylim([0.0, 0.6])
    plt.xlabel(r"Effective zones/dimension on finest level", fontsize=24)
    plt.ylabel(r"$^{56}$Ni generated (M$_{\odot}$)", fontsize=24)
    plt.tick_params(labelsize=20)
    plt.tight_layout()
    plt.legend(loc='best', prop={'size':20})
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, get_diagfile(results_dir), 'diag')

    plt.close()



def rho_T_sliceplots(output_base, results_base, smallplt = True, prefix = "",
                     domain_frac = 1.0, x_ticks = [2.0e9, 4.0e9], y_ticks = [2.0e9, 4.0e9], scale_exp = 9):
    """Create a rho/T sliceplot for every plotfile in a given directory."""

    import os

    if not os.path.isdir(output_base):
        os.makedirs(output_base)

    from PIL import Image

    if (smallplt):
        plt_prefix = 'smallplt'

    param_list = wdmerger.get_parameter_list(results_base)

    if (param_list == []):
        return

    dir_list = [results_base + param + '/output/' for param in param_list]

    for param, directory in zip(param_list, dir_list):

        output_dir = output_base + param

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        plt_list = wdmerger.get_plotfiles(directory, plt_prefix)
        eps_list = [output_dir + "/rho_T_slice" + prefix + "_" + param + "_t_" + str("%.2f" % wdmerger.get_time_from_plotfile(directory + pltfile)) + '.eps' for pltfile in plt_list]

        # Generate the plot

        for eps_file, pltfile in zip(eps_list, plt_list):

            if os.path.isfile(eps_file):
                continue

            print("Generating plot with filename " + eps_file)

            wdmerger.rho_T_sliceplot(eps_file, directory + pltfile,
                                     domain_frac = domain_frac,
                                     x_ticks = x_ticks,
                                     y_ticks = y_ticks,
                                     scale_exp = scale_exp)

        jpg_list = [eps.replace('eps', 'jpg') for eps in eps_list]

        mpg_filename = output_dir + "/rho_T_slice" + prefix + "_" + param + ".mpg"

        if not os.path.isfile(mpg_filename):
            wdmerger.make_movie(output_dir, jpg_list, mpg_filename)



# Main execution: do all of the analysis routines.

if __name__ == "__main__":
    """Generate the plots and tables for the collision tests."""

    results_base = 'results/'
    plots_dir = 'plots/'

    burning_limiter_e(plots_dir + "dtnuc_e_max_Ni56.eps", results_base)
    burning_limiter_X(plots_dir + "dtnuc_X_max_Ni56.eps", results_base)
    burning_limiter(plots_dir + "dtnuc_max_Ni56.eps", results_base)
    small_temp(plots_dir + "small_temp.tbl", results_base)
    burning_mode(plots_dir + "burning_mode.tbl", results_base)
    t_min(plots_dir + "t_min.tbl", results_base)
    rho_min(plots_dir + "rho_min.tbl", results_base)
    ode_tolerances(plots_dir + "ode_tolerance.eps", results_base)
    amr(plots_dir + "amr.eps", results_base)
