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



def burning_limiter_e(eps_filename, results_base, do_ni56 = True):
    """Plot the effect of the burning timestep factor castro.dtnuc_e."""

    if (os.path.isfile(eps_filename)):
        return
    else:
        print("Generating file %s" % eps_filename)

    results_dir = results_base + 'burning_limiter_e/'

    if not os.path.isdir(results_dir):
        return

    # Get the list of parameter values we have tried

    dtnuc_list = wdmerger.get_parameter_list(results_dir)

    if (dtnuc_list == []):
        return

    # We have two options: plot the total amount of nickel generated,
    # or plot the average atomic mass number at the end of the simulation.

    if (do_ni56):
        result_arr = np.array(get_ni56(results_dir))
    else:
        result_arr = np.array(get_abar(results_dir))

    # Don't make the plot if we have fewer than two data points.

    if len(result_arr) < 2:
        return

    dtnuc_list = np.array([float(dtnuc[len('dt'):]) for dtnuc in dtnuc_list])

    # Sort the lists

    result_arr = np.array([x for _,x in sorted(zip(dtnuc_list,result_arr))])
    dtnuc_list = sorted(dtnuc_list)

    plt.plot(dtnuc_list, result_arr, linestyle=linestyles[0], lw = 4.0, color = colors[0])

    plt.xscale('log')
    plt.xlabel(r"Nuclear burning timestep factor $f_{be}$", fontsize=20)
    if (do_ni56):
        plt.ylabel(r"$^{56}$Ni generated (M$_{\odot}$)", fontsize=20)
    else:
        plt.ylabel(r"$\overline{A}$", fontsize=20)
    plt.tick_params(labelsize=16)
    plt.tight_layout()
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, get_diagfile(results_dir), 'diag')
    png_filename = eps_filename[0:-3] + "png"
    plt.savefig(png_filename)

    plt.close()



def burning_limiter_X(eps_filename, results_base, do_ni56 = True):
    """Plot the effect of the burning timestep limiter castro.dtnuc_X."""

    if (os.path.isfile(eps_filename)):
        return
    else:
        print("Generating file %s" % eps_filename)

    results_dir = results_base + 'burning_limiter_X/'

    if not os.path.isdir(results_dir):
        return

    # Get the list of parameter values we have tried

    dtnuc_list = wdmerger.get_parameter_list(results_dir)

    if (dtnuc_list == []):
        return

    if (do_ni56):
        result_arr = np.array(get_ni56(results_dir))
    else:
        result_arr = np.array(get_abar(results_dir))

    # Don't make the plot if we have fewer than two data points.

    if len(result_arr) < 2:
        return

    dtnuc_list = np.array([float(dtnuc[len('dt'):]) for dtnuc in dtnuc_list])

    # Sort the lists

    result_arr = np.array([x for _,x in sorted(zip(dtnuc_list,result_arr))])
    dtnuc_list = sorted(dtnuc_list)

    plt.plot(dtnuc_list, result_arr, linestyle=linestyles[0], lw = 4.0, color = colors[0])

    plt.xscale('log')
    plt.xlabel(r"Nuclear burning timestep factor $f_{bX}$", fontsize=20)
    if (do_ni56):
        plt.ylabel(r"$^{56}$Ni generated (M$_{\odot}$)", fontsize=20)
    else:
        plt.ylabel(r"$\overline{A}$", fontsize=20)
    plt.tick_params(labelsize=16)
    plt.tight_layout()
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, get_diagfile(results_dir), 'diag')
    png_filename = eps_filename[0:-3] + "png"
    plt.savefig(png_filename)

    plt.close()



# Combined plot of both burning timestep factors

def burning_limiter(eps_filename, results_base, do_ni56 = True):
    """Plot the effects of both burning timestep limiters (castro.dtnuc_e and castro.dtnuc_X)."""

    if (os.path.isfile(eps_filename)):
        return
    else:
        print("Generating file %s" % eps_filename)

    results_dir = results_base + 'burning_limiter_e/'

    if not os.path.isdir(results_dir):
        return

    # Get the list of parameter values we have tried

    dtnuc_e_list = wdmerger.get_parameter_list(results_dir)

    if (dtnuc_e_list == []):
        return

    if (do_ni56):
        result_e_arr = np.array(get_ni56(results_dir))
    else:
        result_e_arr = np.array(get_abar(results_dir))

    # Don't make the plot if we have fewer than two data points.

    if len(result_e_arr) < 2:
        return

    dtnuc_e_list = np.array([float(dtnuc[len('dt'):]) for dtnuc in dtnuc_e_list])

    # Sort the lists

    result_e_arr = np.array([x for _,x in sorted(zip(dtnuc_e_list,result_e_arr))])
    dtnuc_e_list = sorted(dtnuc_e_list)



    results_dir = results_base + 'burning_limiter_X/'

    if not os.path.isdir(results_dir):
        return

    # Get the list of parameter values we have tried

    dtnuc_X_list = wdmerger.get_parameter_list(results_dir)

    if (dtnuc_X_list == []):
        return

    if (do_ni56):
        result_X_arr = np.array(get_ni56(results_dir))
    else:
        result_X_arr = np.array(get_abar(results_dir))

    # Don't make the plot if we have fewer than two data points.

    if len(result_X_arr) < 2:
        return

    dtnuc_X_list = np.array([float(dtnuc[len('dt'):]) for dtnuc in dtnuc_X_list])

    # Sort the lists

    result_X_arr = np.array([x for _,x in sorted(zip(dtnuc_X_list,result_X_arr))])
    dtnuc_X_list = sorted(dtnuc_X_list)



    results_dir = results_base + 'burning_limiter_eX/'

    if not os.path.isdir(results_dir):
        return

    # Get the list of parameter values we have tried

    dtnuc_eX_list = wdmerger.get_parameter_list(results_dir)

    if (dtnuc_eX_list == []):
        return

    if (do_ni56):
        result_eX_arr = np.array(get_ni56(results_dir))
    else:
        result_eX_arr = np.array(get_abar(results_dir))

    # Don't make the plot if we have fewer than two data points.

    if len(result_eX_arr) < 2:
        return

    dtnuc_eX_list = np.array([float(dtnuc[len('dt'):]) for dtnuc in dtnuc_eX_list])

    # Sort the lists

    result_eX_arr = np.array([x for _,x in sorted(zip(dtnuc_eX_list,result_eX_arr))])
    dtnuc_eX_list = sorted(dtnuc_X_list)



    plt.plot(dtnuc_e_list, result_e_arr, linestyle=linestyles[2], lw = 4.0, label='Energy limiting', color = colors[0])
    plt.plot(dtnuc_X_list, result_X_arr, linestyle=linestyles[0], lw = 4.0, label='Species limiting', color = colors[1])
    plt.plot(dtnuc_eX_list, result_eX_arr, linestyle=linestyles[1], lw = 4.0, label='Combined limiting', color = colors[2])

    plt.xscale('log')
    plt.xlabel(r"Nuclear burning timestep factor", fontsize=20)
    if (do_ni56):
        plt.ylabel(r"$^{56}$Ni generated (M$_{\odot}$)", fontsize=20)
    else:
        plt.ylabel(r"$\overline{A}$", fontsize=20)
    plt.tick_params(labelsize=16)
    plt.legend(loc='best', prop={'size':20})
    plt.tight_layout()
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, get_diagfile(results_dir), 'diag')
    png_filename = eps_filename[0:-3] + "png"
    plt.savefig(png_filename)

    plt.close()



# Effect of burning timestep limiter mode

def burning_limiter_mode(out_filename, results_base):
    """Create a table with the nickel generation as a function of the mode used for timestep limiting."""

    import os

    if os.path.isfile(out_filename):
        return
    else:
        print("Generating file %s" % out_filename)

    # Get the list of parameter values we have tried

    results_dir = results_base + '/burning_limiter_mode'

    if not os.path.isdir(results_dir):
        return

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
    else:
        print("Generating file %s" % eps_filename)

    results_dir = results_base + 'spec_tol/'

    if not os.path.isdir(results_dir):
        return

    # Get the list of parameter values we have tried

    s_tol_list = wdmerger.get_parameter_list(results_dir)

    if (s_tol_list == []):
        return

    ni56_s_arr = np.array(get_ni56(results_dir))

    tol_s_list = [tol[len('tol'):] for tol in tol_s_list]
    tol_s_list = np.array([float(tol.replace('d','e')) for tol in tol_s_list])

    # Sort the lists

    ni56_s_arr = np.array([x for _,x in sorted(zip(tol_s_list,ni56_s_arr))])
    tol_s_list = sorted(tol_s_list)



    results_dir = results_base + 'temp_tol/'

    if not os.path.isdir(results_dir):
        return

    # Get the list of parameter values we have tried

    tol_t_list = wdmerger.get_parameter_list(results_dir)

    if (tol_t_list == []):
        return

    ni56_t_arr = np.array(get_ni56(results_dir))

    tol_t_list = [tol[len('tol'):] for tol in tol_t_list]
    tol_t_list = np.array([float(tol.replace('d','e')) for tol in tol_t_list])

    # Sort the lists

    ni56_t_arr = np.array([x for _,x in sorted(zip(tol_t_list,ni56_t_arr))])
    tol_t_list = sorted(tol_t_list)



    results_dir = results_base + 'enuc_tol/'

    if not os.path.isdir(results_dir):
        return

    # Get the list of parameter values we have tried

    tol_e_list = wdmerger.get_parameter_list(results_dir)

    if (tol_e_list == []):
        return

    ni56_e_arr = np.array(get_ni56(results_dir))

    tol_e_list = [tol[len('tol'):] for tol in tol_e_list]
    tol_e_list = np.array([float(tol.replace('d','e')) for tol in tol_e_list])

    # Sort the lists

    ni56_e_arr = np.array([x for _,x in sorted(zip(tol_e_list,ni56_e_arr))])
    tol_e_list = sorted(tol_e_list)



    plt.plot(tol_s_list, ni56_s_arr, linestyle=linestyles[0], lw = 4.0, label=r'Species Tolerance')
    plt.plot(tol_t_list, ni56_t_arr, linestyle=linestyles[1], lw = 4.0, label=r'Temperature Tolerance')
    plt.plot(tol_e_list, ni56_e_arr, linestyle=linestyles[2], lw = 4.0, label=r'Energy Tolerance')
    plt.xscale('log')
    plt.xlabel(r"ODE Error Tolerance", fontsize=20)
    plt.ylabel(r"$^{56}$Ni generated (M$_{\odot}$)", fontsize=20)
    plt.tick_params(labelsize=16)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, get_diagfile(results_dir), 'diag')
    png_filename = eps_filename[0:-3] + "png"
    plt.savefig(png_filename)

    plt.close()



# Effect of the minimum temperature for reactions

def t_min(out_filename, results_base):
    """Create a table with the nickel generation as a function of the minimum temperature for reactions."""

    if (os.path.isfile(out_filename)):
        return
    else:
        print("Generating file %s" % out_filename)

    # Get the list of parameter values we have tried

    results_dir = results_base + '/T_min'

    if not os.path.isdir(results_dir):
        return

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
    else:
        print("Generating file %s" % out_filename)

    # Get the list of parameter values we have tried

    results_dir = results_base + '/rho_min'

    if not os.path.isdir(results_dir):
        return

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
    else:
        print("Generating file %s" % out_filename)

    # Get the list of parameter values we have tried

    results_dir = results_base + '/burning_mode'

    if not os.path.isdir(results_dir):
        return

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

    precision = 5 # Because the suppressed burning mode creates so little nickel

    write_ni56_table(results_dir, out_filename, mode_list, ni56_list, comment, col1title, label, precision = precision, title = title, caption = caption)



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



def rho_T_sliceplots(output_base, results_base, smallplt = True, prefix = "",
                     domain_frac = 1.0, x_ticks = [2.0e9, 4.0e9], y_ticks = [2.0e9, 4.0e9], scale_exp = 9):
    """Create a rho/T sliceplot for every plotfile in a given directory."""

    import os

    if not os.path.isdir(output_base):
        os.makedirs(output_base)

    from PIL import Image

    if (smallplt):
        plt_prefix = 'smallplt'

    r_list = wdmerger.get_parameter_list(results_base)

    # Strip out the non-refinement directories

    r_list = [r for r in r_list if r[0] == 'r']

    if (r_list == []):
        return

    dir_list = [results_base + r + '/output/' for r in r_list]

    for param, directory in zip(r_list, dir_list):

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

    mass_P = '0.64'
    mass_S = '0.64'

    results_base = 'results/collision_2D/mass_P_' + mass_P + '/mass_S_' + mass_S + '/'
    plots_dir = 'plots/'

    ncell_list = os.listdir(results_base)

    for ncell in ncell_list:

        r_list = os.listdir(results_base + ncell)

        for r in r_list:

            results_dir = results_base + ncell + '/' + r + '/self-heat/'

            # Plot the effect of the burning timestep limiters

            for do_ni56 in [True, False]:

                if (do_ni56):
                    label = 'max_Ni56'
                else:
                    label = 'abar'

                out_file_base = label + '_m_P_' + mass_P + '_m_S_' + mass_S + '_' + ncell + '_' + r + '.eps'

                burning_limiter_e(plots_dir + 'dtnuc_e_' + out_file_base, results_dir, do_ni56 = do_ni56)
                burning_limiter_X(plots_dir + 'dtnuc_X_' + out_file_base, results_dir, do_ni56 = do_ni56)
                burning_limiter(plots_dir + 'dtnuc_' + out_file_base, results_dir, do_ni56 = do_ni56)



            # Calculate the collision time

            results_dir = results_base + ncell + '/' + r + '/start/'

            collision_time(results_dir)



    #amr_nickel(plots_dir + 'amr_nickel.eps', results_base)
    #amr_detonation(plots_dir + 'amr_detonation.eps', results_base)

    #for n in ncell_list:
        #plot_dir = plots_dir + 'slices/dxnuc/' + n + '/'
        #prefix = '_' + n
        #rho_T_sliceplots(plot_dir, results_dir, prefix = prefix)
