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

    diag_filename_list = [results_dir + '/' + directory + '/output/species_diag.out' for directory in dir_list]

    ni56_arr = [np.amax(spec_diag.get_column('Mass ni56', diag_filename)) - np.amin(spec_diag.get_column('Mass ni56', diag_filename))
                for diag_filename in diag_filename_list]

    return ni56_arr



def get_e_release(results_dir):
    """Return a list of energy generation from all completed sub-directories in results_dir."""

    import wdmerger_grid_diag_analysis as grid_diag
    import numpy as np

    dir_list = wdmerger.get_parameter_list(results_dir)

    diag_filename_list = [results_dir + '/' + directory + '/output/grid_diag.out' for directory in dir_list]

    e_arr = [np.amax(grid_diag.get_column('TOTAL ENERGY', diag_filename)) - np.amin(grid_diag.get_column('TOTAL ENERGY', diag_filename))
             for diag_filename in diag_filename_list]

    # Divide by 10**51

    e_arr = [e / 1.0e51 for e in e_arr]

    return e_arr



def get_diagfile(results_dir):
    """Return the name of a diagnostics file in results_dir that is completed.

       Useful for the functions that insert git commits into files and need a 
       representative diagnostic file to get the git commits from. We arbitrarily
       pick the last directory."""

    dir_list = wdmerger.get_parameter_list(results_dir)

    diag_filename_list = [results_dir + '/' + directory + '/output/species_diag.out' for directory in dir_list]

    diag_filename = diag_filename_list[-1]

    return diag_filename



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

    ni56_arr = np.array(get_ni56(results_dir))

    dtnuc_list = np.array([float(dtnuc[len('dt'):]) for dtnuc in dtnuc_list])

    # Sort the lists

    ni56_arr = np.array([x for _,x in sorted(zip(dtnuc_list,ni56_arr))])
    dtnuc_list = sorted(dtnuc_list)

    plt.plot(dtnuc_list, ni56_arr, linestyle=linestyles[0], lw = 4.0, label=r'$\Delta t_{be}$')



    results_dir = results_base + 'burning_limiter_X/'

    # Get the list of parameter values we have tried

    dtnuc_list = wdmerger.get_parameter_list(results_dir)

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

    mode_list = [param[4:] for param in param_list]

    ni56_list = get_ni56(results_dir)

    comment = '%\n' + \
              '% Summary of the effect of the timestep limiting mode \n' + \
              '% on 56Ni yield in a 2D white dwarf collision.\n' + \
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



# Effect of the dual energy parameter eta2

def eta2(out_filename, results_base):
    """Create a table with the nickel generation as a function of the dual energy eta2 parameter."""

    if (os.path.isfile(out_filename)):
        return

    # Get the list of parameter values we have tried

    results_dir = results_base + '/eta2'

    param_list = wdmerger.get_parameter_list(results_dir)

    eta_list = [float(param[3:]) for param in param_list]

    ni56_list = get_ni56(results_dir)

    # Sort the lists

    ni56_list = np.array([x for _,x in sorted(zip(eta_list,ni56_list))])
    eta_list = sorted(eta_list)

    comment = '%\n' + \
              '% Summary of the effect of the dual energy parameter \n' + \
              '% castro.dual_energy_eta2 on 56Ni yield in a 2D white \n' + \
              '% dwarf collision.\n' + \
              '% The first column is the value of eta2 and the second \n' + \
              '%is the final nickel mass produced, in solar masses.\n' + \
              '%\n'

    col1title = '$\\eta_2$'

    title = 'Effect of dual-energy parameter $\\eta_2$'

    label = 'eta2'

    write_ni56_table(results_dir, out_filename, eta_list, ni56_list, comment, col1title, label, title = title)



# Effect of the dual energy parameter eta3

def eta3(out_filename, results_base):
    """Create a table with the nickel generation as a function of the dual energy eta3 parameter."""

    if (os.path.isfile(out_filename)):
        return

    # Get the list of parameter values we have tried

    results_dir = results_base + '/eta3'

    param_list = wdmerger.get_parameter_list(results_dir)

    eta_list = [float(param[3:]) for param in param_list]

    ni56_list = get_ni56(results_dir)

    # Sort the lists

    ni56_list = np.array([x for _,x in sorted(zip(eta_list,ni56_list))])
    eta_list = sorted(eta_list)

    comment = '%\n' + \
              '% Summary of the effect of the dual energy parameter \n' + \
              '% castro.dual_energy_eta3 on 56Ni yield in a 2D white \n' + \
              '% dwarf collision.\n' + \
              '% The first column is the value of eta3 and the second \n' + \
              '% is the final nickel mass produced, in solar masses.\n' + \
              '%\n'

    title = 'Effect of dual-energy parameter $\\eta_3'

    col1title = '$\\eta_3$'

    label = 'eta3'

    write_ni56_table(results_dir, out_filename, eta_list, ni56_list, comment, col1title, label, title = title)



# Effect of the temperature floor

def small_temp(out_filename, results_base):
    """Create a table with the nickel generation as a function of the temperature floor."""

    if (os.path.isfile(out_filename)):
        return

    # Get the list of parameter values we have tried

    results_dir = results_base + '/small_temp'

    param_list = wdmerger.get_parameter_list(results_dir)

    temp_list = [float(param[1:]) for param in param_list]

    ni56_list = get_ni56(results_dir)

    # Sort the lists

    ni56_list = np.array([x for _,x in sorted(zip(temp_list,ni56_list))])
    temp_list = sorted(temp_list)

    # Return temp_list to scientific notation

    temp_list = np.array(["{:.1e}".format(t) for t in temp_list])

    comment = '%\n' + \
              '% Summary of the effect of the temperature floor \n' + \
              '% castro.small_temp on 56Ni yield in a 2D white \n' + \
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

    temp_list = [float(param[1:]) for param in param_list]

    ni56_list = get_ni56(results_dir)

    # Sort the lists

    ni56_list = np.array([x for _,x in sorted(zip(temp_list,ni56_list))])
    temp_list = sorted(temp_list)

    # Return temp_list to scientific notation

    temp_list = np.array(["{:.1e}".format(t) for t in temp_list])

    comment = '%\n' + \
              '% Summary of the effect of the minimum temperature\n' + \
              '% for reactions on 56Ni yield in a 2D white \n' + \
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

    temp_list = [float(param[3:]) for param in param_list]

    ni56_list = get_ni56(results_dir)

    # Sort the lists

    ni56_list = np.array([x for _,x in sorted(zip(temp_list,ni56_list))])
    temp_list = sorted(temp_list)

    # Return temp_list to scientific notation

    temp_list = np.array(["{:.1e}".format(t) for t in temp_list])


    comment = '%\n' + \
              '% Summary of the effect of the minimum density\n' + \
              '% for reactions on 56Ni yield in a 2D white \n' + \
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

    mode_list = [param for param in param_list]

    ni56_list = get_ni56(results_dir)

    # Sort the lists

    ni56_list = np.array([x for _,x in sorted(zip(mode_list,ni56_list))])
    mode_list = sorted(mode_list)

    # Replace the burning mode numbers with names for output purposes

    mode_list[mode_list.index('0')] = 'Hydrostatic'
    mode_list[mode_list.index('1')] = 'Self-heating'
    mode_list[mode_list.index('2')] = 'Hybrid'
    mode_list[mode_list.index('3')] = 'Suppressed'

    comment = '%\n' + \
              '% Summary of the effect of the burning mode \n' + \
              '% for reactions on 56Ni yield in a 2D white \n' + \
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



# Main execution: do all of the analysis routines.

if __name__ == "__main__":
    """Generate the plots and tables for the 2D and 3D collision tests."""

    mass_P = '0.64'
    mass_S = '0.64'

    mass_string = "_m_P_" + mass_P + "_m_S_" + mass_S

    results_base = 'results/2D/'
    plots_dir = 'plots/'

    burning_limiter_e(plots_dir + "dtnuc_e_max_Ni56" + mass_string + ".eps", results_base)
    burning_limiter_X(plots_dir + "dtnuc_X_max_Ni56" + mass_string + ".eps", results_base)
    burning_limiter(plots_dir + "dtnuc_max_Ni56" + mass_string + ".eps", results_base)
    eta2(plots_dir + "eta2" + mass_string + ".tbl", results_base)
    eta3(plots_dir + "eta3" + mass_string + ".tbl", results_base)
    small_temp(plots_dir + "small_temp" + mass_string + ".tbl", results_base)
    burning_mode(plots_dir + "burning_mode" + mass_string + ".tbl", results_base)
    t_min(plots_dir + "t_min" + mass_string + ".tbl", results_base)
    rho_min(plots_dir + "rho_min" + mass_string + ".tbl", results_base)
    ode_tolerances(plots_dir + "ode_tolerance.eps", results_base)
