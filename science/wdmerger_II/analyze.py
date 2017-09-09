import numpy as np
import wdmerger
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import os

linestyles = ['-', '--', '-.', ':']
markers = ['o', 's', 'D', '^']



def get_ni56(results_dir):
    """Return a list of maximum 56Ni production from all completed sub-directories in results_dir."""

    import wdmerger_spec_diag_analysis as spec_diag

    dir_list = wdmerger.get_parameter_list(results_dir)

    diag_filename_list = [results_dir + '/' + directory + '/output/species_diag.out' for directory in dir_list]

    ni56_arr = [np.amax(spec_diag.get_column('Mass ni56', diag_filename)) for diag_filename in diag_filename_list]

    return ni56_arr



def get_diagfile(results_dir):
    """Return the name of a diagnostics file in results_dir that is completed.

       Useful for the functions that insert git commits into files and need a 
       representative diagnostic file to get the git commits from. We arbitrarily
       pick the last directory."""

    dir_list = wdmerger.get_parameter_list(results_dir)

    diag_filename_list = [results_dir + '/' + directory + '/output/species_diag.out' for directory in dir_list]

    diag_filename = diag_filename_list[-1]

    return diag_filename



def write_ni56_table(results_dir, out_filename, v_list, n_list, comment, caption, col1title, label):
    """Write out a LaTeX table of 56Ni mass as a function of a given variable."""

    wdmerger.insert_commits_into_txt(out_filename, get_diagfile(results_dir), 'diag')

    out_file = open(out_filename, 'a')

    out_file.write(comment)

    out_file.write('\\begin{deluxetable*}{ll}' + '\n')
    out_file.write('  \\tablecaption{ ' + caption + ' }' + '\n')

    out_file.write('  \\tablehead{ \\colhead{ ' + col1title + ' } & \\colhead{ Max. $^{56}$Ni }' + '\n')

    out_file.write('  \\startdata' + '\n')

    for v, n in zip(v_list, n_list):

        line_str = '  ' + str(v) + ' & ' + "%.3f" % float(n)

        if (v != v_list[-1]):
            line_str += " \\\\"  + "\n" 
        else:
            line_str += "\n"

        out_file.write(line_str)

    out_file.write('  \\enddata' + '\n')
    out_file.write('  \\label{table:' + label + '}' + '\n')
    out_file.write('\end{deluxetable*}' + '\n')
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
    plt.tick_params(labelsize=16)
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
    plt.tick_params(labelsize=16)
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

    print ni56_arr
    print dtnuc_list

    # Sort the lists

    ni56_arr = np.array([x for _,x in sorted(zip(dtnuc_list,ni56_arr))])
    dtnuc_list = sorted(dtnuc_list)

    print ni56_arr
    print dtnuc_list

    plt.plot(dtnuc_list, ni56_arr, linestyle=linestyles[1], lw = 4.0, label=r'$\Delta t_{bX}$')

    plt.xscale('log')
    plt.xlabel(r"Nuclear burning timestep factor", fontsize=20)
    plt.ylabel(r"$^{56}$Ni generated (M$_{\odot}$)", fontsize=20)
    plt.tick_params(labelsize=16)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, get_diagfile(results_dir), 'diag')

    plt.close()



# Effect of the C/O mass fraction

def c_fraction(out_filename, results_base):
    """Create a table with the nickel generation as a function of the C/O fraction."""

    if (os.path.isfile(out_filename)):
        return

    # Get the list of parameter values we have tried

    results_dir = results_base + '/co'

    param_list = wdmerger.get_parameter_list(results_dir)

    c_list = np.array([float(param.split('c')[1].split('o')[0]) for param in param_list])

    ni56_list = np.array(get_ni56(results_dir))

    # Sort the lists

    ni56_list = np.array([x for _,x in sorted(zip(c_list,ni56_list))])
    c_list = sorted(c_list)

    comment = '%\n' + \
              '% Summary of the effect of C/O mass fraction \n' + \
              '% on 56Ni yield in a 2D white dwarf collision.\n' + \
              '% The first column is the initial carbon mass \n' + \
              '% fraction and the second is the final nickel \n' + \
              '% mass produced, in solar masses.\n' + \
              '%\n'

    caption = ' Maximum $^{56}Ni$ mass ($M_{\odot}) produced as a function of the ' + \
              'initial mass fraction $C$. Since the white dwarfs are purely ' + \
              'carbon/oxygen the initial oxygen fraction is $1 - C$.'

    col1title = '$C$'

    label = 'co'

    write_ni56_table(results_dir, out_filename, c_list, ni56_list, comment, caption, col1title, label)




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

    caption = 'Maximum $^{56}Ni$ mass ($M_{\odot}) produced as a function of the ' + \
              'dual energy formalism parameter $\\eta_2$.'

    col1title = '$\\eta_2$'

    label = 'eta2'

    write_ni56_table(results_dir, out_filename, eta_list, ni56_list, comment, caption, col1title, label)



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

    caption = 'Maximum $^{56}Ni$ mass ($M_{\odot}) produced as a function of the ' + \
              'dual energy formalism parameter $\\eta_3$.'

    col1title = '$\\eta_3$'

    label = 'eta3'

    write_ni56_table(results_dir, out_filename, eta_list, ni56_list, comment, caption, col1title, label)



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

    caption = 'Maximum $^{56}Ni$ mass ($M_{\odot}) produced as a function of the ' + \
              'temperature floor.'

    col1title = '$T_{\\text{floor}}$'

    label = 'smalltemp'

    write_ni56_table(results_dir, out_filename, temp_list, ni56_list, comment, caption, col1title, label)



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

    caption = 'Maximum $^{56}Ni$ mass ($M_{\odot}) produced as a function of the ' + \
              'minimum temperature for allowing nuclear reactions.'

    col1title = '$T_{\\text{min}}$'

    label = 'tmin'

    write_ni56_table(results_dir, out_filename, temp_list, ni56_list, comment, caption, col1title, label)



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

    caption = 'Maximum $^{56}Ni$ mass ($M_{\odot}) produced as a function of the ' + \
              'minimum density for allowing nuclear reactions.'

    col1title = '$\\rho_{\\text{min}}$'

    label = 'tmin'

    write_ni56_table(results_dir, out_filename, temp_list, ni56_list, comment, caption, col1title, label)



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

    comment = '%\n' + \
              '% Summary of the effect of the burning mode \n' + \
              '% for reactions on 56Ni yield in a 2D white \n' + \
              '% dwarf collision.\n' + \
              '% The first column is the density minimum and the \n' + \
              '% second is the final nickel mass produced, in solar masses.\n' + \
              '%\n'

    caption = 'Maximum $^{56}Ni$ mass ($M_{\odot}) produced as a function of the ' + \
              'burning mode for allowing nuclear reactions. Burning mode 0 is a ' + \
              'hydrostatic burn, mode 1 is a self-heating burn, mode 2 is a hybrid ' + \
              'burn, and mode 3 is a suppressed burn. The meaning of these burning ' + \
              'modes is explained in the text.'

    col1title = 'Burning Mode'

    label = 'burningmode'

    write_ni56_table(results_dir, out_filename, mode_list, ni56_list, comment, caption, col1title, label)



# Effect of the imapct parameter

def impact_parameter(eps_filename, results_dir):
    """Plot the effect of the impact parameter."""

    import os

    if os.path.isfile(eps_filename):
        return

    import numpy as np
    from matplotlib import pyplot as plt

    # Get the list of parameter values we have tried

    dtnuc_list = wdmerger.get_parameter_list(results_dir)

    ni56_arr = get_ni56(results_dir)

    b_list = [dtnuc[len('b'):] for dtnuc in dtnuc_list]

    b_list = sorted(b_list)

    plt.plot(np.array(b_list), np.array(ni56_arr), markers[0], markersize = 12.0)

    xaxis_buffer = 0.025

    plt.xlim([float(b_list[0]) - xaxis_buffer, float(b_list[-1]) + xaxis_buffer])
    plt.xlabel(r"Impact parameter $b$", fontsize=24)
    plt.ylabel(r"$^{56}$Ni generated (M$_{\odot}$)", fontsize=24)
    plt.tick_params(labelsize=20)
    plt.tight_layout()
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, get_diagfile(results_dir), 'diag')

    plt.close()



# Gravitational wave signal

def gravitational_wave_signal(eps_filename, results_dir):
    """Plot the gravitational radiation waveform."""

    import os

    if os.path.isfile(eps_filename):
        return

    if not os.path.exists(results_dir):
        return

    from matplotlib import pyplot as plt

    print "Generating plot with filename " + eps_filename

    diag_file = results_dir + '/output/grid_diag.out'

    time     = wdmerger.get_column('TIME', diag_file)

    strain_p_1 = wdmerger.get_column('h_+ (axis 1)', diag_file)
    strain_x_1 = wdmerger.get_column('h_x (axis 1)', diag_file)

    strain_p_2 = wdmerger.get_column('h_+ (axis 2)', diag_file)
    strain_x_2 = wdmerger.get_column('h_x (axis 2)', diag_file)

    strain_p_3 = wdmerger.get_column('h_+ (axis 3)', diag_file)
    strain_x_3 = wdmerger.get_column('h_x (axis 3)', diag_file)

    plt.plot(time, strain_p_1 / 1.e-22, lw = 4.0, color = colors[0], linestyle = linestyles[0], marker = markers[0], markersize = 12, markevery = 250, label = r'$h^x_+$')
    plt.plot(time, strain_x_1 / 1.e-22, lw = 4.0, color = colors[1], linestyle = linestyles[1], marker = markers[1], markersize = 12, markevery = 250, label = r'$h^x_\times$')

    plt.plot(time, strain_p_2 / 1.e-22, lw = 4.0, color = colors[2], linestyle = linestyles[2], marker = markers[2], markersize = 12, markevery = 250, label = r'$h^y_+$')
    plt.plot(time, strain_x_2 / 1.e-22, lw = 4.0, color = colors[3], linestyle = linestyles[3], marker = markers[3], markersize = 12, markevery = 250, label = r'$h^y_\times$')

    plt.plot(time, strain_p_3 / 1.e-22, lw = 4.0, color = colors[4], linestyle = linestyles[4], marker = markers[4], markersize = 12, markevery = 250, label = r'$h^z_+$')
    plt.plot(time, strain_x_3 / 1.e-22, lw = 4.0, color = colors[5], linestyle = linestyles[5], marker = markers[5], markersize = 12, markevery = 250, label = r'$h^z_\times$')

    plt.tick_params(labelsize=20)

    plt.xlabel(r'$t\ \mathrm{(s)}$', fontsize=24)
    plt.ylabel(r'$h\, /\, 10^{-22}$', fontsize=24)
    plt.legend(loc = 'best', prop = {'size':20})
    plt.tight_layout()

    plt.savefig(eps_filename)

    wdmerger.insert_commits_into_eps(eps_filename, diag_file, 'diag')

    plt.close()



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
    c_fraction(plots_dir + "co_frac" + mass_string + ".tbl", results_base)
    eta2(plots_dir + "eta2" + mass_string + ".tbl", results_base)
    eta3(plots_dir + "eta3" + mass_string + ".tbl", results_base)
    small_temp(plots_dir + "small_temp" + mass_string + ".tbl", results_base)
    burning_mode(plots_dir + "burning_mode" + mass_string + ".tbl", results_base)
    t_min(plots_dir + "t_min" + mass_string + ".tbl", results_base)
    rho_min(plots_dir + "rho_min" + mass_string + ".tbl", results_base)
