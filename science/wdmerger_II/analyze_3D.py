import wdmerger

linestyles = ['-', '--', ':', '-', '--', ':']
markers = ['', '', '', 'o', 's', 'D']
colors = ['b', 'g', 'r', 'c', 'm', 'k']

def get_ni56(results_dir):
    """Return a list of maximum 56Ni production from all completed sub-directories in results_dir."""

    import numpy as np

    dir_list = sorted(wdmerger.get_parameter_list(results_dir))

    diag_filename_list = [results_dir + '/' + directory + '/output/species_diag.out' for directory in dir_list]

    ni56_arr = [np.amax(wdmerger.get_column('Mass ni56', diag_filename)) for diag_filename in diag_filename_list]

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



def write_ni56_table(results_dir, out_filename, v_list, n_list, comment, title, caption, col1title, label):
    """Write out a LaTeX table of 56Ni mass as a function of a given variable."""

    import os

    if os.path.isfile(out_filename):
        return

    wdmerger.insert_commits_into_txt(out_filename, get_diagfile(results_dir), 'diag')

    out_file = open(out_filename, 'a')

    out_file.write(comment)

    out_file.write('\\begin{deluxetable}{cc}' + '\n')
    out_file.write('  \\tablecaption{ ' + title + ' \label{table:' + label + '} }' + '\n')

    # Pad the table columns so that they fix the text column. This is an experimentally derived solution.

    padlabel = '\\tbl' + label + 'pad'
    padwidth = '0.325in'

    out_file.write('  \\newlength{' + padlabel + '}' + '\n')
    out_file.write('  \\setlength{' + padlabel + '}{' + padwidth + '}' + '\n')

    padstr = '\\hspace{\\tbl' + label + 'pad}'

    out_file.write('  \\tablehead{ \\colhead{ ' + padstr + col1title + padstr + ' } & \\colhead{ ' + padstr + 'Max. $^{56}$Ni (${M_{\odot}}$)' + padstr + ' } }' + '\n')

    out_file.write('  \\startdata' + '\n')

    for v, n in zip(v_list, n_list):

        line_str = '  ' + padstr + str(v) + padstr + ' & ' + padstr + "%.3f" % float(n) + padstr

        if (v != v_list[-1]):
            line_str += " \\\\"  + "\n" 
        else:
            line_str += "\n"

        out_file.write(line_str)

    out_file.write('  \\enddata' + '\n')
    out_file.write('  \\tablecomments{ ' + caption + ' }' + '\n')
    out_file.write('\end{deluxetable}' + '\n')
    out_file.close()



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
    """Generate the plots and tables for the 3D collisions."""

    mass_P = '0.64'
    mass_S = '0.64'

    mass_string = "_m_P_" + mass_P + "_m_S_" + mass_S

    c_frac = '0.5'
    o_frac = '0.5'

    results_base = 'results/3D/'
    plots_dir = 'plots/'

    impact_parameter(plots_dir + "impact_parameter.eps", results_base + "mass_P_" + mass_P + "/mass_S_" + mass_S + "/c" + c_frac + "/o" + o_frac + "/impact_parameter")

    mass_P = '0.80'
    mass_S = '0.60'

    b = '0.8'

    gravitational_wave_signal(plots_dir + "gw_signal_3D.eps", results_base + "mass_P_" + mass_P + "/mass_S_" + mass_S + "/c" + c_frac + "/o" + o_frac + "/gw_signal/b" + b)
