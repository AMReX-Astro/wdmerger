import os
import wdmerger
import numpy as np
import matplotlib.pyplot as plt



#
# Get the data column from a wdmerger diagnostic output file whose header name is col_name.
#

def get_column(col_name, diag_filename):

    # Open up the file for reading. Get the names of the columns, as well as a 2D list
    # with all the data.

    diag_file = open(diag_filename,'r')

    vc_line = 'git'

    # Skip the first few lines, they store the version control information

    line = diag_file.readline()

    while (line.split()[2] == vc_line):
        line = diag_file.readline()

    # The very next line will be the column headers

    col_names = line.split('  ')

    # Now read in the data

    diag_list = diag_file.readlines()
    data = []
    for line in diag_list:
        data.append( list(map(float, line[0:-1].split())) )

    diag_file.close()

    # Convert the data list into a 2D numpy array.

    data = np.array(data)

    # Let's do some cleanup.

    col_names.pop(0)                                        # Get rid of the # at the beginning
    col_names = [string.strip() for string in col_names]    # Remove any leading or trailing whitespace
    col_names = list(filter(None, col_names))               # Remove any remaining blank entries

    # Obtain the column index and then return the column with that index.

    col_index = col_names.index(col_name)

    return data[:,col_index]



#
# Plot the absolute magnitude of the fractional change in energy over time.
#

def plot_energy_error(diag_filename, output_filename):

    diag_file = open(diag_filename, 'r')
    lines = diag_file.readlines()

    time = get_column("TIME", diag_filename)
    energy = get_column("TOTAL ENERGY",diag_filename)

    rot_period = 0.0

    dir = os.path.dirname(diag_filename)

    inputs_filename = dir + '/' + wdmerger.get_inputs_filename(dir)
    rot_period = wdmerger.get_inputs_var(inputs_filename, "castro.rotational_period")

    if (rot_period > 0.0):
        time = time / rot_period
        xlabel = "Time / Rotational Period"
    else:
        xlabel = "Time (s)"
    ylabel = "Relative Change in Energy"

    plt.plot(time, abs((energy - energy[0]) / energy[0]), lw = 4.0)
    plt.yscale('log')
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    plt.ylim([1.0e-8, 1.0e0])
    plt.tick_params(labelsize=16)
    plt.tight_layout()

    plt.savefig(output_filename)
    wdmerger.insert_commits_into_eps(output_filename, diag_filename, 'diag')

    plt.close()



#
# Plot the absolute magnitude of the fractional change in angular momentum over time.
#

def plot_angular_momentum_error(diag_filename, output_filename):

    diag_file = open(diag_filename, 'r')
    lines = diag_file.readlines()

    time = get_column("TIME", diag_filename)
    energy = get_column("ANG. MOM. Z",diag_filename)

    rot_period = 0.0

    dir = os.path.dirname(diag_filename)

    inputs_filename = dir + '/' + wdmerger.get_inputs_filename(dir)
    rot_period = wdmerger.get_inputs_var(inputs_filename, "castro.rotational_period")

    if (rot_period > 0.0):
        time = time / rot_period
        xlabel = "Time / Rotational Period"
    else:
        xlabel = "Time (s)"
    ylabel = "Relative Change in Angular Momentum"

    plt.plot(time, abs((energy - energy[0]) / energy[0]), lw = 4.0)
    plt.yscale('log')
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    plt.ylim([1.0e-8, 1.0e0])
    plt.tick_params(labelsize=16)
    plt.tight_layout()

    plt.savefig(output_filename)
    wdmerger.insert_commits_into_eps(output_filename, diag_filename, 'diag')

    plt.close()



#
# Plot the gravitational wave strain over time.
# Optionally we can plot the expected signal for a circular orbit.
#

def plot_gw_signal(diag_filename, output_filename, n_orbits = -1, do_analytical = 0):

    diag_file = open(diag_filename, 'r')

    time   = get_column("TIME", diag_filename)
    hplus  = get_column("h_+ (axis 3)", diag_filename)
    hcross = get_column("h_x (axis 3)", diag_filename)

    # Normalize time by rotational period.

    rot_period = 0.0

    dir = os.path.dirname(diag_filename)

    inputs_filename = dir + '/' + wdmerger.get_inputs_filename(dir)
    probin_filename = dir + '/probin'

    rot_period = wdmerger.get_inputs_var(inputs_filename, "castro.rotational_period")

    if (rot_period > 0.0):
        time = time / rot_period
        xlabel = "Time / Rotational Period"
        if (n_orbits > 0):
            idx = np.where(time < n_orbits)
            time   = time[idx]
            hplus  = hplus[idx]
            hcross = hcross[idx]
    else:
        xlabel = "Time (s)"

    ylabel = "Gravitational Wave Strain"

    markers = ['o', '+', '.', ',', '*']


    # Now work out the analytical result for a circular binary orbit, if desired.

    if (do_analytical == 1):

        # Get relevant CGS constants.

        G_const = wdmerger.get_castro_const('Gconst')
        M_solar = wdmerger.get_castro_const('M_solar')
        c_light = wdmerger.get_castro_const('c_light')
        parsec  = wdmerger.get_castro_const('parsec')

        # Relevant variables from the probin file.

        mass_P  = wdmerger.get_probin_var(probin_filename, 'mass_P') * M_solar
        mass_S  = wdmerger.get_probin_var(probin_filename, 'mass_S') * M_solar

        gw_dist = wdmerger.get_probin_var(probin_filename, 'gw_dist') # In kpc

        mu = mass_P * mass_S / (mass_P + mass_S)
        M_tot = mass_P + mass_S
        dist = gw_dist * 1.e3 * parsec
        omega = 2.0 * np.pi / rot_period
        coeff = -4 * G_const * mu / (c_light**4 * dist) * (G_const * M_tot * omega)**(2.0/3.0)
        hplus_true  = coeff * np.cos(2.0 * omega * rot_period * time)
        hcross_true = coeff * np.sin(2.0 * omega * rot_period * time)

        plt.plot(time,hplus_true,  lw = 4.0, ls='-')
        plt.plot(time,hcross_true, lw = 4.0, ls='--')

        plt.plot(time, hplus,  marker=markers[0], ms=12.0, markevery=200, label=r"$\plus$" + " polarization")
        plt.plot(time, hcross, marker=markers[4], ms=12.0, markevery=200, label=r"$\times$" + " polarization")

    else:

        plt.plot(time, hplus,  lw = 4.0, ls = '-',  label=r"$\plus$" + " polarization")
        plt.plot(time, hcross, lw = 4.0, ls = '--', label=r"$\times$" + " polarization")


    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)

    # Use the 'best' location for the legend, since for a generic function like this
    # it is hard to know ahead of time where the legend ought to go.
    # The alpha value controls the transparency, since we may end up covering some data.

    plt.legend(loc='best', fancybox=True)

    # The padding ensures that the lower-left ticks on the x- and y-axes don't overlap.

    plt.tick_params(labelsize=16, pad=10)

    # We have now increased the size of both the ticks and the axis labels,
    # which may have caused the latter to fall off the plot. Use tight_layout
    # to automatically adjust the plot to fix this.

    plt.tight_layout()

    # Save it into our designated file, which is usually EPS format.

    plt.savefig(output_filename)

    # Insert git commit hashes into this file from the various code sources.

    wdmerger.insert_commits_into_eps(output_filename, diag_filename, 'diag')

    plt.close()
