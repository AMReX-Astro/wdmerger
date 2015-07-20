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
        data.append( map(float, line[0:-1].split()) )

    diag_file.close()

    # Convert the data list into a 2D numpy array.

    data = np.array(data)

    # Let's do some cleanup.

    col_names.pop(0)                                        # Get rid of the # at the beginning
    col_names = [string.strip() for string in col_names]    # Remove any leading or trailing whitespace
    col_names = filter(None, col_names)                     # Remove any remaining blank entries

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
