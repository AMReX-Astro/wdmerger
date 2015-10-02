import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger
import wdmerger_grid_diag_analysis as grid_diag
import wdmerger_star_diag_analysis as star_diag
import os

markers = ['o', '+', '.', ',', '*']

# Long orbits

diag_list = []

for rot in ['0', '1']:
    for ratio in ['equal', 'unequal']:

        diag_filename = "results/" + ratio + "/rot" + rot + "/output/grid_diag.out"

        eps_filename = "plots/" + ratio + "_energy_rot" + rot + ".eps"

        if (not os.path.isfile(eps_filename)):

            grid_diag.plot_energy_error(diag_filename, eps_filename)

        eps_filename = "plots/" + ratio + "_angular_momentum_rot" + rot + ".eps"

        if (not os.path.isfile(eps_filename)):

            grid_diag.plot_angular_momentum_error(diag_filename, eps_filename)

        eps_filename = "plots/" + ratio + "_gw_rot" + rot + ".eps"

        if (not os.path.isfile(eps_filename)):

            grid_diag.plot_gw_signal(diag_filename, eps_filename, n_orbits=2, do_analytical=1)

        diag_filename = "results/" + ratio + "/rot" + rot + "/output/star_diag.out"

        eps_filename = "plots/" + ratio + "_location_rot" + rot + ".eps"

        if (not os.path.isfile(eps_filename)):
                
            star_diag.plot_wd_location(diag_filename, eps_filename)

        # Store this for the multipanel plot comparing all options

        diag_list.append(diag_filename)

eps_filename = "plots/circular_orbit_comparison.eps"

if (not os.path.isfile(eps_filename)):

    star_diag.plot_wd_location(diag_list, eps_filename)

# Source term tests

# Write out the change in energy to a LaTeX deluxetable

output_filename = 'plots/sources.table'

if (not os.path.isfile(output_filename)):

    gs_arr = ['1', '2', '3', '4']
    rs_arr = ['0', '1', '2', '3', '4']

    data_file = open(output_filename, 'w')

    data_file.write('\\begin{deluxetable*}{llllll}' + '\n')
    data_file.write('  \\tablecaption{ Change in energy after a single orbit, i.e. $|\Delta E / E|$. ' +
                                       '``rs\'\' is shorthand for the code parameter {\\tt castro.rot\\_source\\_type} and ' +
                                       '``gs\'\' is shorthand for the code parameter {\\tt castro.grav\\_source\\_type}. ' +
                                       'The parameter meanings are explained in the main text.}' + '\n')

    colHeadString = ' \\colhead{        } &' + '\n'

    j = 0
    for rs in rs_arr:
        colHeadString += '              \\colhead{ rs = ' + rs + ' }'
        if (j < len(rs_arr) - 1):
            colHeadString += ' & ' + '\n'
        else:
            colHeadString += ' } '
        j += 1

    data_file.write('  \\tablehead{'  + colHeadString + '\n')
    data_file.write('  \\startdata' + '\n')

    i = 0

    for gs in gs_arr:

        j = 0

        line_str = '    gs = ' + str(gs) + ' & '

        for rs in rs_arr:

            diag_filename = "results/sources/gs" + gs + "/rs" + rs + "/output/grid_diag.out"

            energy = star_diag.get_column("TOTAL ENERGY", diag_filename)

            dE = abs( (energy[-1] - energy[0]) / energy[0] )

            new_str = "%.1e" % dE
            [mantissa, exponent] = new_str.split('e')

            # Remove leading zeros from exponent, and get rid of plus sign.

            exponent = exponent[0].replace("+","") + exponent[1:].lstrip("0")

            # If our exponent really was zero, the strip will remove all of them, so replace it.

            if (exponent == ''):
                exponent = ' 0'

            # Now replace with the correct LaTeX syntax

            new_str = mantissa + " \\times 10^{" + exponent + "}"

            line_str += "$ " + new_str + " $"

            if (j == len(rs_arr) - 1):
                if (i != len(gs_arr) - 1):
                    line_str += " \\\\" 
                line_str += "\n"
            else:
                line_str += " & "

            j += 1

        data_file.write(line_str)

        i += 1

    data_file.write('  \\enddata' + '\n')
    data_file.write('  \\label{table:sources}' + '\n')

    [castro_hash, boxlib_hash, wdmerger_hash] = wdmerger.get_git_commits_from_diagfile(diag_filename)
    data_file.write('  % Castro git hash: ' + castro_hash + '\n')
    data_file.write('  % BoxLib git hash: ' + boxlib_hash + '\n')
    data_file.write('  % wdmerger git hash: ' + wdmerger_hash + '\n')
    data_file.write('\end{deluxetable*}' + '\n')
    data_file.close()


# Boundary conditions

# Determine the fractional change in distance of the WD binary after one orbit,
# and see whether it depends on boundary conditions.

eps_filename = 'plots/gravity_bcs.eps'

order_arr = [0, 2, 4, 6, 8, 10, 12, 14, 16]
dist_error_arr = np.zeros(len(order_arr))

if (not os.path.isfile(eps_filename)):

    j = 0

    for rot in ['0', '1']:

        i = 0

        for order in order_arr:

            diag_filename = "results/gravity_bcs/rot" + rot + "/order" + str(order) + "/output/star_diag.out"

            dist = star_diag.get_column('WD DISTANCE', diag_filename)
            
            dist_error_arr[i] = abs( (dist[-1] - dist[0]) / dist[0] )

            i += 1

        if (rot == '0'):
            label = "Inertial Frame"
        elif (rot == '1'):
            label = "Rotating Frame"

        plt.plot(order_arr, dist_error_arr, marker=markers[j], ms=12.0, lw=4.0, label=label)

        j += 1

    plt.xlabel('Maximum Multipole Moment', fontsize=20)
    plt.ylabel('Distance Change', fontsize=20)
    plt.yscale('log')
    plt.tick_params(labelsize=16)
    plt.legend()
    plt.savefig(eps_filename)
    wdmerger.insert_commits_into_eps(eps_filename, diag_filename, 'diag')

    plt.close()



# Spatial resolution

# Plot the WD distance over time, for the various resolutions.

for rot in ['0', '1']:

    eps_filename = 'plots/spatial_convergence_rot' + rot + '.eps'

    if (os.path.isfile(eps_filename)):
        continue

    diag_list = []
    label_list = []

    for ncell in ['256', '512', '1024']:
    
        diag_filename = 'results/spatial_convergence/rot' + rot + '/n' + ncell + '/output/star_diag.out'

        diag_list.append(diag_filename)

        label_list.append('n = ' + ncell)

    # Plot the first tenth of an orbit.

    n_orbits = 0.1

    star_diag.plot_wd_distance(diag_list, eps_filename, n_orbits, label_list)



# Time resolution

# Same as above, for various timesteps.

for rot in ['0', '1']:

    eps_filename = 'plots/time_convergence_rot' + rot + '.eps'

    if (os.path.isfile(eps_filename)):
        continue

    diag_list = []
    label_list = []

    for dt in ['0.01', '0.005', '0.0025', '0.00125']:
    
        diag_filename = 'results/time_convergence/rot' + rot + '/dt' + dt + '/output/star_diag.out'

        diag_list.append(diag_filename)

        label_list.append('dt = ' + dt)

    star_diag.plot_wd_distance(diag_list, eps_filename, labels=label_list)
