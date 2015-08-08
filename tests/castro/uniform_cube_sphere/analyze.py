import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger

ncell_arr = [16, 32, 64, 128, 256]

num_problems = 3

l2 = np.zeros((num_problems,len(ncell_arr)))
convergence = np.zeros((num_problems,len(ncell_arr)))

# Get the problem setup data from the probin file

probin_file  = open('source/probin', 'r')
probin_lines = probin_file.readlines()

diameter = None
density = None
ambient_density = None

for line in probin_lines:
    lsplit = line.split()
    if (len(lsplit) >= 3):
        if (lsplit[0] == 'diameter'):
            diameter = lsplit[2]
            diameter = float(diameter.replace("d","e"))
            diameter = yt.YTQuantity(diameter, 'cm')
        elif (lsplit[0] == 'density'):
            density = lsplit[2]
            density = float(density.replace("d","e"))
            density = yt.YTQuantity(density, 'g/cm**3')
        elif (lsplit[0] == 'ambient_dens'):
            ambient_density = lsplit[2]
            ambient_density = float(ambient_density.replace("d","e"))
            ambient_density = yt.YTQuantity(ambient_density, 'g/cm**3')

if (diameter == None):
    print "Density not found in the probin file; exiting."
    exit()

if (density == None):
    print "Diameter not found in the probin file; exiting."
    exit()

if (ambient_density == None):
    print "Ambient density not found in the probin file; exiting."
    exit()

probin_file.close()

radius   = diameter / 2.0e0

# Get the gravitational constant from the CASTRO cgs constants file

Gconst = wdmerger.get_castro_const('Gconst')
Gconst = yt.YTQuantity(Gconst, 'cm**3*g**-1*s**-2')

for p in range(num_problems):

    problem = p + 1

    print "Now doing the analysis for problem", problem

    index = 0

    for ncell in ncell_arr:

        print "Now doing the comparison for ncell =", ncell

        # First we load in the numerical data

        pf_name = "results/p" + str(problem) + "/n" + str(ncell) + "/output/plt00000"

        pf = yt.load(pf_name)

        problo = pf.domain_left_edge
        probhi = pf.domain_right_edge
        dim    = pf.domain_dimensions

        plot_data = pf.covering_grid(level=0,left_edge=problo, dims=dim)['phiGrav']

        dx = (probhi - problo) / dim

        # Set up the grid representing the coordinate space

        x = problo[0] + dx[0] * (np.arange(dim[0]) + 0.5e0)
        y = problo[1] + dx[1] * (np.arange(dim[1]) + 0.5e0)
        z = problo[2] + dx[2] * (np.arange(dim[2]) + 0.5e0)
        xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")

        rr = (xx**2 + yy**2 + zz**2)**0.5

        # Now we'll create the array for the analytical solution

        exact = yt.YTArray(np.zeros(dim), 'erg/g')

        if (problem == 1):

            # Simple uniform sphere problem, with analytical mass estimate

            title = 'Sphere'

            mass = 4.0 / 3.0 * np.pi * radius**3 * density
            innerIdx = np.where(rr <= radius)
            outerIdx = np.where(rr > radius)
            exact[innerIdx] = Gconst * mass * (3 * radius**2 - rr[innerIdx]**2) / (2 * radius**3)
            exact[outerIdx] = Gconst * mass / rr[outerIdx]            

        elif (problem == 2):

            title = 'Normalized sphere'

            # Simple uniform sphere problem, with numerical mass estimate

            densGrid = pf.covering_grid(level=0,left_edge=problo, dims=dim)['density']
            mass = densGrid[np.where(densGrid > 2.0 * ambient_density)].v.sum() * dx[0] * dx[1] * dx[2]
            innerIdx = np.where(rr <= radius)
            outerIdx = np.where(rr > radius)
            exact[innerIdx] = Gconst * mass * (3 * radius**2 - rr[innerIdx]**2) / (2 * radius**3)
            exact[outerIdx] = Gconst * mass / rr[outerIdx]

        elif (problem == 3):

            title = 'Cube'

            # Uniform density cube; solution is courtesy of Waldvogel (1976)

            x = yt.YTArray(np.zeros((2,dim[0],dim[1],dim[2])), 'cm')
            y = yt.YTArray(np.zeros((2,dim[0],dim[1],dim[2])), 'cm')
            z = yt.YTArray(np.zeros((2,dim[0],dim[1],dim[2])), 'cm')

            x[0] = radius + xx
            x[1] = radius - xx
            y[0] = radius + yy
            y[1] = radius - yy
            z[0] = radius + zz
            z[1] = radius - zz

            for ii in range(2):
                for jj in range(2):
                    for kk in range(2):
                        r = (x[ii]**2 + y[jj]**2 + z[kk]**2)**0.5
                        exact += Gconst * density * (x[ii]*y[jj]*np.arctanh(z[kk] / r) + \
                                                     y[jj]*z[kk]*np.arctanh(x[ii] / r) + \
                                                     z[kk]*x[ii]*np.arctanh(y[jj] / r) - \
                                                     x[ii]**2 / 2.0 * np.arctan(y[jj]*z[kk]/(x[ii]*r)) - \
                                                     y[jj]**2 / 2.0 * np.arctan(z[kk]*x[ii]/(y[jj]*r)) - \
                                                     z[kk]**2 / 2.0 * np.arctan(x[ii]*y[jj]/(z[kk]*r)))

        else:
            
            print "This problem does not exist."

        exact_L2 = np.sqrt(dx[0]*dx[1]*dx[2] * (exact**2).sum())

        L2_field = (plot_data - exact)**2

        l2[p,index] = np.sqrt(dx[0]*dx[1]*dx[2] * L2_field.sum())  / exact_L2

        if (index > 0):
            convergence[p,index] = np.log(l2[p,index-1] / l2[p,index]) / np.log(ncell_arr[index] / ncell_arr[index-1])

        index += 1

    if (p == 0):
        lstyle = '--'
        mrkr = 'o'
    elif (p == 1):
        lstyle = '-.'
        mrkr = 'o'
    elif (p == 2):
        lstyle = ':'
        mrkr = 's'

    plt.plot(ncell_arr,l2[p,:], label=title, linestyle=lstyle, lw = 4.0, marker=mrkr, ms = 12.0)

# Write out the data to a text file.

file = open('plots/phi_comparison.dat','w')

file.write('     Problem 1          Problem 2          Problem 3   \n')
file.write('  Error    Order     Error    Order     Error    Order \n')

for n in range(len(ncell_arr)):
    outString = ''
    for p in range(num_problems):
        fmt1 = '{:4.3e}'
        fmt2 = '{: 4.3f}'
        outString += fmt1.format(l2[p,n]) + ' ' + fmt2.format(convergence[p,n]) + '   '

    file.write(outString)
    file.write('\n')

file.close()

# Set up an analytical curve that represents perfect second-order convergence.
# We'll place it a bit below the lowest curve on the grid, and give it a grayscale
# color by using the color = 'float' command where float is some number between 0 and 1.

minerr = min(l2[:,0])
second_order = (minerr / 3.0) * (np.array(ncell_arr) / ncell_arr[0])**(-2.0)

plt.plot(ncell_arr, second_order, label='Second order convergence', linestyle='-', lw = 3.0, color='0.75')

eps_filename = 'plots/phi_comparison.eps'

plt.xlim([ncell_arr[0], ncell_arr[len(ncell_arr)-1]])
plt.xlabel("Number of cells per dimension", fontsize=20)
plt.ylabel(r"Relative L$^2$ error", fontsize=20)
plt.xscale('log', basex=2)
plt.yscale('log')
plt.legend(loc='upper right', fontsize=12, handlelength=5)
plt.tick_params(labelsize=16)
plt.savefig(eps_filename)
wdmerger.insert_commits_into_eps(eps_filename, pf_name, 'plot')
