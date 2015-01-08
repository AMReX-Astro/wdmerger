import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger

ncell_arr = [16, 32, 64, 128, 256]

num_problems = 3

l2 = np.zeros((num_problems,len(ncell_arr)))

diameter = yt.YTQuantity(2.0e0, 'cm')
density  = yt.YTQuantity(1.0e3, 'g/cm**3')

radius   = diameter / 2.0e0

ambient_density = yt.YTQuantity(1.0e-8, 'g/cm**3')

Gconst   = yt.YTQuantity(6.67428e-8, 'cm**3*g**-1*s**-2')

for p in range(num_problems):

    problem = p + 1

    print "Now doing the analysis for problem", problem

    index = 0

    for ncell in ncell_arr:

        print "Now doing the comparison for ncell =", ncell

        # First we load in the numerical data

        pf_name = "results/problem" + str(problem) + "/" + str(ncell) + "/plt00000"

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

        index += 1

    plt.plot(ncell_arr,l2[p,:],label=title)

# Set up an analytical curve that represents perfect second-order convergence

maxerr = max(l2[:,0])
second_order = 3.0 * maxerr * (np.array(ncell_arr) / ncell_arr[0])**(-2.0)

plt.plot(ncell_arr, second_order, label='Second order convergence')

plt.xlim([ncell_arr[0], ncell_arr[len(ncell_arr)-1]])
plt.xlabel("Number of cells per dimension")
plt.ylabel(r"Relative L$^2$ error")
plt.xscale('log', basex=2)
plt.yscale('log')
plt.legend(loc='upper right')
plt.savefig('phi_comparison.png')
