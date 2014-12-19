import yt
import numpy as np
from matplotlib import pyplot as plt

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

            # Uniform density cube; solution is courtesy of Hummer (1996)

            c = yt.YTArray(np.zeros((2,3,dim[0],dim[1],dim[2])), 'cm')

            c[0,2] = -radius - zz
            c[1,2] =  radius - zz
            c[0,1] = -radius - yy
            c[1,1] =  radius - yy
            c[0,0] = -radius - xx
            c[1,0] =  radius - xx

            for ii in range(2):
                for jj in range(2):
                    for ll in range(3):

                        num1 = ( (c[ii,ll]**2 + c[jj,(ll+1) % 3]**2 + c[1,(ll+2) % 3]**2)**(0.5) + c[1,(ll+2) % 3] )**3
                        num2 = ( (c[ii,ll]**2 + c[jj,(ll+1) % 3]**2 + c[0,(ll+2) % 3]**2)**(0.5) - c[0,(ll+2) % 3] )
                        den1 = ( (c[ii,ll]**2 + c[jj,(ll+1) % 3]**2 + c[1,(ll+2) % 3]**2)**(0.5) - c[1,(ll+2) % 3] )
                        den2 = ( (c[ii,ll]**2 + c[jj,(ll+1) % 3]**2 + c[0,(ll+2) % 3]**2)**(0.5) + c[0,(ll+2) % 3] )**3

                        exact = exact + (0.5 * Gconst * density) * 0.5 * (-1)**(ii+jj) * ( c[ii,ll] * c[jj,(ll+1) % 3] * \
                                        np.log( num1 * num2 / (den1 * den2) ) )

            for ii in range(2):
                for jj in range(2):
                    for kk in range(2):
                        for ll in range(3):
                            exact = exact + (0.5 * Gconst * density) * (-1)**(ii+jj+kk+1) * c[ii,ll]**2 * \
                                  np.arctan2( c[ii,ll] * c[kk,(ll+2) % 3], \
                                              c[ii,ll]**2 + c[jj,(ll+1) % 3]**2 + c[jj,(ll+1) % 3] * \
                                              (c[ii,ll]**2 + c[jj,(ll+1) % 3]**2 + c[kk,(ll+2) % 3]**2)**(0.5) )

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
