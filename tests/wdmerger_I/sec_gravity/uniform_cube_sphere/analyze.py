import yt
import numpy as np
import math
from matplotlib import pyplot as plt

ncell_arr = [16, 32, 64, 128, 256]

num_problems = 3

l2 = np.zeros((num_problems,len(ncell_arr)))

diameter = 2.0e0
density  = 1.0e3

radius   = diameter / 2.0

ambient_density = 1.0e-8

Gconst   = 6.67428e-8

c = np.zeros((2,3))

for p in range(num_problems):

    problem = p + 1

    print "Now doing the analysis for problem", problem

    index = 0

    for ncell in ncell_arr:

        print "Now doing the comparison for ncell =", ncell

        dx = 3.2e0 / ncell

        # First we load in the numerical data

        pf_name = "results/problem" + str(problem) + "/" + str(ncell) + "/plt00000"

        pf = yt.load(pf_name)

        plot_data = pf.covering_grid(level=0,left_edge=[-1.6,-1.6,-1.6], dims=pf.domain_dimensions)['phiGrav']

        # Now we'll create the array for the analytical solution

        exact = yt.YTArray(np.zeros((ncell,ncell,ncell)), 'erg/g')

        if (problem == 1):

            mass = 4.0 / 3.0 * math.pi * radius**3 * density

        elif (problem == 2):

            densGrid = pf.covering_grid(level=0,left_edge=[-1.6,-1.6,-1.6], dims=pf.domain_dimensions)['density']
            mass = densGrid[np.where(densGrid > 2.0 * ambient_density)].v.sum() * dx**3

        for k in range(ncell):
            zz = -1.6e0 + dx * (k + 0.5e0)
            for j in range(ncell):
                yy = -1.6e0 + dx * (j + 0.5e0)
                for i in range(ncell):
                    xx = -1.6e0 + dx * (i + 0.5e0)

                    phi = 0.0

                    if (problem == 1 or problem == 2):

                        rr = (xx**2 + yy**2 + zz**2)**0.5

                        if (rr <= radius):
                            exact[i,j,k] = Gconst * mass * (3 * radius**2 - rr**2) / (2 * radius**3)
                        else:
                            exact[i,j,k] = Gconst * mass / rr

                    if (problem == 3):

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

                                    phi = phi + 0.5 * (-1)**(ii+jj) * ( c[ii,ll] * c[jj,(ll+1) % 3] * \
                                                np.log( num1 * num2 / (den1 * den2) ) )

                        for ii in range(2):
                            for jj in range(2):
                                for kk in range(2):
                                    for ll in range(3):
                                        phi = phi + (-1)**(ii+jj+kk+1) * c[ii,ll]**2 * \
                                              np.arctan2( c[ii,ll] * c[kk,(ll+2) % 3], \
                                                          c[ii,ll]**2 + c[jj,(ll+1) % 3]**2 + c[jj,(ll+1) % 3] * \
                                                          (c[ii,ll]**2 + c[jj,(ll+1) % 3]**2 + c[kk,(ll+2) % 3]**2)**(0.5) )

                        exact[i,j,k] = phi * 0.5 * Gconst * density

                    else:
                        
                        print "This is not a valid problem."

        # Now for problem 2, the only difference from problem 1 is that we normalize the mass
        # by the amount of mass actually on the grid for the sphere.

        if (problem == 2):
            densGrid = pf.covering_grid(level=0,left_edge=[-1.6,-1.6,-1.6], dims=pf.domain_dimensions)['density']
            actual_mass = densGrid[np.where(densGrid > 2.0 * ambient_density)].v.sum() * dx**3
            exact = exact * mass / actual_mass

        exact_L2 = np.sqrt(dx**3 * (exact**2).sum())

        L2_field = (plot_data - exact)**2

        l2[p,index] = np.sqrt(dx**3 * L2_field.sum())  / exact_L2

        index += 1

plt.plot(ncell_arr,l2[0,:],ncell_arr,l2[1,:],ncell_arr,l2[2,:])
plt.xlabel("Number of cells per dimension")
plt.ylabel("Relative L2 error")
plt.xscale('log')
plt.yscale('log')
plt.savefig('phi_comparison.png')
