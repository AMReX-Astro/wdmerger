import yt
import numpy as np
from matplotlib import pyplot as plt

ncell_arr = [16, 32, 64, 128, 256]

l2 = np.zeros(len(ncell_arr))

index = 0

cube_a   = 2.0e0
cube_rho = 1.0e3

Gconst   = 6.67428e-8

c = np.zeros((2,3))

for ncell in ncell_arr:

    print "Now doing the comparison for ncell = ", ncell

    dx = 3.2e0 / ncell

    # First we'll create the array for the analytical solution

    exact = yt.YTArray(np.zeros((ncell,ncell,ncell)), 'erg/g')

    for k in range(ncell):
        zz = -1.6e0 + dx * (k + 0.5e0)
        for j in range(ncell):
            yy = -1.6e0 + dx * (j + 0.5e0)
            for i in range(ncell):
                xx = -1.6e0 + dx * (i + 0.5e0)

                c[0,2] = -cube_a/2 - zz
                c[1,2] =  cube_a/2 - zz
                c[0,1] = -cube_a/2 - yy
                c[1,1] =  cube_a/2 - yy
                c[0,0] = -cube_a/2 - xx
                c[1,0] =  cube_a/2 - xx

                phi = 0.0
                
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


                exact[i,j,k] = phi * 0.5 * Gconst * cube_rho

    exact_L2 = np.sqrt(dx**3 * (exact**2).sum())

    pf_name = "results/" + str(ncell) + "/plt00000"

    pf = yt.load(pf_name)

    plot_data = pf.covering_grid(level=0,left_edge=[-1.6,-1.6,-1.6], dims=pf.domain_dimensions)['phiGrav']

    L2_field = (plot_data - exact)**2

    l2[index] = np.sqrt(dx**3 * L2_field.sum())  / exact_L2

    index += 1

plt.plot(ncell_arr,l2)
plt.xlabel("Number of cells per dimension")
plt.ylabel("Relative L2 error")
plt.yscale('log')
plt.savefig('phi_comparison.png')
