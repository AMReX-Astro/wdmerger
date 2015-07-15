import yt
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import wdmerger

# Determine our range of multiple orders by reading in the run_test.sh file

minell = None
maxell = None

file = open('run_test.sh', 'r')
for line in file:
    sline = line.split('=')
    if (len(sline) > 1):
        if (sline[0] == 'maxell'):
            maxell = int(sline[1])
        elif (sline[0] == 'minell'):
            minell = int(sline[1])

# Safety check

if (minell == None):
    print "Could not find minimum l value in run_test.sh; exiting."
    exit()

if (maxell == None):
    print "Could not find maximum l value in run_test.sh; exiting."
    exit()

ell = range(minell, maxell+1)

l2 = np.zeros(len(ell))

plt_name = "plt00000"

# Load exact results, for comparison

pf_exact = yt.load("results/true/output/" + plt_name)
exact_data = pf_exact.covering_grid(level=0, left_edge=pf_exact.domain_left_edge, dims=pf_exact.domain_dimensions)['phiGrav']
exact_L2 = np.sqrt((exact_data**2).sum())

index = 0

for l in ell:

    pf_name = "results/" + str(l) + "/output/" + plt_name

    pf = yt.load(pf_name)

    plot_data = pf.covering_grid(level=0,left_edge=pf_exact.domain_left_edge, dims=pf.domain_dimensions)['phiGrav']                                           

    L2_field = (plot_data - exact_data)**2

    l2[index] = np.sqrt(L2_field.sum())  / exact_L2

    index += 1

eps_filename = 'plots/bc_comparison.eps'

plt.plot(ell,l2, ls='-', lw=4.0, marker='o', ms=12.0)
plt.xlabel("Maximum multipole order", fontsize=20)
plt.ylabel(r"Relative $L^2$ error", fontsize=20)
plt.tick_params(labelsize=16)
plt.yscale('log')
plt.savefig(eps_filename)
wdmerger.insert_commits_into_eps(eps_filename, pf_name, 'plot')
