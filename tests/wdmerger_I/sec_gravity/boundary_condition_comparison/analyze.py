import yt
import numpy as np
from matplotlib import pyplot as plt
import wdmerger

ell = range(21)

l2 = np.zeros(len(ell))

# Load exact results, for comparison

pf_exact = yt.load("results/true/plt00000")
print pf_exact.domain_dimensions
exact_data = pf_exact.covering_grid(level=0, left_edge=[0.0,0.0,0.0], dims=pf_exact.domain_dimensions)['phiGrav']
exact_L2 = np.sqrt((exact_data**2).sum())

index = 0

for l in ell:

    pf_name = "results/" + str(l) + "/plt00000"

    pf = yt.load(pf_name)

    plot_data = pf.covering_grid(level=0,left_edge=[0.0,0.0,0.0], dims=pf.domain_dimensions)['phiGrav']                                           

    L2_field = (plot_data - exact_data)**2

    l2[index] = np.sqrt(L2_field.sum())  / exact_L2

    index += 1

eps_filename = 'plots/bc_comparison.eps'

plt.plot(ell,l2)
plt.xlabel("Maximum multipole order")
plt.ylabel("Relative L2 error")
plt.yscale('log')
plt.savefig(eps_filename)
wdmerger.insert_commits_into_eps(eps_filename, pf_name, 'plot')
