# This Python program is used to create a plot displaying the sponge 
# function we use in the CASTRO hydrodynamics for the wdmerger problem.

import numpy as np
import matplotlib.pyplot as plt

def sponge(r):
    sp

rs = 0.75
rt = 0.85

r = np.linspace(0.0, 1.0, 1000)
f = np.zeros(len(r))

idx = np.where(r < rs)
f[idx] = 0.0

idx = np.where(r < rt)
idx = np.where(r[idx] >= rs)
f[idx] = 0.5 * (1.0 - np.cos(np.pi * (r[idx] - rs) / (rt - rs)))

idx = np.where(r >= rt)
f[idx] = 1.0

plt.plot(r, 1.0 - f, linewidth=4.0)
plt.xlabel('Radius', fontsize=20)
plt.ylabel(r'$1 - f_S$', fontsize=20)

plt.xlim([0.0, 1.0])
plt.ylim([-0.05, 1.05])

plt.tick_params(labelsize=16)

plt.tight_layout()

plt.savefig('sponge.eps')
