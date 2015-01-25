import yt
import numpy as np
from matplotlib import pyplot as plt
import os

ncell_arr = [  16,  128, 1024,  8192]
nproc_arr = [2048, 4096, 8192, 16384]

for ncell in ncell_arr:

    # Open up the standard output for analysis. It will be the numerically last file
    # starting with the designated output string.

    files = os.listdir("results/" + str(nprocs))

    files = sorted(filter(files[0:8] == "wdmerger",files))

    output_filename = "results/" + str(nprocs) + files[len(files)-1]

plt.savefig('weak_scaling.png')
