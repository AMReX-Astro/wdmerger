import yt
import numpy as np
from matplotlib import pyplot as plt
import os

ncell = 64

# First we load in the numerical data. We want to use the last plotfile
# since that is at the time of interest.

dir = "results/" + str(ncell)

dir_contents = os.listdir(dir)

# Strip out non-plotfiles, then select the last one in the list.

plotfiles = sorted(filter(lambda s: s[0:3] == 'plt', dir_contents))

pf_name = "results/" + str(ncell) + "/" + plotfiles[len(plotfiles)-1]

pf = yt.load(pf_name)

plot = yt.SlicePlot(pf, 'z', "density", width=(1.0, 'cm'))

plot.save('density.eps')
