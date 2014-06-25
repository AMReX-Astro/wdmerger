import numpy as np
import pylab

data = np.loadtxt("scaling.txt")

cores = data[:,0]*data[:,1]
time = data[:,2]

pylab.loglog(cores, time, "o-", color="r", 
             label=r"3 levels, jumps of 2$\times$ (0 $\rightarrow$ 1) and 4$\times$ (1 $\rightarrow$ 2)")

ideal = time[0]*cores[0]/cores

pylab.loglog(cores, ideal, "--", color="k")

pylab.ylim(100,3000)

pylab.xlabel("number of cores")
pylab.ylabel("average time to advance the coarse level")
pylab.title("NERSC Edison Scaling for WD Merger", fontsize=16)

pylab.legend(frameon=False)

pylab.tight_layout()

pylab.savefig("edison-wdmerger-scaling.png")
pylab.savefig("edison-wdmerger-scaling.eps")
