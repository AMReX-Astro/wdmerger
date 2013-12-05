# This script creates a 2-D plot of field value versus distance over a number
# of specifiable 1-D rays for your plotfiles.
from yt.mods import *
import matplotlib
import pylab as pl
import sys
import os
import fnmatch


# Specify field as well as  derived field if you are using one. You can choose
# to use your own plotfile list or have the script find them all based on
# their prefix. Also input the range you want to show on the graph, and a
# buffer if you want one. Setting find_fieldrange = True will make the script
# loop over your plotfiles and set the absolute min and max equal to the
# range. It can automatically find the max/min of the log of a field as well,
# but for anything more complex you will have to manually input ranges.
pflist = ['']
find_pfs = True
pfprefix = 'plt'
fld = 'Temp'
derived_field = 'logTemp'
use_derived_field = True
derived_field_log = True
fieldrange = [0,0]
find_fieldrange = True
fieldrange_buffer = [0.5,0.5]

# If you are using a derived field, specify it here. For example, you could
# return np.log10(data['Temp']) if you wanted to plot the log of temperature.
def _new_dfield(field, data):
    return np.log10(data['Temp'])

# Rays specified by starting point and ending point in the domain. Currently
# only up to 3 allowed. Each ray is a list with 2 elements, each having 3
# additional elements.
# Note: rays will be left-justified on the plot regardless of length.
rays = np.array([[[0.0, 5.12e9, 5.12e9],[1.024e10, 5.12e9, 5.12e9]]])

# Colors are b g r c m y k w. See http://matplotlib.org/api/colors_api.html.
# Size_in_inches refers to the size of your image, normal linewidth = 1.
colors = ['-b', '-g', '-r']
x_limits = [0.0, 1.024e10]
title = 'logTemp Along the X-Axis'
xlabel = 'Distance Along Ray (cm)'
ylabel = 'logTemp [log(K)]'
background_grid = True
linewidth = 0.7
size_in_inches = [10.0, 10.0]

# -----Everything below is automatic-----

# Find lengths of rays
raylen = np.zeros(len(rays))
for i, ray in enumerate(rays):
    raylen[i] = np.sqrt(sum(np.square(rays[i][1]-rays[i][0])))

# Find plotfiles
if find_pfs:
    pflist = []
    for p in os.listdir("."):
        if fnmatch.fnmatch(p, "%s*" % (pfprefix)) and not fnmatch.fnmatch(p,'*%s' % ('.yt')):
            pflist.append(p)
    pflist.sort()


# Loop over all plotfiles to find the absolute maximum and minimum of the data
# and set it equal to our range
if find_fieldrange:

    for i, plotf in enumerate(pflist):
        pf = load(plotf)
        if i==0:
            if use_derived_field and derived_field_log:
                fieldrange[0], fieldrange[1] = np.log10(pf.h.all_data().quantities[
                        "Extrema"](fld)[0])
            else: 
                fieldrange[0], fieldrange[1] = pf.h.all_data().quantities[
                    "Extrema"](fld)[0]
        print pf
        if use_derived_field and derived_field_log:
            mi, ma = np.log10(pf.h.all_data().quantities[
                    "Extrema"](fld)[0])
        else:
            mi, ma = pf.h.all_data().quantities["Extrema"](fld)[0]
        if mi < fieldrange[0]:
            fieldrange[0] = mi
        if ma > fieldrange[1]:
            fieldrange[1] = ma
        print 'current fieldrange =', fieldrange

# Increase the range by the buffer
fieldrange[0] -= fieldrange_buffer[0]
fieldrange[1] += fieldrange_buffer[1]

# Loop over all plotfiles to create an image
for plotf in pflist:
    pf = load(plotf)
    print plotf
    pf.h
    
    # If we are using a derived field, create the data for that field and
    # use it
    if use_derived_field:
        add_field(derived_field, function=_new_dfield, units='')
        fld = derived_field
    
    # Create rays with data yt can use from the user-specified coordinates
    # of each ray, then initialize our image
    yt_rays = []
    for i, ray in enumerate(rays):
        yt_rays.append(pf.h.ray(rays[i][0],rays[i][1]))

    fig = pl.figure()

    # Plot each ray using their lengths, data, and specified color.
    # This needs some work in generalizing to an arbitrary number of rays
    if len(rays) == 3:
        plot = pl.plot(yt_rays[0]['t']*raylen[0], yt_rays[0][fld], colors[0],
                       yt_rays[1]['t']*raylen[1], yt_rays[1][fld], colors[1],
                       yt_rays[2]['t']*raylen[2], yt_rays[2][fld], colors[2])
    elif len(rays) == 2:
        plot = pl.plot(yt_rays[0]['t']*raylen[0], yt_rays[0][fld], colors[0],
                       yt_rays[1]['t']*raylen[1], yt_rays[1][fld], colors[1])
    elif len(rays) == 1:
        plot = pl.plot(yt_rays[0]['t']*raylen[0], yt_rays[0][fld], colors[0])
    else: sys.exit('Up to 3 rays only')

    # Set up the limits, labels, title, line thickness, gridding, and image
    # size, then save a unique image based on the plotfile and field.
    pl.ylim(fieldrange[0], fieldrange[1])
    pl.xlim(x_limits[0], x_limits[1])
    pl.title(title)
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)
    pl.setp(plot, linewidth=linewidth)
    if background_grid: pl.grid(True)
    plt = pl.gcf()
    plt.set_size_inches(size_in_inches[0], size_in_inches[1])

    pl.savefig('ray%s.%s.png' %(plotf, fld))
    print 'Wrote ray%s.%s.png' %(plotf, fld)

    
