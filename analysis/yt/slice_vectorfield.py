from yt.mods import *
import sys
import os
import fnmatch
import string

# Type in path to directory for data; 
# .png images will also be saved there by default.
# Otherwise change 'outpath'.
# By default we will assume that you are running this in the directory with the data.

path = './'
outpath = path
pfprefix='plt'

# The list of all fields in the form ('a', 'b') can be found
# by command pf.derived_field_list after data was loaded.

field = ('gas', 'density')

#for log plot

log = True

#----------Automatic----------

pflist = []
for p in os.listdir(path):
    if fnmatch.fnmatch(p, "%s*" % (pfprefix)) and not fnmatch.fnmatch(p,'*%s' % ('.yt')):
        pflist.append(p)
    
pflist.sort()

# Find min and max for colorbar, eventually expand to emulate generalline.py
for i, plotf in enumerate(pflist):
    pf = load(path + plotf)
    print plotf
    if i==0:
        fieldrange = [0,0]
        fieldrange[0],fieldrange[1] = pf.h.all_data().quantities["Extrema"](field)
    mi, ma = pf.h.all_data().quantities["Extrema"](field)
    if mi < fieldrange[0]:
        fieldrange[0] = mi
    if ma > fieldrange[1]:
        fieldrange[1] = ma
    print 'current fieldrange =', fieldrange

for plotf in pflist:
    pf = load(path + plotf)
    print plotf

    pf.field_info[field].take_log = log

    p = SlicePlot(pf, "z", field, width=((1.024e10, "cm"), (1.024e10, "cm")),
                  fontsize=14)
    p.set_zlim(field, fieldrange[0],fieldrange[1])
    p.set_cmap(field, 'hot')

    #p.annotate_contour('density', clim=(1.05e-4, 1.16e-4), ncont=5, label=False)
    p.annotate_grids(alpha=0.2, min_level=2)
    p.annotate_quiver('x_velocity', 'y_velocity', factor=16)
    p.annotate_title('Log density with velocity field')
    p.annotate_point([5.0e7, 1.0e8],
                     'Current time: {a} s'.format(a=pf.current_time),
                     text_args={'size':'large', 'color':'w'})
    #p.annotate_marker([5.0e9, 5.0e9], marker="x")
    #p.annotate_point([5.7e9, 5.1e9], 'pfpf')

    ##These guys are ignored by default in the code so we can't use them. Dunno why.
    ##Can find it in yt-hg/yt/visualization/plot_window.py, line 705.
    #p.annotate_units(unit="g/cm*s^2", factor = 10, text_annotate=True,
     #                text_which=1)
    #p.annotate_axis_label("rtw")

    p.save('{outpath}{field[1]}_velocityfield.{plotf}.png'.format(outpath = outpath, field = field, plotf = plotf))
