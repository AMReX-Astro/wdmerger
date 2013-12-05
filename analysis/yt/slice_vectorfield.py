from yt.mods import *
import sys
import os
import fnmatch
import string

pfprefix='plt'
pflist = []
for p in os.listdir("."):
    if fnmatch.fnmatch(p, "%s*" % (pfprefix)) and not fnmatch.fnmatch(p,'*%s' % ('.yt')):
        pflist.append(p)
    
pflist.sort()
#pflist=['plt00650']#,'plt00400','plt00100']

field = 'logden'

# Find min and max for colorbar.. eventually expand to emulate generalline.py
for i, plotf in enumerate(pflist):
    pf = load(plotf)
    print plotf
    if i==0:
        fieldrange = [0,0]
        fieldrange[0],fieldrange[1] = pf.h.all_data().quantities["Extrema"](field)[0]
    mi, ma = pf.h.all_data().quantities["Extrema"](field)[0]
    if mi < fieldrange[0]:
        fieldrange[0] = mi
    if ma > fieldrange[1]:
        fieldrange[1] = ma
    print 'current fieldrange =', fieldrange

for plotf in pflist:
    pf = load(plotf)
    print plotf
    pf.h

    pf.field_info[field]._units = r"\rm{log}[\rm{g}/\rm{cm}^{3}]"

    p = SlicePlot(pf, "z", field, width=((1.024e10, "cm"), (1.024e10, "cm")),
                  fontsize=14)
    p.set_zlim(field, fieldrange[0],fieldrange[1])
    p.set_cmap(field, 'jet')

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

    p.save('{a}_velocityfield.{b}.png'.format(a=field,b=plotf))
