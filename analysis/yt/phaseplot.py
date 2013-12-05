from yt.mods import *
import sys
import os
import fnmatch
import string

zfield = 'magmom'
pfprefix='plt'
pflist = []
for p in os.listdir("."):
    if fnmatch.fnmatch(p, "%s*" % (pfprefix)) and not fnmatch.fnmatch(p,'*%s' % ('.yt')):
        pflist.append(p)
    
pflist.sort()
#pflist=['plt00650']#,'plt00400','plt00100']
for i, plotf in enumerate(pflist):
    pf = load(plotf)
    print plotf
    if i==0:
        fieldrange = [0,0]
        fieldrange[0],fieldrange[1] = pf.h.all_data().quantities["Extrema"](zfield)[0]
    mi, ma = pf.h.all_data().quantities["Extrema"](zfield)[0]
    if mi < fieldrange[0]:
        fieldrange[0] = mi
    if ma > fieldrange[1]:
        fieldrange[1] = ma
    print 'current fieldrange =', fieldrange



for plotf in pflist:
    pf = load(plotf)
    print plotf
    
    pf.h
    pf.field_info[zfield]._units = r"\rm{g}*\rm{cm}/\rm{s}"
    pf.field_info['Temp']._units = r"\rm{K}"

    def _newDen(field, data):
        return data['density']
    add_field("Density", function=_newDen, units=r"\rm{g}/\rm{cm}^{3}")

    pc = PlotCollection(pf, 'c')

    region = pf.h.region([5.12e9, 5.12e9, 5.12e9], [0.0, 0.0, 0.0],
                         [1.024e10, 1.024e10, 1.024e10])
    phase = pc.add_phase_object(region, ['Temp', 'Density',zfield],
                                weight=None, x_bins=250,
                                y_bins=250)
    phase.set_zlim(fieldrange[0],fieldrange[1])

    pc.save('test'+'%s.png' %plotf)
