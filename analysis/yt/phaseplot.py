from yt.mods import *
import sys
import os
import fnmatch
import string

#type in path to directory for data; 
#.png images will also be saved there by default
#otherwise change 'outpath'
#PS: don't forget to include '/' at the end of path

path = '/home/'
outpath = path
pfprefix='plt'

#the list of all fiels in the form ('a', 'b') can be found
#by command pf.derived_field_list after data was loaded. 

xfield = ('boxlib', 'Temp')
yfield = ('boxlib', 'density')
zfield = ('boxlib', 'magmom')

#by default region includes all data
#if you want to plot specific region, then
#set UseRegion to 'True' and input your region for each direction

UseRegion = False
xregion = [5.12e9, 5.12e9, 5.12e9]
yregion = [0.0, 0.0, 0.0]
zregion = [1.024e10, 1.024e10, 1.024e10]

#also you can change the accuracy of your plot
#by changing x_bins and y_bins.
#higher value -> better accuracy + longer to compute

x_bins = 256
y_bins = 256

#you can choose a weight field for calculating weighted averages
#if None, the profile values are the sum of the field values within the bin

weight_field = None


#----------Automatic----------

fields = [xfield, yfield, zfield]
pflist = []

#Find all data files
for p in os.listdir(path):
    if fnmatch.fnmatch(p, "%s*" % (pfprefix)) and not fnmatch.fnmatch(p,'*%s' % ('.yt')):
        pflist.append(p)
pflist.sort()

#Write out the ranges for zfield    
for i, plotf in enumerate(pflist):
    pf = load(path + plotf)
    print plotf
    if i==0:
        fieldrange = [0,0]
        fieldrange[0],fieldrange[1] = pf.h.all_data().quantities["Extrema"](zfield)
    mi, ma = pf.h.all_data().quantities["Extrema"](zfield)
    if mi < fieldrange[0]:
        fieldrange[0] = mi
    if ma > fieldrange[1]:
        fieldrange[1] = ma
    print 'current fieldrange =', fieldrange

#Load and plot the data
for plotf in pflist:
    pf = load(path + plotf)
    print plotf
    skip = False

    #Checking if fields contain only 0 values; cannot be plotted    
    for field in fields:
        mi, ma = pf.h.all_data().quantities["Extrema"](field)
        if '0.0' in str(ma):
            skip = True
            print str(field) + ' in ' + plotf + ' contains only 0 values'
            
    if skip:
        print '\n---> Skipping ' + plotf + '\n'
        continue

    #Choosing a region
    if UseRegion:
        region = pf.region(xregion,yregion,zregion)
    else:
        region = pf.all_data()

    #Plotting
    phase = PhasePlot(region, xfield, yfield, zfield, weight_field=weight_field, 
		      x_bins=x_bins, y_bins=y_bins)
    
    #phase.set_zlim(fieldrange[0],fieldrange[1])

    phase.save(outpath + 'phase_'+'%s.png' %plotf)
