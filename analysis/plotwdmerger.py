#!/usr/bin/env python

# a simple script to plot 2-d or 3-d BoxLib data using the matplotlib
# library
#
# 2011-12-02 M. Zingale

import fsnapshot
import numpy
import pylab
import matplotlib
import os
import sys
import getopt
import math
import string
import mpl_toolkits.axes_grid1


#==============================================================================
# do_plot
#==============================================================================
def do_plot(plotfile, component, outFile, log, 
            minval, maxval, eps, dpi, origin, annotation, 
            xmin_pass, ymin_pass, zmin_pass, 
            xmax_pass, ymax_pass, zmax_pass):


    #--------------------------------------------------------------------------
    # construct the output file name
    #--------------------------------------------------------------------------
    if (outFile == ""):
        outFile = os.path.normpath(plotfile) + "_" + component

        if (not eps):
            outFile += ".png"

        else:
            outFile += ".eps"

    else:
        # make sure the proper extension is used
        if (not eps):
            if (not string.rfind(outFile, ".png") > 0):
                outFile = outFile + ".png"

        else:
            if (not string.rfind(outFile, ".eps") > 0):
                outFile = outFile + ".eps"


    #--------------------------------------------------------------------------
    # read in the meta-data from the plotfile
    #--------------------------------------------------------------------------
    (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotfile)

    time = fsnapshot.fplotfile_get_time(plotfile)

    (xmin, xmax, ymin, ymax, zmin, zmax) = \
        fsnapshot.fplotfile_get_limits(plotfile)

    dx = (xmax - xmin)/nx
    x = xmin + numpy.arange( (nx), dtype=numpy.float64 )*dx

    dy = (ymax - ymin)/ny
    y = ymin + numpy.arange( (ny), dtype=numpy.float64 )*dy

    if (nz > 0):
        dz = (zmax - zmin)/nz
        z = zmin + numpy.arange( (nz), dtype=numpy.float64 )*dz


    #----------------------------------------------------------------------
    # 3-d plot
    #----------------------------------------------------------------------
        
    # starting points for the figure positions

    # assume that the width of the plotting area is 0.05 to 0.95,
    # leaving 0.9 to be split amongst the 3 plots.  So each has a
    # width of 0.3

    # for the height, we will assume that the colorbar at the
    # bottom gets 0.15, and that we go until 0.95, leaving 0.8 of
    # height for the plots.
            
    pos1 = [0.05, 0.15, 0.3, 0.8]
    pos2 = [0.35, 0.15, 0.3, 0.8]
    pos3 = [0.65, 0.15, 0.3, 0.8]

    fig = pylab.figure()


    # read in the slices
    # x-y
    data_xy = numpy.zeros( (nx, ny), dtype=numpy.float64)

    indir = 3
    (data_xy, err) = \
        fsnapshot.fplotfile_get_data_3d(plotfile, component, indir, 
                                        origin, data_xy)
    if (not err == 0):
        sys.exit(2)

    data_xy = numpy.transpose(data_xy)

    if log:
        if (numpy.min(data_xy) < 0):
            data_xy = numpy.log10(numpy.abs(data_xy))
        else:
            data_xy = numpy.log10(data_xy)




                
    if (not minval == None): 
        if (log):
            minval = math.log10(minval)
    else:
        minval = numpy.min(data_xy)


    if (not maxval == None): 
        if (log):
            maxval = math.log10(maxval)
    else:
        maxval = numpy.max(data_xy)



    # x-y
    extent = [xmin, xmax, ymin, ymax]

    if (not xmin_pass == None):
        extent[0] = xmin_pass

    if (not ymin_pass == None):
        extent[2] = ymin_pass

    if (not xmax_pass == None):
        extent[1] = xmax_pass

    if (not ymax_pass == None):
        extent[3] = ymax_pass

    ix0 = 0
    if (not xmin_pass == None):
        ix0 = int((xmin_pass - xmin)/dx)

    iy0 = 0
    if (not ymin_pass == None):
        iy0 = int((ymin_pass - ymin)/dy)

    ix = nx
    if (not xmax_pass == None):
        ix = int((xmax_pass - xmin)/dx)

    iy = ny
    if (not ymax_pass == None):
        iy = int((ymax_pass - ymin)/dy)



    ax = pylab.subplot(1,1,1)
    #pylab.subplots_adjust(wspace=0.4)

    im=pylab.imshow(data_xy[iy0:iy,ix0:ix],origin='lower', extent=extent, 
                    vmin=minval, vmax=maxval, axes=pos1)

    pylab.xlabel("x")
    pylab.ylabel("y")

    # axis labels in scientific notation with LaTeX
    ax = pylab.gca()
    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    ax.xaxis.offsetText.set_fontsize("small")
    ax.yaxis.offsetText.set_fontsize("small")

    cl = pylab.getp(ax, 'ymajorticklabels')
    pylab.setp(cl, fontsize=10)
    cl = pylab.getp(ax, 'xmajorticklabels')
    pylab.setp(cl, fontsize=10)

    # draw a line through the middle
    pylab.plot([extent[0], extent[1]], [0.5*(ymax + ymin), 0.5*(ymax+ymin)], c="w")


    # do a fixed offset in pixels from the (xmin,ymin) data point
    fig1 = ax.get_figure()
    trans=matplotlib.transforms.offset_copy(ax.transData, x=0, y=-0.5, 
                                            fig=fig1, units='inches')

    pylab.text(xmin, ymin, "time = %7.3g s" % (time), 
               verticalalignment="bottom", transform = trans, 
               clip_on=False, fontsize=10)




    pylab.axis(extent)

    # colorbar
    pylab.subplots_adjust(bottom=0.1, left=0.05, right=0.95)
    
    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
        
    #cax = pylab.axes([0.05, 0.06, 0.9, 0.04])
    pylab.colorbar()#orientation="horizontal", cax=cax, format=formatter)

    pylab.title(component)


    #--------------------------------------------------------------------------
    # save the figure
    #--------------------------------------------------------------------------
    if (not eps):
        pylab.savefig(outFile, bbox_inches='tight', dpi=dpi, pad_inches=0.33)
    else:
        pylab.savefig(outFile, bbox_inches='tight', pad_inches=0.33)



#==============================================================================
# usage
#==============================================================================
def usage():
    usageStr = """
    hacked version of plotsinglevar.py for making wdmerger comparison plots"""

    print usageStr



#==============================================================================
# main
#==============================================================================
if __name__== "__main__":

    outFile = ""
    log = 0
    eps = 0
    minvar = None
    maxvar = None
    dpi = 100
    origin = 0
    annotation = ""
    xmax = None
    ymax = None
    zmax = None
    xmin = None
    ymin = None
    zmin = None

    try: opts, next = getopt.getopt(sys.argv[1:], "o:m:M:x:X:y:Y:z:Z:", 
                                    ["log","eps","dpi=","origin","annotate="])
    except getopt.GetoptError:
        print "invalid calling sequence"
        usage()
        sys.exit(2) 
               

    for o, a in opts:

        if o == "-o":
            outFile = a

        if o == "-m":
            try: minvar = float(a)
            except ValueError:
                print "invalid value for -m"
                sys.exit(2)

        if o == "-M":
            try: maxvar = float(a)
            except ValueError:
                print "invalid value for -M"
                sys.exit(2)

        if o == "-X":
            try: xmax = float(a)
            except ValueError:
                print "invalid value for -X"
                sys.exit(2)            

        if o == "-Y":
            try: ymax = float(a)
            except ValueError:
                print "invalid value for -Y"
                sys.exit(2)            

        if o == "-Z":
            try: zmax = float(a)
            except ValueError:
                print "invalid value for -Z"
                sys.exit(2)            

        if o == "-x":
            try: xmin = float(a)
            except ValueError:
                print "invalid value for -x"
                sys.exit(2)            

        if o == "-y":
            try: ymin = float(a)
            except ValueError:
                print "invalid value for -y"
                sys.exit(2)            

        if o == "-z":
            try: zmin = float(a)
            except ValueError:
                print "invalid value for -z"
                sys.exit(2)            
 
        if o == "--log":
            log = 1

        if o == "--eps":
            eps = 1

        if o == "--dpi":
            try: dpi = int(a)
            except ValueError:
                print "invalid value for --dpi"
                sys.exit(2)

        if o == "--origin":
            origin = 1

        if o == "--annotate":
            annotation = a



    try: plotfile = next[0]
    except IndexError:
        print "ERROR: plotfile not specified"
        usage()
        sys.exit(2)

    try: component = next[1]
    except IndexError:
        print "ERROR: no component specified"
        usage()
        sys.exit(2)    

    do_plot(plotfile, component, outFile, 
            log, minvar, maxvar, eps, dpi, origin, annotation,
            xmin, ymin, zmin, xmax, ymax, zmax)
