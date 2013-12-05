from yt.mods import *
import sys
import os
import fnmatch
import string

#----User-defined variables and specifications----

# Information about the range, frequency, and name of your plotfiles.
# The prefix is everything before the plotfile number

pfprefix = 'plt'

# prefix of the image files you are creating. The plotfile name will be
# appended to this. If you want to do it a different way, you can modify
# the code for this located at the very end of the script.
imprefix = 'den_tempslice'

# Choose the field you are rendering
field = 'Temp'
slice_field = 'density'


# Camera specifications. 

# c is the center of the viewport--the focal point for the camera.
# This is in physical coordinates for the problem domain.
c = [5.12e9, 5.12e9, 5.12e9]

# L is the normal vector between camera position and the center.
# E.g. if the stars are on the x-y plane and we wanted to look
# from above along the negative z direction, we set L = [0, 0, -1]
L = [0.0, 0.4, -1.0]

# Wfrac is the width of the image as a fraction of the domain width.
# This just allows us to zoom in.
Wfrac = 0.67

#  Nvec supplies the dimensions of the output image file in pixels.
Nvec = (1024,1024)

# no_ghost is an optimization option (true=optimized) that can
# sometimes create artifacts at grid edges when used. The speed
# difference is only noticeable for large plotfiles, so it's usually
# fine to opt for a better image.
no_ghost = False


# Extra options

# simulated rotation based on the plotfile time and rotational
# frequency.  rot_vector is the vector we rotate around.  This is
# based on the data, not the current view.  For Omega > 0, we set this
# to be negative ([0, 0, -1]), because it is the camera moving.
rotation = True
rot_vector=[0.0,0.0,-1.0]

# an arbitrary wire-frame box you can use, for example, to enclose
# interesting things or easily show scale.  box_dim uses the
# "lower-left" and "upper-right" coordinates of the box (physical
# coordinates) to define it. box_col specifies the RGBA color of the
# box (alpha ranges from 0-1).

draw_box = False
box_dim = [[3.0e9, 2.0e9, 5.0e9], [4.5e9, 6.0e9, 10.0e9]]
box_col = [1.0, 1.0, 1.0, 0.11]




#----User-defined layers (isocontours)----

# Will your rendering use logarithmic values (best, especially for
# large field ranges)? If so, set this to true. If you do, remember
# that all of your values, widths, and ranges must be logarithmic!
take_log=True


# The grey_opacity option allows layers to be obscured more by other
# layers--good for opaque renderings.
grey_opacity = True
grey_opacity_slice = True


#--Single layers--
# [value, gaussian width, [colormap range], [alpha], 'colormap']

# alpha can be greater than 100--the values are normalized during
# the rendering.
# A list of avaliable colormaps can be found at
# http://yt-project.org/docs/dev/visualizing/colormaps/index.html


c1 = [np.log10(1.4e8), 0.003,[np.log10(4.0e7), np.log10(5.5e8)],
      [27.0], 'YlOrBr_r']
c2 = [np.log10(1.0e4), 0.001, [np.log10(9.0e3), np.log10(5.5e8)],
      [45.0], 'RdBu']


# Put your single layers in this list
single = [c1,c2]
single_s = []


#--Multiple layers--
# [number of layers, width of layers, min value, max value, [colormap
#  range], [alpha list], 'colormap']

# Same as above, but quickly creates multiple layers that are evenly
# spaced in the range from min val to max val. Width of layers is the
# same for all, and can be automatically calculated as 
# 0.001*(max_val - min_val)/(number_of_layers)
# by setting it to None.

l1 = [4, 0.0015, np.log10(3.0e7), np.log10(7.0e7),
      [np.log10(1.0e7), np.log10(9.0e7)], 12*na.ones(4,dtype='float64'),
      'summer']
l2 = [4, 0.0003, np.log10(7.5e7), np.log10(1.1e8),
      [np.log10(4.0e5), np.log10(5.5e8)], [40.0,50.0,60.0,70.0],
      'RdBu']
l3 = [35, 0.023, np.log10(5.0e-3), np.log10(1.5e7),
      [np.log10(5.0e-3), np.log10(1.5e7)], 18*na.ones(35,dtype='float64'),
      'jet']

# Put your multiple-layer sets in this list
multiple = [l1,l2]
multiple_s = [l3]
    



#----Everything below is automatic----


# find all plotfiles with the prefix pfprefix
pflist = []
for p in os.listdir("."):
    if fnmatch.fnmatch(p, "%s*" % (pfprefix)) and not fnmatch.fnmatch(p,'*%s' % ('.yt')):
        pflist.append(p)
    
pflist.sort()
#pflist = ['plt00100', 'plt00400','plt00650']

# Grab the rotational frequency from the first plotfile's job_info file.
if rotation:
    f = open(pflist[0] + '/job_info', 'r')
    rot_period = None
    for line in f:
        if string.find(line, "rotational_period") > 0:
            rot_period = float(string.split(line, "= ")[1])
            break

    f.close()

    if rot_period == None:
        sys.exit("ERROR: rotational_period not found")

    print "rotation period = ", rot_period
else: print "Not performing rotation.."

# Correct vectors
rot_vector = [rot_vector[0],rot_vector[1],-rot_vector[2]]
L = [-L[0],L[1],-L[2]]

# Loop over plotfiles and render us up some isocontours

for plotf in pflist:
    
    # Load plotfile
    pf = load(plotf)
    print plotf
    

    # Get min/max of field, convert to log if necessary
    dd = pf.h.all_data()
    mi, ma = dd.quantities['Extrema'](field)[0]
    tmi, tma = dd.quantities['Extrema'](slice_field)[0]
    W = Wfrac*(pf.domain_right_edge - pf.domain_left_edge)
    if take_log:
        mi, ma = np.log10(mi), np.log10(ma)
        tmi, tma = np.log10(tmi), np.log10(tma)
        pf.field_info[field].take_log=True
        pf.field_info[slice_field].take_log=True
    

    # Initialize transfer function
    tf = ColorTransferFunction((mi-1, ma+1), nbins=8.0e5)
    tf2 = ColorTransferFunction((tmi-1, tma+1), nbins=8.0e5)

    # Initialize  all layers
    for ele in single:
        tf.sample_colormap(ele[0], ele[1], col_bounds=ele[2],
                           alpha=ele[3], colormap=ele[4])
        
    for elem in multiple:
        tf.add_layers(elem[0], w=elem[1], mi=elem[2], ma=elem[3],
                      col_bounds=elem[4], alpha=elem[5], colormap=elem[6])

    for ele in single_s:
        tf2.sample_colormap(ele[0], ele[1], col_bounds=ele[2],
                           alpha=ele[3], colormap=ele[4])
        
    for elem in multiple_s:
        tf2.add_layers(elem[0], w=elem[1], mi=elem[2], ma=elem[3],
                      col_bounds=elem[4], alpha=elem[5], colormap=elem[6])


    # Assign grey_opacity
    tf.grey_opacity=grey_opacity
    tf2.grey_opacity=grey_opacity_slice


    # Initialize camera object
    
    normal_vector = L
    normal_vector /= np.sqrt( np.dot(normal_vector, normal_vector))
    vecs = np.identity(3)
    tbv = np.cross(normal_vector, vecs).sum(axis=1)
    ax = tbv.argmax()
    east_vector = np.cross(vecs[2,:], normal_vector).ravel()
    while np.dot(east_vector,east_vector) == 0.0:
        axnum = 1
        east_vector = np.cross(vecs[axnum,:], normal_vector).ravel()
        axnum -= 1
    north_vector = np.cross(normal_vector, east_vector).ravel()
    north_vector = -north_vector



    cam = pf.h.camera(c, L, W, Nvec, transfer_function=tf,
                      north_vector=north_vector,
                      fields=[field], pf=pf, no_ghost=no_ghost)

    # Add rotation if necessary
    if rotation:
        cam.rotate(2*np.pi*(pf.current_time/rot_period),
                   rot_vector=rot_vector)

    # Render the isocontours
    im = cam.snapshot()

    # Add a box if necessary
    if draw_box:
        cam.draw_box(im, np.array(box_dim[0]), np.array(box_dim[1]),
                     np.array(box_col))

    # Write out a .png image for that plotfile
    im.write_png(imprefix  + '.vol.%s.png' %plotf, background=None)


    # Do the same for the 'slice'
    cam = pf.h.camera(c, L, W, Nvec, transfer_function=tf2,
                      north_vector=north_vector,
                      fields=[slice_field], pf=pf, no_ghost=no_ghost,
                      le=[0.0,0.0,5.0e9], re=[1.024e10,1.024e10,5.24e9])
    if rotation:
        cam.rotate(2*np.pi*(pf.current_time/rot_period),
                   rot_vector=rot_vector)
    im = cam.snapshot()
    if draw_box:
        cam.draw_box(im, np.array(box_dim[0]), np.array(box_dim[1]),
                     np.array(box_col))
    im.write_png(imprefix  + '.slice.%s.png' %plotf)

    os.system('composite '+imprefix+'.vol.{a}.png '.format(a=plotf)+imprefix+'.slice.{b}.png '.format(b=plotf)+imprefix+'{c}.png'.format(c=plotf))

    tf.clear()
    tf2.clear()
