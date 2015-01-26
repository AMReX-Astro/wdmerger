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
imprefix = 'generic2'

# Choose the field you are rendering
field = 'density'


# Camera specifications. c is the center of the viewport--the focal
# point for the camera. L is the normal vector between camera position
# and the center. Wfrac is the width of the image as a fraction of the
# domain width. Nvec supplies the dimensions of your image file in pixels.
# no_ghost is an optimization option (true=optimized) that can sometimes
# create artifacts at grid edges when used. The speed difference is only
# noticeable for large plotfiles, so it's usually fine to opt for a better
# image. north_vector defines an 'up' direction to relate L to.
# Note: yt vectors are automatically normalized and can be picky--you
# may have to play around to get what you want.

c = [5.0e9, 5.0e9, 5.0e9]
L = [0.15, 0.4, 1.0]
Wfrac = 0.67
Nvec = (1024,1024)
no_ghost = False
north_vector=[0.0, -1.0, 0.0]


# Extra options: simulated rotation based on the plotfile time and
# rotational frequency, and an arbitrary wire-frame box you can use,
# for example, to enclose interesting things or easily show scale.
# rot_vector is the vector we rotate around. box_dim uses the "lower-left"
# and "upper-right" coordinates of the box to define it. box_col
# specifies the RGBA color of the box (alpha ranges from 0-1).

rotation = True
rot_vector=[0.0,1.0,0.0]
draw_box = True
box_dim = [[1.6e9, 3.25e9, 3.25e9], [8.4e9, 6.75e9, 6.75e9]]
box_col = [1.0, 1.0, 1.0, 0.11]




#----User-defined layers (isocontours)----

# Will your rendering use logarithmic values (best, especially for
# large field ranges)? If so, set this to true. If you do, remember
# that all of your values, widths, and ranges must be logarithmic!
take_log=True


# The grey_opacity option allows layers to be obscured more by other
# layers--good for opaque renderings.
grey_opacity = True


#--Single layers--
# [value, gaussian width, [colormap range], [alpha], 'colormap']

# alpha can be greater than 100--the values are normalized during
# the rendering.
# A list of avaliable colormaps can be found at
# http://yt-project.org/docs/dev/visualizing/colormaps/index.html

c1 = [np.log10(1.6e-1), 0.008, [np.log10(7.0e-2), np.log10(1.1e0)],
      [13.0], 'YlOrBr']
c2 = [np.log10(1.0e0), 0.005, [np.log10(7.0e-2), np.log10(1.8e0)],
      [26.0], 'YlOrBr']
c3 = [np.log10(1.0e1), 0.008, [np.log10(9.0e0), np.log10(1.0e6)],
      [380.0], 'RdBu_r']
c4 = [np.log10(1.0e2), 0.028, [np.log10(8.0e1), np.log10(1.0e7)],
      [100.0], 'ocean']

# Put your single layers in this list
single = [c1,c2,c3,c4]


#--Multiple layers--
# [number of layers, width of layers, min value, max value, [colormap
#  range], [alpha list], 'colormap']

# Same as above, but quickly creates multiple layers that are evenly
# spaced in the range from min val to max val. Width of layers is the
# same for all, and can be automatically calculated as 
# 0.001*(max_val - min_val)/(number_of_layers)
# by setting it to None.

l1 = [3, 0.0003, np.log10(1.8e-2), np.log10(1.0e-1),
      [np.log10(1.6e-2), np.log10(1.4e-1)], 13*na.ones(4,dtype='float64'),
      'summer']
l2 = [4, 0.014, np.log10(1.0e2), np.log10(1.7e7),
      [np.log10(2.0e1), np.log10(1.75e7)], reversed(np.logspace(1.0,1.15,4)),
      'RdBu']

# Put your multiple-layer sets in this list
multiple = [l1,l2]
    



#----Everything below is automatic----


# find all plotfiles with the prefix pfprefix
pflist = []
for p in os.listdir("."):
    if fnmatch.fnmatch(p, "%s*" % (pfprefix)):
        pflist.append(p)
    
pflist.sort()


# Grab the rotational frequency from the first plotfile's job_info file.
f = open(pflist[0] + '/job_info', 'r')
rot_freq = None
for line in f:
    
    if string.find(line, "rotational_freq") > 0:
        rot_freq = float(string.split(line, "=")[1])
        break

f.close()

if rot_freq == None:
    sys.exit("ERROR: rotational_frequency not found")

print "rotation frequency = ", rot_freq


# Loop over plotfiles and render us up some isocontours

for plotf in pflist:
    
    # Load plotfile
    pf = load(plotf)
    print plotf
    

    # Get min/max of field, convert to log if necessary
    dd = pf.h.all_data()
    mi, ma = dd.quantities['Extrema'](field)[0]
    W = Wfrac*(pf.domain_right_edge - pf.domain_left_edge)
    if take_log:
        mi, ma = np.log10(mi), np.log10(ma)
        pf.field_info[field].take_log=True
    

    # Initialize transfer function
    tf = ColorTransferFunction((mi-1, ma+1), nbins=8.0e5)

    # Initialize  all layers
    for i in single:
        tf.sample_colormap(i[0], i[1], col_bounds=i[2],
                           alpha=i[3], colormap=i[4])
        
    for j in multiple:
        tf.add_layers(j[0], w=j[1], mi=j[2], ma=j[3],
                      col_bounds=j[4], alpha=j[5], colormap=j[6])

    # Assign grey_opacity
    tf.grey_opacity=grey_opacity
    

    # Initialize camera object
    cam = pf.h.camera(c, L, W, Nvec, transfer_function = tf,
                      north_vector=north_vector,
                      fields=[field], pf=pf, no_ghost=no_ghost)

    # Add rotation if necessary
    if rotation:
        cam.rotate(2*np.pi*(pf.current_time*rot_freq),
                   rot_vector=rot_vector)

    # Render the isocontours
    im = cam.snapshot()
    im.add_background_color('black', inline=True)

    # Add a box if necessary
    if draw_box:
        cam.draw_box(im, np.array(box_dim[0]), np.array(box_dim[1]),
                     np.array(box_col))

    # Write out a .png image for that plotfile
    im.write_png(imprefix  + '%s.png' %plotf)

    tf.clear()
