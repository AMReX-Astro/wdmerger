import os
import numpy as np


#
# Return the name of the latest wdmerger output file in the directory dir.
#

def get_last_output(dir):

    # Open up the standard output for analysis. It will be the numerically last file
    # starting with the designated output string.

    files = os.listdir(dir)

    files = sorted(filter(lambda s: s[0:9] == "wdmerger.",files))

    if (len(files) == 0):
        print "Error: No wdmerger output files in directory " + dir
        exit()

    return "dir" + "/" + files[len(files)-1]


#
# Get the timing data for a wdmerger output file.
#

def timing(output_filename):

    # Read in the file and return the following data in arrays:
    # - The coarse timestep time for each step
    # - The plotfile time for each step
    # For the sake of generality, we will count up how many
    # steps we have taken after the fact.

    output = open(output_filename, 'r')
    lines = output.readlines()
    coarseSteps = filter(lambda s: s[0:6] == "Coarse",lines)
    coarseSteps = [float(s.split()[3]) for s in coarseSteps] # Extract out the time only

    avg_timestep = np.average(coarseSteps)

    # Now subtract out the gravity solve time in each case

    grav_time = filter(lambda s: s[0:7] == "Gravity",lines)
    grav_time = [float(s.split()[3]) for s in grav_time]

    # Remove the first two, since they are related to the initial multilevel solve

    grav_time = grav_time[2:]

    # For each coarse timestep, there are two BC calculations and two Poisson solves. Let's
    # sum these for each timestep. For each refined level, there's two Poisson solves per subcycle.
    # For the coarse grid, there's two Poisson solves and also two BC fills.
    # Therefore, the number of gravity calculations per timestep is equal to
    # (nlevs - 1) * ref_ratio + 4

    ref_ratio = 4
    nlevs = 2

    grav_per_timestep = (nlevs - 1) * 2 * ref_ratio + 4

    for n in range(len(coarseSteps)):

        for i in range(grav_per_timestep):
            coarseSteps[n] -= float(grav_time[n+i])

    avg_timestep_no_grav = np.average(coarseSteps)

    return [avg_timestep, avg_timestep_no_grav]



#
# Get a sorted list of all the plotfiles in directory dir.
#

def get_plotfiles(dir):

    # Check to make sure the directory exists.

    if (not os.path.isdir(dir)):
        print "Error: Directory " + dir + " does not exist, exiting."
        exit()

    # List all of the files in the directory.

    dir_contents = os.listdir(dir)

    # Strip out non-plotfiles.

    plotfiles = sorted(filter(lambda s: s[0:3] == 'plt', dir_contents))

    # Remove any redundant plotfiles. These are generated with names like "old" or "temp"
    # and an easy way to get rid of them is just to check on all plotfiles that don't have
    # the right number of characters.

    plotfiles = filter(lambda s: len(s) == 8 or len(s) == 9, plotfiles)

    # This step strips out 'plt' from each plotfile, sorts them numerically,
    # then recasts them as strings. This is necessary because we may have some
    # plotfiles with six digits and some with five, and the default sort will
    # not sort them properly numerically.

    plotfiles = ['plt' + str(y).zfill(5) for y in sorted([int(x[3:]) for x in plotfiles])]

