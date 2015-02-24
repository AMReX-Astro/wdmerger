import os
import numpy as np
import yt
from yt.analysis_modules.level_sets.api import *
import string

#
# Returns the current git hash of the wdmerger repo.
# Credit: http://stackoverflow.com/questions/14989858/get-the-current-git-hash-in-a-python-script
#

def get_wdmerger_git_commit_hash():
    import subprocess
    cwd = os.getcwd()
    WDMERGER_HOME = os.getenv('WDMERGER_HOME')
    os.chdir(WDMERGER_HOME)
    hash = subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip()
    os.chdir(cwd)

    return hash



#
# Return the CASTRO directory on the current machine, based on your environment variables.
#

def get_castro_dir():
    CASTRO_DIR = os.getenv('CASTRO_DIR')
    return CASTRO_DIR



#
# Given a plotfile directory, return the CASTRO and BoxLib git commit hashes.
#

def get_git_commits_from_plotfile(plotfile):
    job_info = open(plotfile + "/job_info", 'r')
    lines = job_info.readlines()
    lines = [line.split() for line in lines]
    castro_hash   = ""
    boxlib_hash   = ""
    wdmerger_hash = ""
    for line in lines:
        if (len(line) == 4):
            if (line[0] == "Castro" and line[1] == "git" and line[2] == "hash:"):
                castro_hash = line[3]
            elif (line[0] == "BoxLib" and line[1] == "git" and line[2] == "hash:"):
                boxlib_hash = line[3]
            elif (line[0] == "wdmerger" and line[1] == "git" and line[2] == "hash:"):
                wdmerger_hash = line[3]

    job_info.close()

    return [castro_hash, boxlib_hash, wdmerger_hash]



#
# Given a diagnostic output file, return the CASTRO and BoxLib git commit hashes.
#

def get_git_commits_from_diagfile(diagfile):
    diagfile = open(diagfile, 'r')

    castro_hash   = ""
    boxlib_hash   = ""
    wdmerger_hash = ""

    line = diagfile.readline().split()

    if (line[0] == "Castro" and line[1] == "git" and line[2] == "hash:"):
        castro_hash = line[4]

    line = diagfile.readline().split()

    if (line[0] == "BoxLib" and line[1] == "git" and line[2] == "hash:"):
        boxlib_hash = line[4]

    line = diagfile.readline().split()

    if (line[0] == "wdmerger" and line[1] == "git" and line[2] == "hash:"):
        wdmerger_hash = line[4]

    diagfile.close()

    return [castro_hash, boxlib_hash, wdmerger_hash]



#
# Given the stdout from a Castro run, return the CASTRO and BoxLib git commit hashes.
#

def get_git_commits_from_infofile(infofile):
    infofile = open(infofile, 'r')
    lines = infofile.readlines()
    lines = [line.split() for line in lines]
    castro_hash   = ""
    boxlib_hash   = ""
    wdmerger_hash = ""
    for line in lines:
        if (len(line) == 4):
            if (line[0] == "Castro" and line[1] == "git" and line[2] == "hash:"):
                castro_hash = line[3]
            elif (line[0] == "BoxLib" and line[1] == "git" and line[2] == "hash:"):
                boxlib_hash = line[3]
            elif (line[0] == "wdmerger" and line[1] == "git" and line[2] == "hash:"):
                wdmerger_hash = line[3]

    infofile.close()

    return [castro_hash, boxlib_hash, wdmerger_hash]



#
# Given CASTRO and BoxLib hashes that were used to create the plot for a given plotfile,
# insert these and the current wdmerger hash into the EPS file.
# Credit: http://stackoverflow.com/questions/1325905/inserting-line-at-specified-position-of-a-text-file-in-python
# A comma in a print statement effectively prevents a newline from being appended to the print. This is valid in Python 2.x,
# but is not valid in Python 3.x. Instead, one should use the print() function in that case.
# Source: http://stackoverflow.com/questions/11266068/python-avoid-new-line-with-print-command
#

def insert_commits_into_eps(eps_file, data_file, data_file_type):
    import fileinput

    if (data_file_type == 'plot'):
        [castro_hash, boxlib_hash, wdmerger_hash] = get_git_commits_from_plotfile(data_file)
    elif (data_file_type == 'diag'):
        [castro_hash, boxlib_hash, wdmerger_hash] = get_git_commits_from_diagfile(data_file)
    elif (data_file_type == 'info'):
        [castro_hash, boxlib_hash, wdmerger_hash] = get_git_commits_from_infofile(data_file)
    else:
        print "Error: Data file type not recognized."

    input = fileinput.input(eps_file, inplace=True)

    for line in input:
        print line,
        if line.startswith('%%CreationDate:'):
            print "%%CASTRO git hash: " + castro_hash + "\n" + \
                  "%%BoxLib git hash: " + boxlib_hash + "\n" + \
                  "%%wdmerger git hash: " + wdmerger_hash



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

    return dir + "/" + files[len(files)-1]


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
# Add up all the energy and momentum losses reported from castro.print_energy_diagnostics.
#

def energy_momentum_diagnostics(output_filename):

    # Read in the file and return the following data in arrays:
    # - All energy added from the gravitational source terms
    # - All energy added from the gravitational correction terms
    # - All momentum added from the gravitational correction terms

    output = open(output_filename, 'r')
    lines = output.readlines()
    rho_E_lines = filter(lambda s: s[0:7] == "(rho E)",lines)
    grav_E_lines = filter(lambda s: s.split()[4] == "grav.",rho_E_lines)
    rot_E_lines = filter(lambda s: s.split()[4] == "rot.",rho_E_lines)
    flux_E_lines    = filter(lambda s: s.split()[4] == "fluxes",rho_E_lines)
    mass_lines = filter(lambda s: s[0:7] == "   Mass",lines)
    xmom_lines = filter(lambda s: s[0:4] == "xmom",lines)
    ymom_lines = filter(lambda s: s[0:4] == "ymom",lines)
    zmom_lines = filter(lambda s: s[0:4] == "zmom",lines)

    neg_rho_m_lines = filter(lambda s: s.split()[3] == "negative",mass_lines)
    neg_rho_E_lines = filter(lambda s: s.split()[4] == "negative",mass_lines)

    rot_xmom_lines = filter(lambda s: s.split()[3] == "rot.",xmom_lines)
    rot_ymom_lines = filter(lambda s: s.split()[3] == "rot.",ymom_lines)
    rot_zmom_lines = filter(lambda s: s.split()[3] == "rot.",zmom_lines)
    grav_xmom_lines = filter(lambda s: s.split()[3] == "grav.",xmom_lines)
    grav_ymom_lines = filter(lambda s: s.split()[3] == "grav.",ymom_lines)
    grav_zmom_lines = filter(lambda s: s.split()[3] == "grav.",zmom_lines)
    flux_xmom_lines = filter(lambda s: s.split()[3] == "fluxes",xmom_lines)
    flux_ymom_lines = filter(lambda s: s.split()[3] == "fluxes",ymom_lines)
    flux_zmom_lines = filter(lambda s: s.split()[3] == "fluxes",zmom_lines)
    reset_E_lines   = filter(lambda s: s.split()[4] == "reset",rho_E_lines)

    mass_added_neg_reset = sum([float(s.split()[7]) for s in neg_rho_m_lines])
    E_added_neg_reset    = sum([float(s.split()[8]) for s in neg_rho_E_lines])

    E_added_flux = sum([float(s.split()[6]) for s in flux_E_lines])
    E_added_grav = sum([float(s.split()[8]) for s in grav_E_lines])
    E_added_rot  = sum([float(s.split()[8]) for s in rot_E_lines])

    xmom_added_flux = sum([float(s.split()[5]) for s in flux_xmom_lines])
    ymom_added_flux = sum([float(s.split()[5]) for s in flux_ymom_lines])
    zmom_added_flux = sum([float(s.split()[5]) for s in flux_zmom_lines])
    xmom_added_grav = sum([float(s.split()[7]) for s in grav_xmom_lines])
    ymom_added_grav = sum([float(s.split()[7]) for s in grav_ymom_lines])
    zmom_added_grav = sum([float(s.split()[7]) for s in grav_zmom_lines])
    xmom_added_rot = sum([float(s.split()[7]) for s in rot_xmom_lines])
    ymom_added_rot = sum([float(s.split()[7]) for s in rot_ymom_lines])
    zmom_added_rot = sum([float(s.split()[7]) for s in rot_zmom_lines])

    E_added_reset   = sum([float(s.split()[7]) for s in reset_E_lines])

    print ""
    print "Analysis of output file " + output_filename
    print ""

    print "Mass added from negative density resets = " + str(mass_added_neg_reset)
    print "Energy added from negative density resets = " + str(E_added_neg_reset)
    print "Energy added from gravitational sources = " + str(E_added_grav)
    print "Energy added from rotation sources = " + str(E_added_rot)
    print "Energy added from hydro fluxes = " + str(E_added_flux)
    print "Energy added from resets = " + str(E_added_reset)
    print "xmom added from hydro fluxes = " + str(xmom_added_flux)
    print "ymom added from hydro fluxes = " + str(ymom_added_flux)
    print "zmom added from hydro fluxes = " + str(zmom_added_flux)
    print "xmom added from gravitational sources = " + str(xmom_added_grav)
    print "ymom added from gravitational sources = " + str(ymom_added_grav)
    print "zmom added from gravitational sources = " + str(zmom_added_grav)
    print "xmom added from rotation sources = " + str(xmom_added_rot)
    print "ymom added from rotation sources = " + str(ymom_added_rot)
    print "zmom added from rotation sources = " + str(zmom_added_rot)

    print ""
    print "Final diagnostics:"
    print ""

    print "Mass added = " + str(mass_added_neg_reset)
    print "Energy added = " + str(E_added_grav + E_added_flux + E_added_rot + E_added_reset + E_added_neg_reset)
    print "xmom added = " + str(xmom_added_flux + xmom_added_grav + xmom_added_rot)
    print "ymom added = " + str(ymom_added_flux + ymom_added_grav + ymom_added_rot)
    print "zmom added = " + str(zmom_added_flux + zmom_added_grav + zmom_added_rot)

    print ""

#    return [avg_timestep, avg_timestep_no_grav]



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

    return plotfiles



#
# Get the data column from a wdmerger diagnostic output file whose header name is col_name.
#

def get_column(col_name, diag_filename):

    # Open up the file for reading. Get the names of the columns, as well as a 2D list
    # with all the data.

    diag_file = open(diag_filename,'r')

    vc_line = 'git'

    # Skip the first few lines, they store the version control information

    line = diag_file.readline()

    while (line.split()[2] == vc_line):
        line = diag_file.readline()

    # The very next line will be the column headers

    col_names = line.split('  ')

    # Now read in the data

    diag_list = diag_file.readlines()
    data = []
    for line in diag_list:
        data.append( map(float, line[0:-1].split()) )

    diag_file.close()

    # Convert the data list into a 2D numpy array.

    data = np.array(data)

    # Let's do some cleanup.

    col_names.pop(0)                                        # Get rid of the # at the beginning
    col_names = [string.strip() for string in col_names]    # Remove any leading or trailing whitespace
    col_names = filter(None, col_names)                     # Remove any remaining blank entries

    # Obtain the column index and then return the column with that index.

    col_index = col_names.index(col_name)

    return data[:,col_index]



#
# Given a plotfile, return the location of the primary and the secondary.
#

def get_star_locs(plotfile):

    ds = yt.load(plotfile)

    # Get a numpy array corresponding to the density.

    problo = ds.domain_left_edge.v
    probhi = ds.domain_right_edge.v
    dim    = ds.domain_dimensions

    dx = (probhi - problo) / dim

    dens = (ds.covering_grid(level=0, left_edge=[0.0, 0.0, 0.0], dims=ds.domain_dimensions)['density']).v

    # Calculate the orbital parameters

    M_solar = 1.99e33
    Gconst = 6.67e-8

    M_P = 0.90
    M_S = 0.60

    M_P = M_P * M_solar
    M_S = M_S * M_solar

    # Get a numpy array corresponding to the density.

    a = (Gconst * (M_P + M_S) * rot_period**2 / (4.0 * np.pi**2))**(1.0/3.0)

    a_2 = a / (1 + M_S / M_P)
    a_1 = (M_S / M_P) * a_2

    # Guess the locations of the stars based on perfect circular rotation

    f = open(plotfile + '/job_info', 'r')

    for line in f:
        if string.find(line, "rotational_period") > 0:
            rot_period = float(string.split(line, "= ")[1])
            break

    f.close()

    t = (ds.current_time).v

    center = (probhi + problo) / 2.0

    loc_P = [-a_1 * np.cos(2 * np.pi * t / rot_period) + center[0], -a_1 * np.sin(2 * np.pi * t / rot_period) + center[1], 0.0 + center[2]]
    loc_S = [ a_2 * np.cos(2 * np.pi * t / rot_period) + center[0],  a_2 * np.sin(2 * np.pi * t / rot_period) + center[1], 0.0 + center[2]]

    loc_P = np.array(loc_P)
    loc_S = np.array(loc_S)

    # Create an array of the zone positions

    x = problo[0] + dx[0] * (np.arange(dim[0]) + 0.5e0)
    y = problo[1] + dx[1] * (np.arange(dim[1]) + 0.5e0)
    z = problo[2] + dx[2] * (np.arange(dim[2]) + 0.5e0)
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")

    rr = (xx**2 + yy**2 + zz**2)**0.5

    # Now what we'll do is to split up the grid into two parts.
    # zones that are closer to the primary's expected location and 
    # zones that are closer to the secondary's expected location.

    rr_P = ( (xx - loc_P[0])**2 + (yy - loc_P[1])**2 + (zz - loc_P[2])**2 )**0.5
    rr_S = ( (xx - loc_S[0])**2 + (yy - loc_S[1])**2 + (zz - loc_S[2])**2 )**0.5

    P_idx = np.where( rr_P < rr_S )
    S_idx = np.where( rr_S < rr_P )

    # Now, do a center of mass sum on each star.

    xx_P_com = np.sum( dens[P_idx] * xx[P_idx] ) / np.sum(dens[P_idx])
    yy_P_com = np.sum( dens[P_idx] * yy[P_idx] ) / np.sum(dens[P_idx])
    zz_P_com = np.sum( dens[P_idx] * zz[P_idx] ) / np.sum(dens[P_idx])

    xx_S_com = np.sum( dens[S_idx] * xx[S_idx] ) / np.sum(dens[S_idx])
    yy_S_com = np.sum( dens[S_idx] * yy[S_idx] ) / np.sum(dens[S_idx])
    zz_S_com = np.sum( dens[S_idx] * zz[S_idx] ) / np.sum(dens[S_idx])

    return [xx_P_com, yy_P_com, zz_P_com, xx_S_com, yy_S_com, zz_S_com]
