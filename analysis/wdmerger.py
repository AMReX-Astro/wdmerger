import os
import numpy as np
import string
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

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
    CASTRO_HOME = os.getenv('CASTRO_HOME')
    return CASTRO_HOME



# Return the last inputs file in a directory.

def get_inputs_filename(directory):

    # If there's a file named 'inputs', then this was a standard run.

    if (os.path.isfile(directory + '/inputs')):

        return 'inputs'

    # If there are files named inputs_*, corresponding to a chain, return the last one.

    elif (os.path.isfile(directory + '/inputs_1')):
        
        inputs_list = filter(lambda s: s[0:6] == "inputs", os.listdir(directory))

        # Now sort them numerically
        
        inputs_list.sort(key=lambda s: int(s[7:]))

        return inputs_list[-1]
 
       

# Get a CASTRO variable value from an inputs file.

def get_inputs_var(inputs, var):

    # Read in all the inputs file lines and search for the one
    # that starts with the desired variable name.

    inputs_file = open(inputs, 'r')
    lines = inputs_file.readlines()

    lines = filter(lambda s: s.split() != [], lines)
    line = filter(lambda s: s.split()[0] == var, lines)

    # The variable is the last item in a line before the comment.
    # This should work correctly even if there is no comment.

    var = (line[0].split('=')[1]).split('#')[0]
    var = var.strip()

    # Now, convert it into a list if it has multiple entries.

    if (var.split() != []):
        var = var.split()

    # Convert this to a floating point array, if possible.
    # If this fails, we'll just leave it as a string.

    try:
        var = np.array(var,dtype='float')    

        # Now convert this to an integer array, if possible.

        if (var[0].is_integer()):
            var = np.array(var,dtype='int')
    except:
        pass

    inputs_file.close()

    return var



# Get a variable value from a probin file.

def get_probin_var(probin, var):

    # Read in all the probin file lines and search for the one
    # that starts with the desired variable name.

    probin_file = open(probin, 'r')
    lines = probin_file.readlines()

    lines = filter(lambda s: s.split() != [], lines)
    line = filter(lambda s: s.split()[0] == var, lines)

    # The variable is the last item in a line before the comment.
    # This should work correctly even if there is no comment.

    var = (line[0].split('=')[1]).split('!')[0]
    var = var.strip()

    # Now, convert it into a list if it has multiple entries.

    if (var.split() != []):
        var = var.split()

    # Convert this to a floating point array, if possible.
    # If this fails, we'll just leave it as a string.

    try:
        var = np.array(var,dtype='float')    

        # Now convert this to an integer array, if possible.

        if (var[0].is_integer()):
            var = np.array(var,dtype='int')
    except:
        pass

    probin_file.close()

    return var



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
    microphysics_hash = ""
    for line in lines:
        if (len(line) == 4):
            if (line[0] == "Castro" and line[1] == "git" and line[2] == "hash:"):
                castro_hash = line[3]
            elif (line[0] == "BoxLib" and line[1] == "git" and line[2] == "hash:"):
                boxlib_hash = line[3]
            elif (line[0] == "wdmerger" and line[1] == "git" and line[2] == "hash:"):
                wdmerger_hash = line[3]
            elif (line[0] == "Microphysics" and line[1] == "git" and line[2] == "hash:"):
                microphysics_hash = line[3]

    job_info.close()

    return [castro_hash, boxlib_hash, wdmerger_hash, microphysics_hash]



#
# Given a diagnostic output file, return the CASTRO and BoxLib git commit hashes.
#

def get_git_commits_from_diagfile(diagfile):
    diagfile = open(diagfile, 'r')

    castro_hash   = ""
    boxlib_hash   = ""
    wdmerger_hash = ""
    microphysics_hash = ""

    line = diagfile.readline().split()

    if (line[1] == "Castro" and line[2] == "git" and line[3] == "hash:"):
        castro_hash = line[4]

    line = diagfile.readline().split()

    if (line[1] == "BoxLib" and line[2] == "git" and line[3] == "hash:"):
        boxlib_hash = line[4]

    line = diagfile.readline().split()

    if (line[1] == "Microphysics" and line[2] == "git" and line[3] == "hash:"):
        microphysics_hash = line[4]

    line = diagfile.readline().split()

    if (line[1] == "wdmerger" and line[2] == "git" and line[3] == "hash:"):
        wdmerger_hash = line[4]

    diagfile.close()

    return [castro_hash, boxlib_hash, wdmerger_hash, microphysics_hash]



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
    microphysics_hash = ""
    for line in lines:
        if (len(line) == 4):
            if (line[0] == "Castro" and line[1] == "git" and line[2] == "hash:"):
                castro_hash = line[3]
            elif (line[0] == "BoxLib" and line[1] == "git" and line[2] == "hash:"):
                boxlib_hash = line[3]
            elif (line[0] == "wdmerger" and line[1] == "git" and line[2] == "hash:"):
                wdmerger_hash = line[3]
            elif (line[0] == "Microphysics" and line[1] == "git" and line[2] == "hash:"):
                microphysics_hash = line[4]

    infofile.close()

    return [castro_hash, boxlib_hash, wdmerger_hash, microphysics_hash]



#
# Given CASTRO and BoxLib hashes that were used to create the plot for a given plotfile,
# insert these and the current Microphysics and wdmerger hashes into an EPS file.
# Credit: http://stackoverflow.com/questions/1325905/inserting-line-at-specified-position-of-a-text-file-in-python
# A comma in a print statement effectively prevents a newline from being appended to the print. This is valid in Python 2.x,
# but is not valid in Python 3.x. Instead, one should use the print() function in that case.
# Source: http://stackoverflow.com/questions/11266068/python-avoid-new-line-with-print-command
#

def insert_commits_into_eps(eps_file, data_file, data_file_type):
    """Insert git commit hashes into an EPS file."""

    import fileinput

    if (data_file_type == 'plot'):
        [castro_hash, boxlib_hash, wdmerger_hash, microphysics_hash] = get_git_commits_from_plotfile(data_file)
    elif (data_file_type == 'diag'):
        [castro_hash, boxlib_hash, wdmerger_hash, microphysics_hash] = get_git_commits_from_diagfile(data_file)
    elif (data_file_type == 'info'):
        [castro_hash, boxlib_hash, wdmerger_hash, microphysics_hash] = get_git_commits_from_infofile(data_file)
    else:
        print "Error: Data file type not recognized."

    input = fileinput.input(eps_file, inplace=True)

    for line in input:
        print line,
        if line.startswith('%%CreationDate:'):
            print "%%CASTRO git hash: " + castro_hash + "\n" + \
                  "%%BoxLib git hash: " + boxlib_hash + "\n" + \
                  "%%Microphysics git hash: " + microphysics_hash + "\n" + \
                  "%%wdmerger git hash: " + wdmerger_hash



#
# Given CASTRO and BoxLib hashes that were used to create the plot for a given plotfile,
# insert these and the current Microphysics and wdmerger hashes into a text file.
# We will append these to the end of the file, so that the code calling this should wait
# until just before it wants the commit hashes to appear to call it, and should probably
# not have the file actively open at the time.
#
# The default comment character will be chosen for use in LaTeX but can be changed if
# another comment character is needed.
#

def insert_commits_into_txt(txt_file, data_file, data_file_type, comment_char='%'):
    """Insert git commit hashes into a text file."""

    if (data_file_type == 'plot'):
        [castro_hash, boxlib_hash, wdmerger_hash, microphysics_hash] = get_git_commits_from_plotfile(data_file)
    elif (data_file_type == 'diag'):
        [castro_hash, boxlib_hash, wdmerger_hash, microphysics_hash] = get_git_commits_from_diagfile(data_file)
    elif (data_file_type == 'info'):
        [castro_hash, boxlib_hash, wdmerger_hash, microphysics_hash] = get_git_commits_from_infofile(data_file)
    else:
        print "Error: Data file type not recognized."

    string = comment_char + " CASTRO git hash: " + castro_hash + "\n" + \
             comment_char + " BoxLib git hash: " + boxlib_hash + "\n" + \
             comment_char + " Microphysics git hash: " + microphysics_hash + "\n" + \
             comment_char + " wdmerger git hash: " + wdmerger_hash + "\n"

    in_file = open(txt_file, 'a')
    in_file.write(string)
    in_file.close()



#
# Return the name of the latest wdmerger output file in a directory.
#

def get_last_output(directory):

    # Open up the standard output for analysis. It will be the numerically last file
    # starting with the designated output string.

    files = os.listdir(directory)

    files = sorted(filter(lambda s: s[0:9] == "wdmerger.",files))

    if (len(files) == 0):
        print "Error: No wdmerger output files in directory " + directory
        exit()

    return directory + "/" + files[len(files)-1]



#
# Return the name of the latest checkpoint in a directory.
#

def get_last_checkpoint(directory):

    if (not os.path.isdir(directory)):
        print "Error: directory " + directory + " does not exist in get_last_checkpoint()."
        raise

    # Doing a search this way will treat first any checkpoint files 
    # with seven digits, and then will fall back to ones with six and then five digits.
    # We want to be smart about this and list the ones in the current directory first,
    # before checking any output directories where the data is archived, because
    # the former are the most likely to be recently created checkpoints.

    checkpointList=[]
    checkpointNums=[]

    dirList = os.listdir(directory)

    from glob import glob

    def add_to_list(chkList, chkNums, chk_string):
        tempList = glob(chk_string)
        chkList += [temp_dir for temp_dir in tempList]
        chkNums += [temp_dir.split('/')[-1] for temp_dir in tempList]

    add_to_list(checkpointList, checkpointNums, directory + '/*chk???????')
    add_to_list(checkpointList, checkpointNums, directory + '/*/*chk???????')
    add_to_list(checkpointList, checkpointNums, directory + '/*chk??????')
    add_to_list(checkpointList, checkpointNums, directory + '/*/*chk??????')
    add_to_list(checkpointList, checkpointNums, directory + '/*chk?????')
    add_to_list(checkpointList, checkpointNums, directory + '/*/*chk?????')

    if not checkpointList or not checkpointNums:
        print "Error: no checkpoints found in directory " + directory
        raise

    # Match up the last checkpoint number with the actual file path location. 

    checkpoint = ""

    for chkNum in checkpointNums:

        for chkFile in checkpointList:

            currBaseName = chkFile.split('/')[-1]

            if currBaseName == chkNum:

                # The Header is the last thing written -- check if it's there, otherwise,
                # we can skip this iteration, because it means the latest checkpoint file 
                # is still being written.

                if os.path.isfile(chkFile + '/Header'):
                    checkpoint = chkFile
                    break

        if checkpoint:
            break

    if not checkpoint:
        print "Error: no completed checkpoint found in directory " + directory
        raise

    # Extract out the search directory from the result.

    checkpoint = checkpoint.split('/')[-1]

    return checkpoint



#
# Get the timing data for a wdmerger output file.
# Return the median timestep wall clock length.
#

def timing(output_filename):

    # Read in the file and return the following data in arrays:
    # - The coarse timestep time for each step
    # - The plotfile time for each step
    # For the sake of generality, we will count up how many
    # steps we have taken after the fact.

    output = open(output_filename, 'r')
    lines = output.readlines()
    coarseSteps = filter(lambda s: "Coarse" in s, lines)
    coarseSteps = [float(s.split('Coarse TimeStep time:')[1]) for s in coarseSteps] # Extract out the time only

    med_timestep = np.median(coarseSteps)

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

    if (len(grav_time) > 0):
        for n in range(len(coarseSteps)):
            for i in range(grav_per_timestep):
                coarseSteps[n] -= float(grav_time[n+i])

    med_timestep_no_grav = np.median(coarseSteps)

    output.close()

    return [med_timestep, med_timestep_no_grav]



# Extract timing data from a series of output files 
# and save it in a file using the BoxLib convention.
# See, e.g., MAESTRO/Docs/managing_jobs/scaling.

def boxlib_timing(output_filenames, data_filename):

    data = open(data_filename, 'w')

    data.write("    PROC     AVG     MIN     MAX\n")

    for output_filename in output_filenames:

        # Read in the file and return the following data in arrays:
        # - The coarse timestep time for each step
        # - The plotfile time for each step
        # For the sake of generality, we will count up how many
        # steps we have taken after the fact.

        output = open(output_filename, 'r')
        lines = output.readlines()
        coarseSteps = filter(lambda s: "Coarse" in s, lines)
        coarseSteps = [float(s.split('Coarse TimeStep time:')[1]) for s in coarseSteps] # Extract out the time only

        max_timestep = np.max(coarseSteps)
        avg_timestep = np.average(coarseSteps)
        min_timestep = np.min(coarseSteps)

        mpi = filter(lambda s: "MPI initialized" in s, lines)
        mpi = int(mpi[0].split()[3])

        omp = filter(lambda s: "OMP initialized" in s, lines)
        omp = int(omp[0].split()[3])

        nprocs = mpi * omp

        data.write("{:8d}".format(nprocs))
        data.write("{:8.1f}".format(avg_timestep))
        data.write("{:8.1f}".format(min_timestep))
        data.write("{:8.1f}".format(max_timestep))
        data.write("\n")

    data.write("\n")
    data.close()



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

    output.close()

#    return [avg_timestep, avg_timestep_no_grav]



#
# Get a sorted list of all the plotfiles in a directory.
#

def get_plotfiles(directory, prefix='plt'):

    # Check to make sure the directory exists.

    if (not os.path.isdir(directory)):
        print "Error: Directory " + directory + " does not exist, exiting."
        exit()

    # List all of the files in the directory.

    dir_contents = os.listdir(directory)

    # Strip out non-plotfiles.

    namlen = len(prefix)

    plotfiles = sorted(filter(lambda s: s[0:namlen] == prefix, dir_contents))

    # Remove any redundant plotfiles. These are generated with names like "old" or "temp"
    # and an easy way to get rid of them is just to check on all plotfiles that don't have
    # the right number of characters.

    plotfiles = filter(lambda s: len(s) == namlen + 5 or len(s) == namlen + 6, plotfiles)

    # This step strips out 'plt' from each plotfile, sorts them numerically,
    # then recasts them as strings. This is necessary because we may have some
    # plotfiles with six digits and some with five, and the default sort will
    # not sort them properly numerically.

    plotfiles = [prefix + str(y).zfill(5) for y in sorted([int(x[namlen:]) for x in plotfiles])]

    return plotfiles



#
# Given a plotfile, return the location of the primary and the secondary.
#

def get_star_locs(plotfile):

    import yt

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



# Get a variable from the CASTRO constants file.

def get_castro_const(var_name):

    CASTRO_HOME = get_castro_dir()

    file  = open(CASTRO_HOME + '/constants/constants_cgs.f90')
    lines = file.readlines()

    const = None

    for line in lines:
        lsplit = line.split()
        if (len(lsplit) >= 4):
            if (lsplit[3] == var_name):
                const = lsplit[5]

                # Remove the kind specification at the end of the string, if it exists.
                const = const.split('_')[0]

                # Convert from e notation to d notation.
                const = float(const.replace("d","e"))

    if (const == None):
        print "Constant " + var_name + " not found in CGS constants file; exiting."
        exit()

    file.close()

    return const



# Given a plotfile directory, return the simulation time.

def get_time_from_plotfile(pltfile):

    f = open(pltfile + '/Header', 'r')
    
    # Skip the first line, then read in the 
    # number of variables in the plotfile.

    line = f.readline()
    nspec = int(f.readline())

    # Now skip all the species, plus one more line.
    
    for n in range(nspec+2):
        line = f.readline()

    time = float(line[:-1])

    return time



# Given a run directory, determine whether the simulation has completed.

def is_dir_done(directory):

    if not os.path.isdir(directory):
        print "Error: directory " + directory + " does not exist in is_dir_done()."
        raise

    # If the directory already exists, check to see if we've reached the desired stopping point.
    # There are two places we can look: in the last checkpoint file, or in the last stdout file. 

    checkpoint = ""
    try:
        checkpoint = get_last_checkpoint(directory)
    except:
        pass

    last_output = ""
    try:
        last_output = get_last_output(directory)
    except:
        pass

    if not checkpoint and not last_output:
        print "Error: no checkpoint or stdout in directory " + directory
        done_status = 0
        return done_status

    # Get the desired stopping time and max step from the inputs file in the directory.
    # Alternatively, we may have set this from the calling script, so prefer that.

    inputs_filename = get_inputs_filename(directory)

    stop_time = -1.0
    max_step = -1

    if os.path.isfile(inputs_filename):
        stop_time = get_inputs_var(inputs_filename, "stop_time")
        max_step = get_inputs_var(inputs_filename, "max_step")

    # Assume we're not done, by default.

    time_flag = 0
    step_flag = 0

    done_status = 0

    if os.path.isfile(directory + '/' + checkpoint + '/jobIsDone'):

        # The problem has explicitly signalled that the simulation is complete; we can stop here.

        done_status = 1

    elif os.path.isfile(directory + '/' + checkpoint + '/jobIsNotDone'):

        # The problem has explicitly signalled that the simulation is NOT complete; again, we can stop here.

        done_status = 0

    elif os.path.isfile(directory + '/' + checkpoint + '/Header'):

        # Extract the checkpoint time. It is stored in row 3 of the Header file.

        chk_file = open(directory + '/' + checkpoint + '/Header')

        for i in range(3):
            line = chk_file.readline()

        # Convert to floating point, since it might be in exponential format.

        chk_time = float(line)

        # Extract the current timestep number.

        chk_step = checkpoint.split('chk')[-1]

        if stop_time >= 0.0:
            time_flag = chk_time >= stop_time

        if max_step >= 0:
            step_flag = chk_step >= max_step

        chk_file.close()

    elif last_output and os.path.isfile(directory + '/' + last_output):

        # Get the details of the last finished timestep.

        for line in reversed(open(directory + '/' + last_output)):
            if "STEP =" in line:
                output_time = float(line.split(' ')[5])
                output_step = int(line.split(' ')[2])
                break

        time_flag = output_time >= stop_time
        step_flag = output_step >= max_step



    if time_flag and step_flag:

        done_status = 1

    return done_status



def get_parameter_list(directory):
    """Return the list of parameters used in a directory.

    Skip all directories that do not contain a completed simulation.
    """

    param_list = []
    dirList = os.listdir(directory)

    for sub_directory in dirList:
        if is_dir_done(directory + '/' + sub_directory):
            param_list.append(sub_directory)

    return param_list



# Make a scatter-plot in the rho-T plane of all of the zones 
# in a plotfile, optionally coloring the points by the dominant
# species in that zone. We'll use a yt covering grid defined at 
# the coarse level for simplicity of access to the data.
# The inspiration for this is the upper-left panel of Figure 3 
# in Hawley et al., 2012.

def rho_T_scatterplot(output_filename, pltfile):

    import yt

    ds = yt.load(pltfile)

    grid = ds.covering_grid(level=0, left_edge = ds.domain_left_edge, dims = ds.domain_dimensions)

    # NumPy arrays of density and temperature.

    dens = grid['density'].v
    temp = grid['Temp'].v

    # Set the data limits. The upper limit in density is 10**8 g cm**-3,
    # which is a good upper limit for white dwarfs at our masses of interest. 
    # The lower limit is 10**4 g cm**-3, since we're not interested in the 
    # low-density material. In temperature, the upper limit will be 1e10 K 
    # because we won't burn hotter than that. For the lower limits, we will 
    # take 10**7 K because we aren't interested in material colder than that.

    min_dens = 1.0e4
    max_dens = 1.0e8

    min_temp = 1.0e7
    max_temp = 1.0e10

    # Create a mask that filters out all data we're not interested in.
    # This is done now because it has the benefit of speeding up the
    # plot creation because there's less data to work with.

    rho_T_mask = np.where(np.logical_and(temp > min_temp, dens > min_dens))

    dens = dens[rho_T_mask]
    temp = temp[rho_T_mask]

    # Get 12C, 28Si, and 56Ni masses.

    rho_c12  = grid['rho_c12'].v
    rho_si28 = grid['rho_si28'].v
    rho_ni56 = grid['rho_ni56'].v

    # Convert to mass fractions.

    xc12  = rho_c12[rho_T_mask]  / dens
    xsi28 = rho_si28[rho_T_mask] / dens
    xni56 = rho_ni56[rho_T_mask] / dens    

    # What we want is to make a mask array that determines which 
    # zones are dominated by each of these three species.

    c12_mask  = np.where(np.logical_and(xc12  > xsi28, xc12  > xni56))
    si28_mask = np.where(np.logical_and(xsi28 > xc12,  xsi28 > xni56))
    ni56_mask = np.where(np.logical_and(xni56 > xc12,  xni56 > xsi28))

    # Now create separate scatterplots for each of these three species.

    markersize = 40

    # The \! fixes an issue with spacing after a superscript in 
    # the version of matplotlib that is shipped with yt, though
    # it looks like it was fixed as of November 2015:
    # https://github.com/matplotlib/matplotlib/pull/4873

    plt.rcParams['mathtext.default'] = 'regular' # So the exponent is the same font as the text

    plt.scatter(dens[c12_mask],  temp[c12_mask],  color='g', s=markersize, marker='o', label=r'${}^{12\!}$C')
    plt.scatter(dens[si28_mask], temp[si28_mask], color='b', s=markersize, marker='s', label=r'${}^{28\!}$Si')
    plt.scatter(dens[ni56_mask], temp[ni56_mask], color='r', s=markersize, marker='d', label=r'${}^{56\!}$Ni')

    # Insert a buffer at the bottom of the plot since there will 
    # likely be a lot of points at the floor.

    min_temp = 0.7 * min_temp
    min_dens = 0.7 * min_dens

    plt.xlim([min_dens, max_dens])
    plt.ylim([min_temp, max_temp])

    plt.xscale('log')
    plt.yscale('log')

    # Axis labels and legend.

    plt.xlabel(r'Density (g / cm$^{-3\!}$)', fontsize=20)
    plt.ylabel(r'Temperature (K)', fontsize=20)
    plt.tick_params(labelsize=16)

    plt.legend(loc='upper left', fontsize=16, scatterpoints=1)

    # Save the plotfile.

    plt.tight_layout()
    plt.savefig(output_filename)
    insert_commits_into_eps(output_filename, pltfile, 'plot')

    plt.close()



def rho_T_sliceplot(output_filename, pltfile,
                    negate_left = False, scale_exp = 9, npix = 2000,
                    dens_range = [1.e-4, 1.e8], temp_range = [1.e7, 1.e10],
                    n_dens_ticks = 4, n_temp_ticks = 4):
    """Create a slice plot of rho and T using yt.

    negate_left: use negative tick marks for the left panel.
    scale_exp: the exponent (in base 10) used for scaling the axes.
    npix: number of pixels to use in the buffer for generating the image.
    dens_range: range of densities to plot.
    temp_range: range of temperatures to plot.
    n_dens_ticks: number of tick marks on the density colorbar.
    n_temp_ticks: number of tick marks on the temperature colorbar.
    """

    import yt
    import matplotlib
    from yt.visualization.api import get_multi_plot
    from matplotlib.colors import LogNorm

    fig, axes, colorbars = get_multi_plot(2, 1, colorbar = 'horizontal', bw = 6)

    ds = yt.load(pltfile)

    dim = ds.dimensionality

    buff = [npix, npix]

    bounds = [ds.domain_left_edge[0], ds.domain_right_edge[0], ds.domain_left_edge[1], ds.domain_right_edge[1]]

    if dim == 1:
        print "This slice plot routine is not implemented in one dimension."
        exit
    elif dim == 2:
        slc = ds.slice(2, 0.0)
    elif dim == 3:
        slc = ds.slice(2, 0.0)

    frb = yt.visualization.fixed_resolution.FixedResolutionBuffer(slc, bounds=bounds, buff_size=buff)

    dens_axis = axes[0][0]
    temp_axis = axes[0][1]

    if negate_left:
        dens_axis.set_xticks([-4.0, -2.0])
    else:
        dens_axis.set_xticks([4.0, 2.0])
    dens_axis.set_yticks([-4.0, -2.0, 0.0, 2.0, 4.0])

    temp_axis.set_xticks([0.0, 2.0, 4.0])
    temp_axis.set_yticks([-4.0, -2.0, 0.0, 2.0, 4.0])

    dens_axis.yaxis.tick_left()
    dens_axis.yaxis.set_label_position("left")
    temp_axis.yaxis.tick_right()
    temp_axis.yaxis.set_label_position("right")

    dens = np.array(frb['density'])
    temp = np.array(frb['Temp'])

    plots = []

    aspect = 1.0

    scale = 10**scale_exp

    if negate_left:
        left_bound = [ bounds[0] / scale, -bounds[1] / scale, bounds[2] / scale, bounds[3] / scale]
    else:
        left_bound = [ bounds[0] / scale, bounds[1] / scale, bounds[2] / scale, bounds[3] / scale]
    plots.append(dens_axis.imshow(dens, norm=LogNorm(), extent=left_bound, aspect=aspect))
    plots[-1].set_clim(dens_range[0], dens_range[1])
    plots[-1].set_cmap('bone')

    right_bound = [ bounds[0] / scale, bounds[1] / scale, bounds[2] / scale, bounds[3] / scale ]
    plots.append(temp_axis.imshow(temp, norm=LogNorm(), extent=right_bound, aspect=aspect))
    plots[-1].set_clim(temp_range[0], temp_range[1])
    plots[-1].set_cmap("hot")

    dens_ticks = np.logspace(np.log10(dens_range[0]), np.log10(dens_range[1]), num=n_dens_ticks)

    cb_dens = fig.colorbar(plots[0], cax=colorbars[0], ax=dens_axis, orientation='horizontal', ticks=dens_ticks)
    cb_dens.solids.set_rasterized(True)
    cb_dens.set_label(r'$\mathrm{Density}\ (\mathrm{g\ cm^{-3}})$')

    temp_ticks = np.logspace(np.log10(temp_range[0]), np.log10(temp_range[1]), num=n_temp_ticks)

    cb_temp = fig.colorbar(plots[1], cax=colorbars[1], ax=temp_axis, orientation='horizontal', ticks=temp_ticks)
    cb_temp.solids.set_rasterized(True)
    cb_temp.set_label(r'$\mathrm{Temperature}\ (\mathrm{K})$')

    dens_axis.set_xlim(dens_axis.get_xlim()[::-1])

    dens_axis.set_position([0.125+0.0575, 0.075, 0.375, 0.75])
    temp_axis.set_position([0.500-0.0575, 0.075, 0.375, 0.75])

    colorbars[0].set_position([0.275, 0.92, 0.2, 0.075])
    colorbars[1].set_position([0.525, 0.92, 0.2, 0.075])

    dens_axis.set_xlabel(r'$r\ (10^{}\ \mathrm{{cm}})$'.format('{' + str(scale_exp) + '}'))
    temp_axis.set_xlabel(r'$r\ (10^{}\ \mathrm{{cm}})$'.format('{' + str(scale_exp) + '}'))

    dens_axis.set_ylabel(r'$z\ (10^9\ \mathrm{cm})$')
    temp_axis.set_ylabel(r'$z\ (10^9\ \mathrm{cm})$')

    matplotlib.rcParams.update({'font.size': 20})

    # Save as EPS

    plt.savefig(output_filename)

    insert_commits_into_eps(output_filename, pltfile, 'plot')

    # Save as PNG

    file_base = output_filename[:-4]

    fig.savefig(file_base + '.png', bbox_inches='tight')

    # Convert to JPG

    try:
        os.system('convert ' + file_base + '.png ' + file_base + '.jpg')
    except:
        pass

    plt.close()



# A routine for doing axis-aligned slice plots over a given field.

def slice_plot(field, output_filename, pltfile, idir=3):
    """Create an axis-aligned slice plot over a given field with yt."""

    import yt

    ds = yt.load(pltfile)

    dim = ds.dimensionality

    fields = [field]

    if dim == 1:

        print "This slice plot routine is not implemented in one dimension."
        exit

    elif dim == 2:

        if idir == 1:
            axis = 'r'
        elif idir == 2:
            axis = 'z'
        elif idir == 3:
            axis = 'theta'
        else:
            print "Unknown direction for slicing in slice_plot."
            exit


    elif dim == 3:

        if idir == 1:
            axis = 'x'
        elif idir == 2:
            axis = 'y'
        elif idir == 3:
            axis = 'z'
        else:
            print "Unknown direction for slicing in slice_plot."


    sp = yt.SlicePlot(ds, 'z', fields=fields)

    sp.set_cmap(field, 'hot')

    plot = sp.plots[field]

    cb = plot.cb
    cb.solids.set_rasterized(True)

    sp._setup_plots()

    plt.savefig(output_filename)

    insert_commits_into_eps(output_filename, pltfile, 'plot')

    plt.savefig(output_filename[:-4] + '.png')

    plt.close()



# Make a movie from the JPG files in a directory.

def make_movie(output_dir, mpg_filename):
    """Make an MPEG movie from the JPG files in a directory."""

    import glob

    curr_dir = os.getcwd()

    os.chdir(output_dir)

    # Assume that we can sort them chronologically with a standard
    # alphanumeric sort.

    jpg_list = sorted(glob.glob('*.jpg'))

    # Exclude any existing symbolic links.

    jpg_list = filter(lambda jpg: 'link' not in jpg, jpg_list)

    # Make a list of files indexed by monotonically increasing integers.

    idx = 0

    for jpg in jpg_list:

        idx_str = '{:05d}'.format(idx)
        link_name = 'link_' + idx_str + '.jpg'

        try:
            os.system('rm -f ' + link_name)
            os.system('ln -s ' + jpg + ' ' + link_name)
        except:
            pass

        idx += 1

    os.chdir(curr_dir)

    try:
        os.system('ffmpeg -i ' + output_dir + 'link_\%05d.jpg ' + mpg_filename)
    except:
        print "Error: could not successfully make a movie with ffmpeg."
        pass

    os.system('rm -f link_*.jpg')
