
#
# Return the CASTRO directory on the current machine, based on your environment variables.
#

def get_castro_dir():
    """Return the location of the CASTRO directory."""

    import os

    CASTRO_HOME = os.getenv('CASTRO_HOME')

    return CASTRO_HOME



# Return the name of the current inputs file in a directory.

def get_inputs_filename(directory):
    """Return the name of the inputs file in a directory."""

    import os

    # At present we have no reason to look for anything other than inputs.

    if os.path.isfile(directory + '/inputs'):

        return 'inputs'

    elif os.path.isfile(directory + '/inputs_2d'):

        return 'inputs_2d'

    elif os.path.isfile(directory + '/inputs_3d'):

        return 'inputs_3d'

    else:

        print("Error: no inputs file found in " + directory + ".")
        exit



# Get a CASTRO variable value from an inputs file.

def get_inputs_var(inputs, var):
    """Retrieve a CASTRO variable value from an inputs file."""

    import numpy as np

    # Read in all the inputs file lines and search for the one
    # that starts with the desired variable name.

    inputs_file = open(inputs, 'r')
    lines = inputs_file.readlines()

    lines = list(filter(lambda s: s.split() != [], lines))
    line = list(filter(lambda s: s.split()[0] == var, lines))

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
    """Retrieve a variable value from a probin file."""

    import numpy as np

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
# Given a plotfile directory, return the git commit hashes.
#

def get_git_commits_from_plotfile(plotfile):
    """Retrieve git commit hashes from a plotfile."""

    job_info = open(plotfile + "/job_info", 'r')
    lines = job_info.readlines()
    lines = [line.split() for line in lines]

    castro_hash   = ""
    amrex_hash   = ""
    microphysics_hash = ""

    for line in lines:
        if (len(line) == 4):
            if (line[0] == "Castro" and line[1] == "git" and line[2] == "hash:"):
                castro_hash = line[3]
            elif (line[0] == "AMReX" and line[1] == "git" and line[2] == "hash:"):
                amrex_hash = line[3]
            elif (line[0] == "Microphysics" and line[1] == "git" and line[2] == "hash:"):
                microphysics_hash = line[3]

    job_info.close()

    return [castro_hash, amrex_hash, microphysics_hash]



#
# Given a diagnostic output file, return the git commit hashes.
#

def get_git_commits_from_diagfile(diagfile):
    """Retrieve git commit hashes from a diagnostic file."""

    diagfile = open(diagfile, 'r')

    castro_hash   = ""
    amrex_hash   = ""
    microphysics_hash = ""

    line = diagfile.readline().split()

    if (line[1] == "Castro" and line[2] == "git" and line[3] == "hash:"):
        castro_hash = line[4]

    line = diagfile.readline().split()

    if (line[1] == "AMReX" and line[2] == "git" and line[3] == "hash:"):
        amrex_hash = line[4]

    line = diagfile.readline().split()

    if (line[1] == "Microphysics" and line[2] == "git" and line[3] == "hash:"):
        microphysics_hash = line[4]

    diagfile.close()

    return [castro_hash, amrex_hash, microphysics_hash]



#
# Given the stdout from a Castro run, return the git commit hashes.
#

def get_git_commits_from_infofile(infofile):
    """Retrieve git commit hashes from a stdout file."""

    infofile = open(infofile, 'r')
    lines = infofile.readlines()
    lines = [line.split() for line in lines]

    castro_hash   = ""
    amrex_hash   = ""
    microphysics_hash = ""

    for line in lines:
        if (len(line) == 4):
            if (line[0] == "Castro" and line[1] == "git" and line[2] == "hash:"):
                castro_hash = line[3]
            elif (line[0] == "AMReX" and line[1] == "git" and line[2] == "hash:"):
                amrex_hash = line[3]
            elif (line[0] == "Microphysics" and line[1] == "git" and line[2] == "hash:"):
                microphysics_hash = line[4]

    infofile.close()

    return [castro_hash, amrex_hash, microphysics_hash]



#
# Given CASTRO and AMReX hashes that were used to create the plot for a given plotfile,
# insert these and the current Microphysics and wdmerger hashes into an EPS file.
# Credit: http://stackoverflow.com/questions/1325905/inserting-line-at-specified-position-of-a-text-file-in-python
#

def insert_commits_into_eps(eps_file, data_file, data_file_type):
    """Insert git commit hashes into an EPS file."""

    import fileinput

    if (data_file_type == 'plot'):
        [castro_hash, amrex_hash, microphysics_hash] = get_git_commits_from_plotfile(data_file)
    elif (data_file_type == 'diag'):
        [castro_hash, amrex_hash, microphysics_hash] = get_git_commits_from_diagfile(data_file)
    elif (data_file_type == 'info'):
        [castro_hash, amrex_hash, microphysics_hash] = get_git_commits_from_infofile(data_file)
    else:
        print("Error: Data file type not recognized.")

    input = fileinput.input(eps_file, inplace=True)

    for line in input:
        print(line, end="") # No additional newline
        if line.startswith('%%CreationDate:'):
            print("%%CASTRO git hash: " + castro_hash + "\n" + \
                  "%%AMReX git hash: " + amrex_hash + "\n" + \
                  "%%Microphysics git hash: " + microphysics_hash)



#
# Given CASTRO and AMReX hashes that were used to create the plot for a given plotfile,
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
        [castro_hash, amrex_hash, microphysics_hash] = get_git_commits_from_plotfile(data_file)
    elif (data_file_type == 'diag'):
        [castro_hash, amrex_hash, microphysics_hash] = get_git_commits_from_diagfile(data_file)
    elif (data_file_type == 'info'):
        [castro_hash, amrex_hash, microphysics_hash] = get_git_commits_from_infofile(data_file)
    else:
        print("Error: Data file type not recognized.")

    string = comment_char + " CASTRO git hash: " + castro_hash + "\n" + \
             comment_char + " AMReX git hash: " + amrex_hash + "\n" + \
             comment_char + " Microphysics git hash: " + microphysics_hash + "\n"

    in_file = open(txt_file, 'a')
    in_file.write(string)
    in_file.close()



#
# Return the name of the latest wdmerger output file in a directory.
#

def get_last_output(directory):
    """Obtain the name of the the last wdmerger output file in a directory."""

    import os

    # Open up the standard output for analysis. It will be the numerically last file
    # starting with the designated output string.

    files = os.listdir(directory + '/output')

    files = sorted(filter(lambda s: s[0:9] == "wdmerger.",files))

    if (len(files) == 0):
        print("Error: No wdmerger output files in directory " + directory)
        exit()

    return directory + "/" + files[len(files)-1]



#
# Return the name of the latest checkpoint in a directory.
#

def get_last_checkpoint(directory):
    """Obtain the name of the last checkpoint in a directory."""

    import os

    if (not os.path.isdir(directory)):
        print("Error: directory " + directory + " does not exist in get_last_checkpoint().")
        return

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
        print("Error: no checkpoints found in directory " + directory)
        return

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
        print("Error: no completed checkpoint found in directory " + directory)
        return

    # Extract out the search directory from the result.

    checkpoint = checkpoint.split('/')[-1]

    return checkpoint



#
# Get the data column from a wdmerger diagnostic output file whose header name is col_name.
#

def get_column(col_name, diag_filename):
    """Get a column of data from a diagnostic file."""

    import numpy as np

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
        data.append( list(map(float, line[0:-1].split())) )

    diag_file.close()

    # Convert the data list into a 2D numpy array.

    data = np.array(data)

    # Let's do some cleanup.

    col_names.pop(0)                                        # Get rid of the # at the beginning
    col_names = [string.strip() for string in col_names]    # Remove any leading or trailing whitespace
    col_names = list(filter(None, col_names))               # Remove any remaining blank entries

    # Obtain the column index and then return the column with that index.

    col_index = col_names.index(col_name)

    return data[:,col_index]



#
# Get the timing data for a wdmerger output file.
# Return the median timestep wall clock length.
#

def timing(output_filename):
    """Return the median wall time per simulation timestep."""

    import numpy as np

    output = open(output_filename, 'r')
    lines = output.readlines()

    coarseSteps = filter(lambda s: "Coarse" in s, lines)

    # Extract out the time only

    coarseSteps = [float(s.split('Coarse TimeStep time:')[1]) for s in coarseSteps]

    med_timestep = np.median(coarseSteps)

    output.close()

    return med_timestep



# Extract timing data from a series of output files 
# and save it in a file using the AMReX convention.
# See, e.g., MAESTRO/Docs/managing_jobs/scaling.

def amrex_timing(output_filenames, data_filename):
    """Extract timing data and save it to an output file."""

    import numpy as np

    data = open(data_filename, 'w')

    data.write("    PROC     AVG     MIN     MAX\n")

    for output_filename in output_filenames:

        output = open(output_filename, 'r')
        lines = output.readlines()

        coarseSteps = filter(lambda s: "Coarse" in s, lines)

        # Extract out the time only

        coarseSteps = [float(s.split('Coarse TimeStep time:')[1]) for s in coarseSteps]

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
    """Add up all the energy and momentum losses reported from castro.print_energy_diagnostics."""

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

    print("")
    print("Analysis of output file " + output_filename)
    print("")

    print("Mass added from negative density resets = " + str(mass_added_neg_reset))
    print("Energy added from negative density resets = " + str(E_added_neg_reset))
    print("Energy added from gravitational sources = " + str(E_added_grav))
    print("Energy added from rotation sources = " + str(E_added_rot))
    print("Energy added from hydro fluxes = " + str(E_added_flux))
    print("Energy added from resets = " + str(E_added_reset))
    print("xmom added from hydro fluxes = " + str(xmom_added_flux))
    print("ymom added from hydro fluxes = " + str(ymom_added_flux))
    print("zmom added from hydro fluxes = " + str(zmom_added_flux))
    print("xmom added from gravitational sources = " + str(xmom_added_grav))
    print("ymom added from gravitational sources = " + str(ymom_added_grav))
    print("zmom added from gravitational sources = " + str(zmom_added_grav))
    print("xmom added from rotation sources = " + str(xmom_added_rot))
    print("ymom added from rotation sources = " + str(ymom_added_rot))
    print("zmom added from rotation sources = " + str(zmom_added_rot))

    print("")
    print("Final diagnostics:")
    print("")

    print("Mass added = " + str(mass_added_neg_reset))
    print("Energy added = " + str(E_added_grav + E_added_flux + E_added_rot + E_added_reset + E_added_neg_reset))
    print("xmom added = " + str(xmom_added_flux + xmom_added_grav + xmom_added_rot))
    print("ymom added = " + str(ymom_added_flux + ymom_added_grav + ymom_added_rot))
    print("zmom added = " + str(zmom_added_flux + zmom_added_grav + zmom_added_rot))

    print("")

    output.close()

    return



#
# Get a sorted list of all the plotfiles in a directory.
#

def get_plotfiles(directory, prefix = 'plt'):
    """Get a sorted list of all the plotfiles in a directory."""

    import os

    # Check to make sure the directory exists.

    if not os.path.isdir(directory):
        print("Error: Directory " + directory + " does not exist, exiting.")
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
    """Given a plotfile, return the location of the primary and the secondary."""

    import numpy as np
    import yt
    import string

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
    """Get a variable from the CASTRO constants file."""

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
        print("Constant " + var_name + " not found in CGS constants file; exiting.")
        exit()

    file.close()

    return const



# Given a plotfile directory, return the simulation time.

def get_time_from_plotfile(pltfile):
    """Given a plotfile, return the simulation time."""

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



# Given a list of plotfiles, return the plotfile whose simulation
# time is closest to the desired time.

def get_nearest_plotfile(pltfiles, time):
    """Given a list of plotfiles, return the one whose time is nearest to a given target."""

    import numpy as np

    times = np.array([get_time_from_plotfile(pltfile) for pltfile in pltfiles])

    idx = (np.abs(times - time)).argmin()

    return pltfiles[idx]



# Given a run directory, determine whether the simulation has completed.

def is_dir_done(directory):
    """Given a run directory, determine whether the simulation has completed."""

    import os

    if not os.path.isdir(directory):
        print("Error: directory " + directory + " does not exist in is_dir_done().")
        return

    if os.path.isfile(directory + '/jobIsDone'):

        # The user has explicitly signalled that the simulation is complete; we can stop here.

        done_status = 1

    else:

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
            print("Error: no checkpoint or stdout in directory " + directory)
            done_status = 0
            return done_status

        # Get the desired stopping time and max step from the inputs file in the directory.
        # Alternatively, we may have set this from the calling script, so prefer that.

        inputs_filename = directory + '/' + get_inputs_filename(directory)

        stop_time = -1.0
        max_step = -1

        if os.path.isfile(inputs_filename):
            stop_time = get_inputs_var(inputs_filename, "stop_time")[0]
            max_step = get_inputs_var(inputs_filename, "max_step")[0]

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
                time_flag = float(chk_time) >= float(stop_time)

            if max_step >= 0:
                step_flag = int(chk_step) >= int(max_step)

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

        if time_flag or step_flag:

            done_status = 1

    if done_status == 1 and not os.path.isfile(directory + '/jobIsDone'):
        os.system('touch %s' % directory + '/jobIsDone')

    return done_status



def get_parameter_list(directory):
    """Return the list of parameters used in a directory.

    Skip all directories that do not contain a completed simulation.
    """

    import os

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
    """Plot density and temperature on a scatterplot."""

    import numpy as np
    import yt
    import matplotlib.pyplot as plt

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
                    x_ticks = [2.0e9, 4.0e9], y_ticks = [2.0e9, 4.0e9],
                    n_dens_ticks = 4, n_temp_ticks = 4,
                    domain_frac = 1.0):
    """Create a slice plot of rho and T using yt.

    negate_left: use negative tick marks for the left panel.
    scale_exp: the exponent (in base 10) used for scaling the axes.
    npix: number of pixels to use in the buffer for generating the image.
    dens_range: range of densities to plot.
    temp_range: range of temperatures to plot.
    n_dens_ticks: number of tick marks on the density colorbar.
    n_temp_ticks: number of tick marks on the temperature colorbar.
    """

    import numpy as np
    import yt
    import matplotlib
    import matplotlib.pyplot as plt
    from yt.visualization.api import get_multi_plot
    from matplotlib.colors import LogNorm

    fig, axes, colorbars = get_multi_plot(2, 1, colorbar = 'horizontal', bw = 6)

    ds = yt.load(pltfile)

    dim = ds.dimensionality

    buff = [npix, npix]

    bounds = [domain_frac * ds.domain_left_edge[0], domain_frac * ds.domain_right_edge[0], 
              domain_frac * ds.domain_left_edge[1], domain_frac * ds.domain_right_edge[1]]

    if dim == 1:
        print("This slice plot routine is not implemented in one dimension.")
        exit
    elif dim == 2:
        slc = ds.slice(2, 0.0)
    elif dim == 3:
        slc = ds.slice(2, 0.0)

    frb = yt.visualization.fixed_resolution.FixedResolutionBuffer(slc, bounds=bounds, buff_size=buff)

    dens_axis = axes[0][0]
    temp_axis = axes[0][1]

    scale = 10**scale_exp

    if negate_left:
        dens_axis.set_xticks([-x_ticks[1] / scale, -x_ticks[0] / scale])
    else:
        dens_axis.set_xticks([x_ticks[1] / scale, x_ticks[0] / scale])
    dens_axis.set_yticks([-y_ticks[1] / scale, -y_ticks[0] / scale, 
                           0.0, 
                           y_ticks[0] / scale,  y_ticks[1] / scale])

    temp_axis.set_xticks([0.0, x_ticks[0] / scale, x_ticks[1] / scale])
    temp_axis.set_yticks([-y_ticks[1] / scale, -y_ticks[0] / scale, 
                           0.0, 
                           y_ticks[0] / scale,  y_ticks[1] / scale])

    dens_axis.yaxis.tick_left()
    dens_axis.yaxis.set_label_position("left")
    temp_axis.yaxis.tick_right()
    temp_axis.yaxis.set_label_position("right")

    dens = np.array(frb['density'])
    temp = np.array(frb['Temp'])

    plots = []

    aspect = 1.0

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

    dens_axis.set_ylabel(r'$z\ (10^{}\ \mathrm{{cm}})$'.format('{' + str(scale_exp) + '}'))
    temp_axis.set_ylabel(r'$z\ (10^{}\ \mathrm{{cm}})$'.format('{' + str(scale_exp) + '}'))

    matplotlib.rcParams.update({'font.size': 20})

    # Save as EPS

    fig.savefig(output_filename, bbox_inches='tight')

    insert_commits_into_eps(output_filename, pltfile, 'plot')

    # Save as JPG

    fig.savefig(output_filename.replace('eps', 'jpg'), bbox_inches='tight')

    plt.close()



# A routine for doing axis-aligned slice plots over a given field.

def slice_plot(field, output_filename, pltfile, idir = 3):
    """Create an axis-aligned slice plot over a given field with yt."""

    import yt
    import matplotlib.pyplot as plt

    ds = yt.load(pltfile)

    dim = ds.dimensionality

    fields = [field]

    if dim == 1:

        print("This slice plot routine is not implemented in one dimension.")
        exit

    elif dim == 2:

        if idir == 1:
            axis = 'r'
        elif idir == 2:
            axis = 'z'
        elif idir == 3:
            axis = 'theta'
        else:
            print("Unknown direction for slicing in slice_plot.")
            exit


    elif dim == 3:

        if idir == 1:
            axis = 'x'
        elif idir == 2:
            axis = 'y'
        elif idir == 3:
            axis = 'z'
        else:
            print("Unknown direction for slicing in slice_plot.")


    sp = yt.SlicePlot(ds, axis, fields=fields)

    sp.set_cmap(field, 'hot')

    plot = sp.plots[field]

    cb = plot.cb
    cb.solids.set_rasterized(True)

    sp._setup_plots()

    plt.savefig(output_filename)

    insert_commits_into_eps(output_filename, pltfile, 'plot')

    plt.savefig(output_filename[:-4] + '.png')

    plt.close()



# A routine for doing axis-aligned multi-panel slice plots over a given field.

def multipanel_slice_plot(field, output_filename, pltfiles, idir = 3,
                          zlim = None, colormap = 'hot', scale_exp = 9,
                          nrows = 2, ncols = 2, axes_pad = 0.10,
                          zoom = 1.0,
                          xticks = None, yticks = None,
                          rect = [0.03,0.075,0.92,0.90],
                          annotate_time = False):
    """Create an axis-aligned multi-panel slice plot over a given field with yt."""

    import yt
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import AxesGrid

    # Sanity check: are there enough plots to fill the axes grid?

    if nrows * ncols != len(pltfiles):
        print("Error: not enough plots for the multipanel slice plot.")
        exit

    fig = plt.figure()

    # Set up the AxesGrid.

    grid = AxesGrid(fig,
                    rect,
                    nrows_ncols = (nrows, ncols),
                    axes_pad = axes_pad,
                    label_mode = "L",
                    share_all = False,
                    cbar_location = "right",
                    cbar_mode = "single",
                    cbar_size = "5%",
                    cbar_pad = "1%")

    for i, pltfile in enumerate(pltfiles):

        ds = yt.load(pltfile)

        dim = ds.dimensionality

        fields = [field]

        if dim == 1:

            print("This slice plot routine is not implemented in one dimension.")
            exit

        elif dim == 2:

            if idir == 1:
                axis = 'r'
            elif idir == 2:
                axis = 'z'
            elif idir == 3:
                axis = 'theta'
            else:
                print("Unknown direction for slicing in multipanel_slice_plot.")
                exit

        elif dim == 3:

            if idir == 1:
                axis = 'x'
            elif idir == 2:
                axis = 'y'
            elif idir == 3:
                axis = 'z'
            else:
                print("Unknown direction for slicing in multipanel_slice_plot.")
                exit

        else:

            print("Error: Invalid dataset dimensionality.")
            exit

        # Set up the SlicePlot

        sp = yt.SlicePlot(ds, axis, fields = fields, axes_unit = '1e' + str(scale_exp) + '*cm')
        sp.set_xlabel(r'$\rm{{x}}\ (10^{}\ \rm{{cm}})$'.format('{' + str(scale_exp) + '}'))
        sp.set_ylabel(r'$\rm{{y}}\ (10^{}\ \rm{{cm}})$'.format('{' + str(scale_exp) + '}'))
        sp.set_minorticks(field, 'off')
        sp.zoom(zoom)

        sp.set_cmap(field, colormap)
        sp.set_cbar_minorticks(field, 'off')

        if zlim is not None:
            sp.set_zlim(field, zlim[0], zlim[1])

        if annotate_time:
            sp.annotate_timestamp(corner='upper_left')

        plot = sp.plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]

        cb = plot.cb
        cb.solids.set_rasterized(True)

        sp._setup_plots()

        if xticks is not None:
            grid.axes_all[i].set_xticks([xtick / 10.0**scale_exp for xtick in xticks])

        if yticks is not None:
            grid.axes_all[i].set_yticks([ytick / 10.0**scale_exp for ytick in yticks])

    plt.savefig(output_filename)

    insert_commits_into_eps(output_filename, pltfile, 'plot')

    plt.close()



# Make a movie from a list of input image files in a directory.

def make_movie(output_dir, img_list, mpg_filename):
    """Make an MPEG movie from a list of image files."""

    import os
    import shutil
    import glob

    cwd = os.getcwd()

    # Get file extension; assume all files have the same format.

    file_format = '.' + img_list[0].split('.')[-1]

    # Copy them to a series of files that are indexed by
    # monotonically increasing integers (this is needed for ffmpeg).

    in_list = [output_dir + '/temp_' + '{:05d}'.format(i) + file_format for i, img in enumerate(img_list)]

    for img_new, img_old in zip(in_list, img_list):
        try:
            shutil.copyfile(img_old, img_new)
        except:
            print("Error: source file " + img_old + " does not exist.")
            pass

    try:
        os.system('ffmpeg -i ' + output_dir + '/temp_\%05d' + file_format + ' -b:v 20M ' + mpg_filename)
    except:
        print("Error: could not successfully make a movie with ffmpeg.")
        pass

    for img in in_list:
        os.remove(img)



def vol_render_density(outfile, ds):
    """Volume render the density given a yt dataset."""

    import numpy as np
    import yt
    import matplotlib
    matplotlib.use('agg')
    from yt.visualization.volume_rendering.api import Scene, VolumeSource
    import matplotlib.pyplot as plt

    ds.periodicity = (True, True, True)

    field = ('boxlib', 'density')
    ds._get_field_info(field).take_log = True

    sc = Scene()

    # Add a volume: select a sphere

    vol = VolumeSource(ds, field=field)
    vol.use_ghost_zones = True

    sc.add_source(vol)

    # Transfer function

    vals = [-1, 0, 1, 2, 3, 4, 5, 6, 7]
    sigma = 0.1

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()

    cm = "spectral"

    for v in vals:
        if v < 3:
            alpha = 0.1
        else:
            alpha = 0.5
        tf.sample_colormap(v, sigma**2, colormap=cm, alpha=alpha)

    sc.get_source(0).transfer_function = tf

    cam = sc.add_camera(ds, lens_type="perspective")
    cam.resolution = (1920, 1080)

    center = 0.5 * (ds.domain_left_edge + ds.domain_right_edge)
    width = ds.domain_width

    # Set the camera so that we're looking down on the xy plane from a 45
    # degree angle. We reverse the y-coordinate since yt seems to use the
    # opposite orientation convention to us (where the primary should be
    # on the left along the x-axis). We'll scale the camera position based
    # on a zoom factor proportional to the width of the domain.

    zoom_factor = 0.75

    cam_position = np.array([center[0], center[1] - zoom_factor * width[1], center[2] + zoom_factor * width[2]])

    cam.position = zoom_factor * ds.arr(cam_position, 'cm')

    # Set the normal vector so that we look toward the center.

    normal = (center - cam.position)
    normal /= np.sqrt(normal.dot(normal))

    cam.switch_orientation(normal_vector = normal, north_vector = [0.0, 0.0, 1.0])
    cam.set_width(width)

    # Render the image.

    sc.render()

    # Save the image without annotation.

    sc.save(outfile, sigma_clip=6.0)

    # Save the image with a colorbar.

    sc.save_annotated(outfile.replace(".png", "_colorbar.png"), sigma_clip=6.0)

    # Save the image with a colorbar and the current time.

    sc.save_annotated(outfile.replace(".png", "_colorbar_time.png"), sigma_clip=6.0,
                      text_annotate=[[(0.05, 0.925),
                                      "t = {:.2f} s".format(float(ds.current_time.d)),
                                      dict(horizontalalignment="left",
                                           fontsize="20")]])
