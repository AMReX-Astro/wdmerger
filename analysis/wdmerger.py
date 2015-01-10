import os
import numpy as np



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
# Given a plotfile directory, return the CASTRO and BoxLib git commit hashes.
#

def get_git_commits_from_plotfile(plotfile):
    job_info = open(plotfile + "/job_info", 'r')
    lines = job_info.readlines()
    lines = [line.split() for line in lines]
    for line in lines:
        if (len(line) == 4):
            if (line[0] == "Castro" and line[1] == "git" and line[2] == "hash:"):
                castro_hash = line[3]
            elif (line[0] == "BoxLib" and line[1] == "git" and line[2] == "hash:"):
                boxlib_hash = line[3]

    job_info.close()

    return [castro_hash, boxlib_hash]



#
# Given a diagnostic output file, return the CASTRO and BoxLib git commit hashes.
#

def get_git_commits_from_diagfile(diagfile):
    diagfile = open(diagfile, 'r')
    castro_hash = (diagfile.readline().split())[4]
    boxlib_hash = (diagfile.readline().split())[4]
    diagfile.close()

    return [castro_hash, boxlib_hash]



#
# Given the stdout from a Castro run, return the CASTRO and BoxLib git commit hashes.
#

def get_git_commits_from_infofile(infofile):
    infofile = open(infofile, 'r')
    lines = infofile.readlines()
    lines = [line.split() for line in lines]
    for line in lines:
        if (len(line) == 4):
            if (line[0] == "Castro" and line[1] == "git" and line[2] == "hash:"):
                castro_hash = line[3]
            elif (line[0] == "BoxLib" and line[1] == "git" and line[2] == "hash:"):
                boxlib_hash = line[3]

    infofile.close()

    return [castro_hash, boxlib_hash]



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
        [castro_hash, boxlib_hash] = get_git_commits_from_plotfile(data_file)
    elif (data_file_type == 'diag'):
        [castro_hash, boxlib_hash] = get_git_commits_from_diagfile(data_file)
    elif (data_file_type == 'info'):
        [castro_hash, boxlib_hash] = get_git_commits_from_infofile(data_file)
    else:
        print "Error: Data file type not recognized."

    wdmerger_hash = get_wdmerger_git_commit_hash();

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
    castro_hash_line = diag_file.readline()
    boxlib_hash_line = diag_file.readline()
    col_names = diag_file.readline().split('  ')
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
