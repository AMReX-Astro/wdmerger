#!/bin/bash 

source $WDMERGER_HOME/job_scripts/run_utils.sh

#----------------------------------------------------------------------------
# User modifiable variables:

# pidfile is a lock file that is used to make sure that only one instance 
# of this script is working on the current directory
pidfile=process.pid

# Set the prefix of the plotfiles and checkpoint files
plt_prefix=*plt
chk_prefix=*chk

# Directory to archive to on the storage system -- set this to the working directory.
# What this does is to extract the part of the directory name that is everything but
# the variable $WORKDIR, which is assumed to be the location you're running in.
# This will exactly reproduce the file structure on the storage system as it appears
# on the compute. It relies on the bash trick that # removes a prefix from a variable:
# http://mintaka.sdsu.edu/GF/bibliog/latex/debian/bash.html
cwd=$(pwd)
storage_dir=$(echo ${cwd#$workdir})

#----------------------------------------------------------------------------
# Initialization stuff

# Check to make sure that the lock file does not already exist.
if [ -f $pidfile ]; then
  echo 2>&1 "Process lock file " $pidfile " already exists, exiting."
  exit -1
fi

# Create the lock file
echo $$ > $pidfile

# If our process is killed, remove the lock file first
trap '/bin/rm -f $pidfile' EXIT HUP TERM XCPU KILL

# Number of seconds to sleep before checking again.
N=60

#----------------------------------------------------------------------------
# Make storage directory

# Once we process a checkpoint or plotfile, we will move it into the output 
# directory.  This then hides them from the script, so if the system
# later purges file and their .processed counterparts, 
# we don't overwrite our archived data with a tarred empty
# directory structure. 

if [ ! -d output ]; then
  mkdir output
fi

#----------------------------------------------------------------------------
# The processing function

# Process Files.  Once a plotfile is successfully processed, we will output
# a file pltXXXXX.processed (checkpoint files are only archived, with a
# chkXXXXX.processed file appearing once the archiving is successful).  
# Subsequent invocations of this routine will skip over any plotfiles or
# checkpoint files that have a corresponding .processed file.


function process_files
{

  if [ ! -f $pidfile ]; then
    echo "Process: $pidfile has been removed, exiting."
    exit
  fi

  archive_all $cwd

}


#----------------------------------------------------------------------------
# The main function. Loop, waiting for plt and chk directories to appear.

while true
do
  process_files
  sleep $N
done
