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


  # Plotfiles

  # Take all but the final plt file -- we want to ensure they're completely 
  # written to disk.  Strip out any tar files that are lying around as well 
  # as pltXXXXX.processed files.  We restrict the find command to a depth of 
  # 1 to avoid catching any already-processed files in the plotfiles/
  # directory
  pltlist5=$(find . -maxdepth 1 -type d -name "${plt_prefix}?????"  -print | sort)
  pltlist6=$(find . -maxdepth 1 -type d -name "${plt_prefix}??????" -print | sort)

  pltlist="$pltlist5 $pltlist6"

  if [ "$pltlist" ]; then
    nl=$(echo "$pltlist" | wc -l)
    nl=$(expr $nl - 1)
    if [ $nl -eq 0 ]; then
      pltlist=""
    else
      pltlist=$(echo "$pltlist" | head -$nl)
    fi
  fi

  for dir in ${pltlist}
  do
    if [ -d ${dir} ]; then

      # Only work on the file if there is not a .processed file in the
      # main directory or the plotfiles/ directory.

      if [ ! -f ${dir}.processed ] && [ ! -f output/${dir}.processed ]; then
        # Store the file on the archive system
        archive $cwd/$dir
      fi

    fi 
  done


  # Checkpoint files

  # Take all but the final chk file -- we want to ensure they're completely 
  # written to disk.  Strip out any tar files that are lying around as well 
  # as chkXXXXX.processed files.  We restrict the find command to a depth of
  # 1 to avoid catching any already-processed files in the checkfiles/ 
  # directory
  chklist5=$(find . -maxdepth 1 -type d -name "${chk_prefix}?????"  -print | sort)
  chklist6=$(find . -maxdepth 1 -type d -name "${chk_prefix}??????" -print | sort)

  chklist="$chklist5 $chklist6"

  if [ "$chklist" ]; then
    nl=$(echo "$chklist" | wc -l)
    nl=$(expr $nl - 1)
    if [ $nl -eq 0 ]; then
      chklist=""
    else
      chklist=$(echo "$chklist" | head -$nl)
    fi
  fi

  for dir in ${chklist}
  do
    if [ -d ${dir} ]; then

      if [ ! -f ${dir}.processed ] && [ ! -f checkfiles/${dir}.processed ]; then
	# Store the file on the archive system
	archive $cwd/$dir
      fi

    fi
  done

  # Diagnostic files

  # Here our strategy will be different: we want to make a copy of the diagnostic file
  # and append the date to that, and then archive the tagged copy, because the base file
  # needs to remain there for the duration of the run.

  diag_files=$(find . -maxdepth 1 -name "*_diag.out")
  datestr=$(date +"%Y%m%d_%H%M_%S_")

  for file in $diag_files
  do
      filebasename=$(basename $file)
      archivefile=$datestr$filebasename
      cp $file $archivefile
      archive $cwd/$archivefile
  done

}


#----------------------------------------------------------------------------
# The main function. Loop, waiting for plt and chk directories to appear.

while true
do
  process_files
  sleep $N
done
