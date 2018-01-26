# probin.sh: functions for working with probin files.

# Given a probin variable name in argument 1 and 
# a directory in argument 2, replace the probin variable 
# with the value of argument 1 using indirect references.

function replace_probin_var {

    # Check if the variable name and directory exist,
    # as well as the probin file in that directory.

    if [ ! -z $1 ]; then
	var=$1
    else
	echo "No variable found in arguments list for replace_probin_var; exiting."
	return
    fi

    if [ -z $dir ]; then
	echo "No directory found in replace_probin_var; exiting."
	return
    fi

    if [ -z $probin ]; then
	probin=probin
    fi

    if [ ! -e $dir/$probin ]; then
	echo "No probin file exists in directory "$dir"; exiting."
	return
    fi

    # We only replace in cases that there's a space both before and 
    # after the variable, signifying a full word match.
    # This will require that your probin variables have
    # at least one space behind them and after them.
    # The equals sign cannot be adjacent to the name.

    if (grep -q "[[:space:]]$var[[:space:]]*=" $dir/$probin); then
	sed -i "s/ $var .*=.*/ $var = ${!var}/g" $dir/$probin
        return
    fi

    # Also allow for the possibility of a newly added probin variable.
    # We will detect this using syntax where the variable must be
    # prepended with probin_<namelist>_, where <namelist> is the name
    # of the specific Fortran namelist you want to add to.

    if [[ $var == "probin_"* ]]; then

        # Extract the namelist
        namelist=$(echo $var | cut -d'_' -f 2)

        # Extract the actual name of the variable
        actual_var=$(echo $var | cut -d'_' -f 3-)

        # Insert the variable into the probin file inside the matching namelist
        sed -i "/&$namelist/ a \ \ $actual_var = ${!var}" $dir/$probin

        # Insert newline before it for aesthetics
        sed -i "/&$namelist/ a \ \ " $dir/$probin

        return

    fi

}
