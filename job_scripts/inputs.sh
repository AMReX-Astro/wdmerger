# inputs.sh: functions for working with inputs files.

# Obtain the value of a variable in the main inputs file.
# Optionally we can provide a second argument to look 
# in the inputs file in another directory.

function get_inputs_var {

    # Check to make sure that the variable name was given to us.

    if [ ! -z $1 ]; then
	var_name=$1
    else
	return
    fi

    if [ ! -z $2 ]; then
	directory=$2
    elif [ ! -z $dir ]; then
	directory=$dir
    else
	echo "No directory found in replace_inputs_var; exiting."
	return
    fi

    if [ -z $inputs ]; then
	inputs=inputs
    fi

    if [ ! -e $directory/$inputs ]; then
	echo "No inputs file exists in directory "$dir"; exiting."
	return
    fi

    inputs_var_name=$(fix_inputs_var_name $var_name)

    var=""

    # To get the value of the inputs variable, 
    # we need to get everything after the equals sign.

    if (grep -q "$inputs_var_name[[:space:]]*=" $directory/$inputs); then

	var=$(grep "$inputs_var_name" $directory/$inputs | awk -F"=" '{print $2}')

    fi

    echo $var

}



# If a variable corresponds to a valid CASTRO inputs variable,
# return the name of this variable with the proper formatting:
# replace the underscore after the namespace with a period.

function fix_inputs_var_name {

    if [ ! -z $1 ]; then
	var=$1
    else
	return
    fi

    namespace=$(echo $var | cut -d _ -f 1)

    inputs_var_name=""

    # We need to be careful with stop_time and max_step,
    # since these don't have an associated namespace.
    # Otherwise, check on the namespaces we know from the inputs file.

    if ([ $var == "stop_time" ] || [ $var == "max_step" ])
    then

	inputs_var_name=$var	    

    elif ( [ $namespace == "amr" ] || [ $namespace == "castro" ] || 
	   [ $namespace == "geometry" ] || [ $namespace == "gravity" ] || 
	   [ $namespace == "mg" ] || [ $namespace == "DistributionMapping" ] ||
	   [ $namespace == "fab" ] )
    then

	# Remove the namespace from the variable, then
	# replace it with the correct period.

	inputs_var_end=$(echo ${var#$namespace\_})
	inputs_var_name="$namespace.$inputs_var_end"

    fi

    echo $inputs_var_name

}



# Given an inputs variable name (in bash-friendly form) in argument 1 
# and a directory in argument 2, replace the inputs variable 
# with the value of argument 1 using indirect references.

function replace_inputs_var {

    # Check if the variable name and directory exist,
    # as well as the inputs file in that directory.

    if [ ! -z $1 ]; then
	var=$1
    else
	echo "No variable found in arguments list for replace_inputs_var; exiting."
	return
    fi

    if [ -z $dir ]; then
	echo "No directory found in replace_inputs_var; exiting."
	return
    fi

    if [ -z $inputs ]; then
	inputs=inputs
    fi

    if [ ! -e $dir/$inputs ]; then
	echo "No inputs file exists in directory "$dir"; exiting."
	return
    fi

    inputs_var_name=$(fix_inputs_var_name $var)

    if [ -z $inputs_var_name ]; then
	return
    fi

    # Make sure this inputs variable actually exists in the inputs file.
    # The -q option to grep means be quiet and only 
    # return a flag indicating whether $var is in inputs.
    # The [[:space:]]* tells grep to ignore any whitespace between the 
    # variable and the equals sign. See:
    # http://www.linuxquestions.org/questions/programming-9/grep-ignoring-spaces-or-tabs-817034/

    if (grep -q "$inputs_var_name[[:space:]]*=" $dir/$inputs)
    then

	# OK, so the parameter name does exist in the inputs file;
	# now we just need to get its value in there. We can do this using
	# bash indirect references -- if $var is a variable name, then
	# ${!var} is the value held by that variable. See:
	# http://www.tldp.org/LDP/abs/html/bashver2.html#EX78
	# http://stackoverflow.com/questions/10955479/name-of-variable-passed-to-function-in-bash

	old_string=$(grep "$inputs_var_name[[:space:]]" $dir/$inputs)
	new_string="$inputs_var_name = ${!var}"

	sed -i "s/$inputs_var_name.*/$new_string/g" $dir/$inputs

    # If the inputs variable doesn't exist in the inputs file, but we got to this point,
    # it had a legitimate BoxLib or CASTRO namespace. We assume the user wants to add 
    # this variable to the inputs file, so we'll add it on at the end of the file.

    else

	string="$inputs_var_name = ${!var}"
	echo $string >> $dir/$inputs

    fi

}
