# math.sh: functions for doing useful mathematical operations.

# Determine the median of a space-separated list of floating point numbers.

function median {

    if [ ! -z "$1" ]; then
	numbers=$1
    else
	return
    fi

    # Determine the number of values we have.

    Nvals=0

    for num in $numbers
    do
      Nvals=$((Nvals+1))
    done

    # Sort them numerically. To use sort, we have to 
    # convert the space-separated list into a newline-separated list.

    numbers=$(echo $numbers | tr " " "\n" | sort -n)

    # Now, find the median value. We have N values. If N is odd, we use the middle value;
    # if N is even, we compute the average of the two middle values.

    idx=1

    median=0.0

    for num in $numbers
    do
	if (( $Nvals % 2 == 0 && ( $idx == $Nvals / 2 || $idx == $Nvals / 2 + 1 ) )); then
	    median=$(echo "$median + ($num / 2.0)" | bc -l)
	elif (( $idx == $Nvals / 2 + 1 )); then
	    median=$(echo "$median + $num" | bc -l)
	fi
	idx=$((idx+1))
    done

    echo $median

}



# Convert a time in (-)HH:MM:SS format to a number of seconds.

function hours_to_seconds {

    if [ ! -z $1 ]; then
	hours=$1
    else
	return
    fi

    # If we start off with a negative sign,
    # create a variable equal to -1 and negate
    # the whole thing at the end.

    if [ ${hours#-} != $hours ]; then
	negative_fac="-1"
	hours=${hours#-}
    else
	negative_fac="1"
    fi

    seconds=$(echo $hours | awk -F: '{ print $1 * 3600 + $2 * 60 + $3 }')
    seconds=$(echo "$seconds * $negative_fac" | bc -l)

    echo $seconds

}
