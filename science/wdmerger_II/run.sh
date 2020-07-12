#!/bin/bash

source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

    if [ -z $ncell ]; then
	echoerr "ncell not set; exiting."
	exit
    fi

    to_run=1

    # Set the processor count and walltime.

    if [ "$MACHINE" == "SUMMIT" ]; then

        queue="batch"
        nprocs=48

        if [ $reactions_refinement -eq 2 ]; then
            nprocs=96
        elif [ $reactions_refinement -eq 4 ]; then
            nprocs=192
        elif [ $reactions_refinement -eq 8 ]; then
            nprocs=384
        elif [ $reactions_refinement -eq 16 ]; then
            nprocs=768
        fi

        walltime="2:00:00"

        amr_blocking_factor="128 128 256 256 512"
        amr_max_grid_size="256 512 1024"

    else

        echoerr "This machine is not set up for this job."

    fi

    # Set up the geometry appropriately.

    amr_n_cell="$(echo "$ncell / 2" | bc) $ncell"
    geometry_prob_lo="0.0 $prob_lo"
    geometry_prob_hi="$prob_hi $prob_hi"

    if [ ! -z $reactions_refinement ]; then

        if [ $reactions_refinement -eq 1 ]; then

            castro_reactions_max_solve_level="0"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=0
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=1
                    amr_ref_ratio="2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=1
                    amr_ref_ratio="4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=2
                    amr_ref_ratio="4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=2
                    amr_ref_ratio="4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 2"
                elif [ $hydro_refinement -eq 1024 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 4"
                elif [ $hydro_refinement -eq 2048 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 4096 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 8192 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 16384 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 32768 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 65536 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 131072 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 262144 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 524288 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1048576 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 2097152 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        elif [ $reactions_refinement -eq 2 ]; then

            castro_reactions_max_solve_level="1"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=1
                    amr_ref_ratio="2"
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=2
                    amr_ref_ratio="2 2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=2
                    amr_ref_ratio="2 4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=3
                    amr_ref_ratio="2 4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=3
                    amr_ref_ratio="2 4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=4
                    amr_ref_ratio="2 4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=4
                    amr_ref_ratio="2 4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=5
                    amr_ref_ratio="2 4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=5
                    amr_ref_ratio="2 4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=6
                    amr_ref_ratio="2 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1024 ]; then
                    amr_max_level=6
                    amr_ref_ratio="2 4 4 4 4 4"
                elif [ $hydro_refinement -eq 2048 ]; then
                    amr_max_level=7
                    amr_ref_ratio="2 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 4096 ]; then
                    amr_max_level=7
                    amr_ref_ratio="2 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 8192 ]; then
                    amr_max_level=8
                    amr_ref_ratio="2 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 16384 ]; then
                    amr_max_level=8
                    amr_ref_ratio="2 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 32768 ]; then
                    amr_max_level=9
                    amr_ref_ratio="2 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 65536 ]; then
                    amr_max_level=9
                    amr_ref_ratio="2 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 131072 ]; then
                    amr_max_level=10
                    amr_ref_ratio="2 4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 262144 ]; then
                    amr_max_level=10
                    amr_ref_ratio="2 4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 524288 ]; then
                    amr_max_level=11
                    amr_ref_ratio="2 4 4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1048576 ]; then
                    amr_max_level=11
                    amr_ref_ratio="2 4 4 4 4 4 4 4 4 4 4"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        elif [ $reactions_refinement -eq 4 ]; then

            castro_reactions_max_solve_level="1"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=1
                    amr_ref_ratio="4"
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=2
                    amr_ref_ratio="4 2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=2
                    amr_ref_ratio="4 4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1024 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 2048 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 4096 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 8192 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 16384 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 32768 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 65536 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 131072 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 262144 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 524288 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        elif [ $reactions_refinement -eq 8 ]; then

            castro_reactions_max_solve_level="2"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=2
                    amr_ref_ratio="4 2"
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 2 2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 2 4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 2 4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 2 4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 2 4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 2 4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 2 4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 2 4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 2 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1024 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 2 4 4 4 4 4"
                elif [ $hydro_refinement -eq 2048 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 2 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 4096 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 2 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 8192 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 2 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 16384 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 2 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 32768 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 2 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 65536 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 2 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 131072 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 2 4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 262144 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 2 4 4 4 4 4 4 4 4 4"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        elif [ $reactions_refinement -eq 16 ]; then

            castro_reactions_max_solve_level="2"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=2
                    amr_ref_ratio="4 4"
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1024 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 2048 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 4096 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 8192 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 16384 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 32768 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 65536 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 131072 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        elif [ $reactions_refinement -eq 32 ]; then

            castro_reactions_max_solve_level="3"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 2"
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 2 2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 2 4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 2 4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 2 4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 2 4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 2 4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 2 4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 2 4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 2 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1024 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 2 4 4 4 4 4"
                elif [ $hydro_refinement -eq 2048 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 2 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 4096 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 2 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 8192 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 2 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 16384 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 2 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 32768 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 2 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 65536 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 2 4 4 4 4 4 4 4 4"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        elif [ $reactions_refinement -eq 64 ]; then

            castro_reactions_max_solve_level="3"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=3
                    amr_ref_ratio="4 4 4"
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1024 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 2048 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 4096 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 8192 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 16384 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 32768 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        elif [ $reactions_refinement -eq 128 ]; then

            castro_reactions_max_solve_level="4"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 2"
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 2 2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 2 4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 2 4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 2 4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 2 4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 2 4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 2 4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 2 4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 2 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1024 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 2 4 4 4 4 4"
                elif [ $hydro_refinement -eq 2048 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 2 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 4096 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 2 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 8192 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 2 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 16384 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 2 4 4 4 4 4 4 4"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        elif [ $reactions_refinement -eq 256 ]; then

            castro_reactions_max_solve_level="4"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=4
                    amr_ref_ratio="4 4 4 4"
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1024 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 2048 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 4096 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 8192 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        elif [ $reactions_refinement -eq 512 ]; then

            castro_reactions_max_solve_level="5"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 2"
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 2 2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 2 4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 2 4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 2 4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 2 4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 2 4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 2 4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 2 4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 2 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1024 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 2 4 4 4 4 4"
                elif [ $hydro_refinement -eq 2048 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 4 2 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 4096 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 4 2 4 4 4 4 4 4"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        elif [ $reactions_refinement -eq 1024 ]; then

            castro_reactions_max_solve_level="5"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=5
                    amr_ref_ratio="4 4 4 4 4"
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1024 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 2048 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        elif [ $reactions_refinement -eq 2048 ]; then

            castro_reactions_max_solve_level="6"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 2 2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 2 4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 2 4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 2 4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 2 4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 2 4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 2 4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 2 4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 4 4 2 4 4 4 4 2"
                elif [ $hydro_refinement -eq 1024 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 4 4 2 4 4 4 4 4"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        elif [ $reactions_refinement -eq 4096 ]; then

            castro_reactions_max_solve_level="6"

            if [ ! -z $hydro_refinement ]; then

                if [ $hydro_refinement -eq 1 ]; then
                    amr_max_level=6
                    amr_ref_ratio="4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 2 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 4 ]; then
                    amr_max_level=7
                    amr_ref_ratio="4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 8 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 16 ]; then
                    amr_max_level=8
                    amr_ref_ratio="4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 32 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 64 ]; then
                    amr_max_level=9
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 128 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                elif [ $hydro_refinement -eq 256 ]; then
                    amr_max_level=10
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                elif [ $hydro_refinement -eq 512 ]; then
                    amr_max_level=11
                    amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                else
                    echo "Unknown refinement factor: "$hydro_refinement"; exiting."
                    exit
                fi

            fi

        fi

        # Only solve gravity on the coarse level. The additional accuracy
        # from solving on the fine levels would provide no benefit on this problem.

        gravity_max_solve_level=0

        # Tag all zones in the star for maximum refinement.

        max_stellar_tagging_level=$amr_max_level

    fi

}

function copy_checkpoint() {

    # Copy the initial checkpoint to the directory if it's not there yet.

    if [ -e "$dir/jobIsDone" ]; then
        return
    fi

    job_flag=$(is_job_running $dir)

    if [ $job_flag -ne 1 ]; then

        done_flag=$(is_dir_done)

        if [ $done_flag -ne 1 ]; then

            no_submit=1

            if [ $to_run -eq 1 ]; then
                run
            fi

            unset no_submit

            # Skip the saved checkpoint copy if there's already a checkpoint
            # in the target directory (either it's the saved checkpoint, or it
            # is a later one because we've already started the run).

            checkpoint=$(get_last_checkpoint $dir)

            if [ "$checkpoint" == "" ]; then

                # By default the starting checkpoint is the last checkpoint
                # from the start directory, but this can be overridden if
                # the user prefers a different starting checkpoint.

                if [ -z $start_checkpoint ]; then
                    start_checkpoint=$(get_last_checkpoint $start_dir)
                fi

                if [ -z $checkpoint_dir ]; then
                    checkpoint_dir=$start_dir
                fi

                echo "Copying initial checkpoint to" $dir"."

                if [ -d "$checkpoint_dir/$start_checkpoint" ]; then

                    cp -r "$checkpoint_dir/$start_checkpoint" $dir
                    rm -f "$dir/$start_checkpoint/jobIsDone"

                    # Also, copy the diagnostic output files
                    # so there is an apparently continuous record.

                    cp $checkpoint_dir/*diag.out $dir

                    # Strip out any data in the diagnostics files
                    # that is after the copied checkpoint, to
                    # ensure apparent continuity.

                    chk_step=$(echo $start_checkpoint | cut -d"k" -f2)
                    chk_step=$(echo $chk_step | bc) # Strip leading zeros

                    for diag in $(find $dir -name "*diag.out");
                    do
                        sed -i "/\ $chk_step \ /q" $diag
                    done

                else
                    echoerr "No initial checkpoint available, exiting."
                    exit
                fi

            fi

        fi

    fi

    if [ ! -z $start_checkpoint ]; then
        unset start_checkpoint
    fi

    if [ ! -z $checkpoint_dir ]; then
        unset checkpoint_dir
    fi

}

# Specify the problem directory.

problem_dir=$CASTRO_HOME/Exec/science/wdmerger

use_first_castro_ex="1"

# Needed for the makefile: we want to compile in 2D for these tests.

DIM="2"

# Use the aprox13 network for all following tests.

NETWORK_DIR="aprox13"

# Get the right inputs and probin files.

inputs="inputs_2d"
probin="probin"

# Abort if we run out of GPU memory.

amrex_abort_on_out_of_gpu_memory="1"

# Variables we need to set up the collision.

problem="0"
collision_separation="2.0"
collision_impact_parameter="0.0"

# Disable refinement based on non-burning criteria.
# These can be re-enabled later for specific tests.

max_stellar_tagging_level="0"
max_temperature_tagging_level="0"
max_center_tagging_level="0"

# Allow first-order interpolations to fine levels.

castro_state_interp_order="1"

# Make a full plotfile every second.

amr_plot_per="1.0"
amr_plot_vars="ALL"
amr_derive_plot_vars="ALL"

# Make small plotfiles rapidly.

amr_small_plot_per="0.1"
amr_small_plot_vars="density Temp rho_e rho_c12 rho_o16 rho_si28 rho_ni56 enuc"
amr_derive_small_plot_vars="pressure soundspeed x_velocity y_velocity t_sound_t_enuc"

castro_plot_per_is_exact="0"
castro_small_plot_per_is_exact="0"

# Save checkpoints every second.

amr_check_per="1.0"

# Initial timestep shortening factor.

castro_init_shrink="0.1"

# CFL number.

castro_cfl="0.8"

# Maximum number of subcycles.

castro_max_subcycles="128"

# Use post-timestep regrids.

castro_use_post_step_regrid="1"

# Enable reactions.

castro_do_react="1"

# Disable rotation.

castro_do_rotation="0"

# Many of the collision papers in the literature use an equal
# C/O ratio  by mass in the initial white dwarfs. We will do
# this too for comparison purposes.

co_wd_c_frac="0.5d0"
co_wd_o_frac="0.5d0"

# Allow the timestep to change by up to 25% per advance.

castro_change_max="1.25"

# Some defaults.

ncell_default="1024"
dtnuc_e_default="1.e200"
dtnuc_X_default="1.e200"
mass_P_default="0.64"
mass_S_default="0.64"

ncell=$ncell_default
mass_P=$mass_P_default
mass_S=$mass_S_default
castro_dtnuc_e=$dtnuc_e_default
castro_dtnuc_X=$dtnuc_X_default
castro_react_T_min="5.0e8"
castro_react_rho_min="1.0e6"

# The flag we will use to determine whether to run the job.

to_run=1

# Simulate 0.64 + 0.64 M_solar WDs.

mass_list="0.64"

# Run four helium shell mass cases. These are the same
# values from Holcomb and Kushnir (2016), so that we
# can make a direct comparison (see Figure 2).

helium_shell_mass_list="0.00 0.04 0.08 0.16"

stop_time="4.0"
prob_lo="-5.12e9"
prob_hi="5.12e9"

reactions_refinement_list="1 2 4 8 16"

hydro_refinement_list="1 4 16 64"

for mass in $mass_list
do

    mass_P=$mass
    mass_S=$mass

    for helium_shell_mass in $helium_shell_mass_list
    do

        co_wd_he_shell_mass=$helium_shell_mass

        for reactions_refinement in $reactions_refinement_list
        do

            effective_ncell=$(echo "$ncell * $reactions_refinement" | bc)

            for hydro_refinement in $hydro_refinement_list
            do

                dir=$results_dir/m$mass/he$helium_shell_mass/reactions_r$reactions_refinement/hydro_r$hydro_refinement
                set_run_opts

                if [ $to_run -eq 1 ]; then
                    run
                fi

            done

        done

    done

done
