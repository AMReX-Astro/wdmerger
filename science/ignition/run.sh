#!/bin/bash

source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

    if [ -z $ncell ]; then
        echoerr "ncell not set; exiting."
        exit
    fi

    # The following assumes we are using Lassen.

    if [ $MACHINE == "LASSEN" ]; then

        queue="pbatch"

        threads_per_task=1

        nodes="1"
        walltime="2:00:00"

        # Set up the geometry and gridding appropriately.

        amr_n_cell="$ncell"
        amr_ref_ratio="4"
        amr_blocking_factor="64"
        amr_max_grid_size="1048576"

        if [ ! -z $refinement ]; then

            # Only refinement factors of 4 are allowed.

            if [ $refinement -eq 1 ]; then
                amr_max_level=0
                amr_ref_ratio="2"
            elif [ $refinement -eq 2 ]; then
                amr_max_level=1
                amr_ref_ratio="2"
            elif [ $refinement -eq 4 ]; then
                amr_max_level=1
                amr_ref_ratio="4"
            elif [ $refinement -eq 8 ]; then
                amr_max_level=2
                amr_ref_ratio="4 2"
            elif [ $refinement -eq 16 ]; then
                amr_max_level=2
                amr_ref_ratio="4 4"
            elif [ $refinement -eq 32 ]; then
                amr_max_level=3
                amr_ref_ratio="4 4 2"
            elif [ $refinement -eq 64 ]; then
                amr_max_level=3
                amr_ref_ratio="4 4 4"
            elif [ $refinement -eq 128 ]; then
                amr_max_level=4
                amr_ref_ratio="4 4 4 2"
            elif [ $refinement -eq 256 ]; then
                amr_max_level=4
                amr_ref_ratio="4 4 4 4"
            elif [ $refinement -eq 512 ]; then
                amr_max_level=5
                amr_ref_ratio="4 4 4 4 2"
            elif [ $refinement -eq 1024 ]; then
                amr_max_level=5
                amr_ref_ratio="4 4 4 4 4"
            elif [ $refinement -eq 2048 ]; then
                amr_max_level=6
                amr_ref_ratio="4 4 4 4 4 2"
            elif [ $refinement -eq 4096 ]; then
                amr_max_level=6
                amr_ref_ratio="4 4 4 4 4 4"
            elif [ $refinement -eq 8192 ]; then
                amr_max_level=7
                amr_ref_ratio="4 4 4 4 4 4 2"
            elif [ $refinement -eq 16384 ]; then
                amr_max_level=7
                amr_ref_ratio="4 4 4 4 4 4 4"
            elif [ $refinement -eq 32768 ]; then
                amr_max_level=8
                amr_ref_ratio="4 4 4 4 4 4 4 2"
            elif [ $refinement -eq 65536 ]; then
                amr_max_level=8
                amr_ref_ratio="4 4 4 4 4 4 4 4"
            else
                echo "Unknown refinement factor: "$refinement"; exiting."
                exit
            fi

        fi

    else

        echoerr "This machine is not set up for this job."

    fi

}

DIM=1

# Specify the problem directory.

compile_dir=$CASTRO_HOME/Exec/science/Detonation
problem_dir=$CASTRO_HOME/Exec/science/Detonation

use_first_castro_ex="1"

# Use the aprox13 network for all the tests.

NETWORK_DIR="aprox13"

# Get the right inputs and probin files.

inputs="inputs-collision"
probin="probin-collision"

# The flag we will use to determine whether to run the job.

to_run=1

results_dir="results"

ncell_list="128 256 512 1024 2048 4096 8192 16384 32768 65536 131072"

burning_mode="1"

ofrac_list="0.00d0 0.45d0 0.50d0"

for ofrac in $ofrac_list
do

    for ncell in $ncell_list
    do

        refinement_list="1"

        if [ $ncell -eq 128 ]; then
            refinement_list="1"
        elif [ $ncell -eq 256 ]; then
            refinement_list="1"
        elif [ $ncell -eq 512 ]; then
            refinement_list="1"
        elif [ $ncell -eq 1024 ]; then
            refinement_list="1"
        elif [ $ncell -eq 2048 ]; then
            refinement_list="1"
        elif [ $ncell -eq 4096 ]; then
            refinement_list="1"
        elif [ $ncell -eq 8192 ]; then
            refinement_list="1"
        elif [ $ncell -eq 16384 ]; then
            refinement_list="1"
        elif [ $ncell -eq 32768 ]; then
            refinement_list="1"
        elif [ $ncell -eq 65536 ]; then
            refinement_list="1"
        elif [ $ncell -eq 131072 ]; then
            refinement_list="1 2 4 8 16 32 64 128 256 512"
        fi

        for refinement in $refinement_list
        do

            dir=$results_dir/ofrac$ofrac/n$ncell/r$refinement
            set_run_opts

            if [ $to_run -eq 1 ]; then
                run
            fi

            to_run=1

        done # refinement

    done # ncell

done # ofrac
