#!/bin/bash

source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

    if [ -z $ncell ]; then
        echoerr "ncell not set; exiting."
        exit
    fi

    # The following assumes we are using Lassen.

    if [ $MACHINE == "PERLMUTTER" ]; then

        nprocs=4
        walltime="120"

        # Set up the geometry and gridding appropriately.

        amr__n_cell="$ncell"
        amr__ref_ratio="4"
        amr__blocking_factor="64"
        amr__max_grid_size="1048576"

        if [ ! -z $refinement ]; then

            if [ $refinement -eq 1 ]; then
                amr__max_level=0
                amr__ref_ratio="2"
            elif [ $refinement -eq 2 ]; then
                amr__max_level=1
                amr__ref_ratio="2"
            elif [ $refinement -eq 4 ]; then
                amr__max_level=1
                amr__ref_ratio="4"
            elif [ $refinement -eq 8 ]; then
                amr__max_level=2
                amr__ref_ratio="4 2"
            elif [ $refinement -eq 16 ]; then
                amr__max_level=2
                amr__ref_ratio="4 4"
            elif [ $refinement -eq 32 ]; then
                amr__max_level=3
                amr__ref_ratio="4 4 2"
            elif [ $refinement -eq 64 ]; then
                amr__max_level=3
                amr__ref_ratio="4 4 4"
            elif [ $refinement -eq 128 ]; then
                amr__max_level=4
                amr__ref_ratio="4 4 4 2"
            elif [ $refinement -eq 256 ]; then
                amr__max_level=4
                amr__ref_ratio="4 4 4 4"
            elif [ $refinement -eq 512 ]; then
                amr__max_level=5
                amr__ref_ratio="4 4 4 4 2"
            elif [ $refinement -eq 1024 ]; then
                amr__max_level=5
                amr__ref_ratio="4 4 4 4 4"
            elif [ $refinement -eq 2048 ]; then
                amr__max_level=6
                amr__ref_ratio="4 4 4 4 4 2"
            elif [ $refinement -eq 4096 ]; then
                amr__max_level=6
                amr__ref_ratio="4 4 4 4 4 4"
            elif [ $refinement -eq 8192 ]; then
                amr__max_level=7
                amr__ref_ratio="4 4 4 4 4 4 2"
            elif [ $refinement -eq 16384 ]; then
                amr__max_level=7
                amr__ref_ratio="4 4 4 4 4 4 4"
            elif [ $refinement -eq 32768 ]; then
                amr__max_level=8
                amr__ref_ratio="4 4 4 4 4 4 4 2"
            elif [ $refinement -eq 65536 ]; then
                amr__max_level=8
                amr__ref_ratio="4 4 4 4 4 4 4 4"
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
exec_dir=$CASTRO_HOME/Exec/science/Detonation

use_first_castro_ex="1"

# Get the right inputs and probin files.

inputs="inputs-collision"
probin="probin-collision"

# Don't subcycle.

amr__subcycling_mode="None"
castro__use_post_step_regrid=0

# Limit where we refine.

amr__refine__dxnucerr__value_greater="1.e-8"

# The flag we will use to determine whether to run the job.

to_run=1

results_dir="results"

ncell_list="128 256 512 1024 2048 4096 8192 16384 32768 65536 131072"
ncell_list="128 256 512 1024 2048 4096 8192 16384 32768 65536"

burning_mode="1"

ofrac_list="0.00d0 0.45d0 0.50d0"
ofrac_list="0.50d0"

network_list="aprox13 He-C-Fe-group"

for NETWORK_DIR in $network_list
do

    executable_keyword=$NETWORK_DIR

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

                dir=$results_dir/$NETWORK_DIR/ofrac$ofrac/n$ncell/r$refinement
                set_run_opts

                if [ $to_run -eq 1 ]; then
                    run
                fi

                to_run=1

            done # refinement

        done # ncell

    done # ofrac

done # network
