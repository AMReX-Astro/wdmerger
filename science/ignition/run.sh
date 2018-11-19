#!/bin/bash

source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

    if [ -z $ncell ]; then
        echoerr "ncell not set; exiting."
        exit
    fi

    # The following assumes we are using Titan.

    if [ $MACHINE == "TITAN" ]; then

        queue="batch"

        threads_per_task=4

        if [ $DIM -eq 1 ]; then

            nprocs="16"
            walltime="2:00:00"

        else

            if   [ $ncell -eq 16384 ]; then
                nprocs="2048"
                walltime="6:00:00"
            elif [ $ncell -eq 8192 ]; then
                nprocs="2048"
                walltime="6:00:00"
            elif [ $ncell -eq 4096 ]; then
                nprocs="1024"
                walltime="2:00:00"

                if [ ! -z $stellar_refinement ]; then
                    if [ $stellar_refinement -eq 8 ]; then
                        nprocs="2048"
                        walltime="6:00:00"
                    elif [ $stellar_refinement -eq 16 ]; then
                        nprocs="2048"
                        walltime="6:00:00"
                    fi
                fi
            elif [ $ncell -eq 2048 ]; then
                nprocs="512"
                walltime="2:00:00"

                if [ ! -z $stellar_refinement ]; then
                    if [ $stellar_refinement -eq 16 ]; then
                        nprocs="1024"
                        walltime="2:00:00"
                    fi
                fi
            elif [ $ncell -eq 1024 ]; then
                nprocs="256"
                walltime="2:00:00"

                if [ ! -z $stellar_refinement ]; then
                    if [ $stellar_refinement -eq 16 ]; then
                        nprocs="512"
                        walltime="2:00:00"
                    fi
                fi
            elif [ $ncell -eq 512 ]; then
                nprocs="128"
                walltime="2:00:00"

                if [ ! -z $stellar_refinement ]; then
                    if [ $stellar_refinement -eq 8 ]; then
                        nprocs="512"
                        walltime="2:00:00"
                    fi
                fi
            elif [ $ncell -eq 256 ]; then
                nprocs="32"
                walltime="2:00:00"

                if [ ! -z $stellar_refinement ]; then
                    if [ $stellar_refinement -eq 16 ]; then
                        nprocs="128"
                        walltime="2:00:00"
                    elif [ $stellar_refinement -eq 64 ]; then
                        nprocs="1024"
                        walltime="2:00:00"
                    fi
                fi
            elif [ $ncell -eq 128 ]; then
                nprocs="32"
                walltime="2:00:00"

                if [ ! -z $stellar_refinement ]; then
                    if [ $stellar_refinement -eq 1 ]; then
                        if [ ! -z $refinement ]; then
                            if [ $refinement -eq 4 ]; then
                                nprocs="32"
                                walltime="2:00:00"
                            elif [ $refinement -eq 16 ]; then
                                nprocs="64"
                                walltime="2:00:00"
                            elif [ $refinement -eq 64 ]; then
                                nprocs="256"
                                walltime="2:00:00"
                            elif [ $refinement -eq 256 ]; then
                                nprocs="1024"
                                walltime="2:00:00"
                            fi
                        fi
                    elif [ $stellar_refinement -eq 4 ]; then
                        nprocs="32"
                        walltime="2:00:00"
                    elif [ $stellar_refinement -eq 16 ]; then
                        nprocs="128"
                        walltime="2:00:00"
                    elif [ $stellar_refinement -eq 64 ]; then
                        nprocs="1024"
                        walltime="2:00:00"
                    fi
                fi
            elif [ $ncell -eq 64 ]; then
                nprocs="16"
                walltime="2:00:00"

                if [ ! -z $stellar_refinement ]; then
                    if [ $stellar_refinement -eq 4 ]; then
                        nprocs="64"
                        walltime="2:00:00"
                    elif [ $stellar_refinement -eq 8 ]; then
                        nprocs="128"
                        walltime="2:00:00"
                    elif [ $stellar_refinement -eq 16 ]; then
                        nprocs="256"
                        walltime="2:00:00"
                    elif [ $stellar_refinement -eq 32 ]; then
                        nprocs="256"
                        walltime="2:00:00"
                    elif [ $stellar_refinement -eq 64 ]; then
                        nprocs="256"
                        walltime="2:00:00"
                    elif [ $stellar_refinement -eq 128 ]; then
                        nprocs="512"
                        walltime="2:00:00"
                    elif [ $stellar_refinement -eq 256 ]; then
                        nprocs="1024"
                        walltime="2:00:00"
                    elif [ $stellar_refinement -eq 512 ]; then
                        nprocs="1024"
                        walltime="2:00:00"
                    elif [ $stellar_refinement -eq 1024 ]; then
                        nprocs="2048"
                        walltime="6:00:00"
                    elif [ $stellar_refinement -eq 2048 ]; then
                        nprocs="4096"
                        walltime="6:00:00"
                    fi
                fi
            else
                echoerr "Unknown number of cells per dimension."
            fi

        fi




        # Set up the geometry and gridding appropriately.

        if [ $DIM -eq 1 ]; then

            amr_n_cell="$ncell"
            amr_ref_ratio="4"

            if [ ! -z $refinement ]; then

                # Only refinement factors of 4 are allowed.

                if [ $refinement -eq 1 ]; then
                    amr_max_level=0
                elif [ $refinement -eq 4 ]; then
                    amr_max_level=1
                elif [ $refinement -eq 16 ]; then
                    amr_max_level=2
                elif [ $refinement -eq 64 ]; then
                    amr_max_level=3
                elif [ $refinement -eq 256 ]; then
                    amr_max_level=4
                elif [ $refinement -eq 1024 ]; then
                    amr_max_level=5
                elif [ $refinement -eq 4096 ]; then
                    amr_max_level=6
                elif [ $refinement -eq 16384 ]; then
                    amr_max_level=7
                elif [ $refinement -eq 65536 ]; then
                    amr_max_level=8
                elif [ $refinement -eq 262144 ]; then
                    amr_max_level=9
                elif [ $refinement -eq 1048576 ]; then
                    amr_max_level=10
                elif [ $refinement -eq 4194304 ]; then
                    amr_max_level=11
                else
                    echo "Unknown refinement factor: "$refinement"; exiting."
                    exit
                fi

                if [ ! -z $dxnuc_r ]; then
                    castro_max_dxnuc_lev=$(log_base_4 $dxnuc_r)
                else
                    castro_max_dxnuc_lev=0
                fi

                if [ ! -z $tempgrad_r ]; then
                    max_tempgrad_rel_lev=$(log_base_4 $tempgrad_r)
                else
                    max_tempgrad_rel_lev=0
                fi

            fi

        else

            amr_n_cell="$ncell $ncell"
            geometry_prob_lo="$prob_lo $prob_lo"
            geometry_prob_hi="$prob_hi $prob_hi"
            castro_lo_bc="3 3"
            castro_hi_bc="2 2"

            if   [ $ncell -eq 64 ]; then

                amr_blocking_factor="16"
                amr_blocking_factor="4"
                amr_max_grid_size="16 16 16 16 32 128 128 256 256"

            elif [ $ncell -eq 128 ]; then

                amr_blocking_factor="16"
                amr_max_grid_size="16 16 32 32 64 64 128 128 256 256"

            elif [ $ncell -eq 256 ]; then

                amr_blocking_factor="32"
                amr_blocking_factor="4"
                amr_max_grid_size="64 128 256 512 1024 1024 1024 1024"

            elif [ $ncell -eq 512 ]; then

                amr_blocking_factor="32"
                amr_max_grid_size="64 64 64 64 128 128 256 256 512 512"

            elif [ $ncell -eq 1024 ]; then

                amr_blocking_factor="32"
                amr_max_grid_size="64 64 128 128 256 256 512 512 1024 1024"

            elif [ $ncell -eq 2048 ]; then

                amr_blocking_factor="32"
                amr_max_grid_size="64 64 128 128 256 256 512 512 1024 1024"

            elif [ $ncell -eq 4096 ]; then

                amr_blocking_factor="64"
                amr_max_grid_size="128 128 256 256 512 512 1024 1024 2048 2048"

            elif [ $ncell -eq 8192 ]; then

                amr_blocking_factor="64"
                amr_max_grid_size="128 128 256 256 512 512 1024 1024 2048 2048"

            elif [ $ncell -eq 16384 ]; then

                amr_blocking_factor="64"
                amr_max_grid_size="128 128 256 256 512 512 1024 1024 2048 2048"

            fi



            if [ ! -z $stellar_refinement ]; then

                if [ $stellar_refinement -eq 1 ]; then

                    max_stellar_tagging_level="0"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=0
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
                        elif [ $refinement -eq 131072 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 262144 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 524288 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 1048576 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 2097152 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 4194304 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 4"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                elif [ $stellar_refinement -eq 2 ]; then

                    max_stellar_tagging_level="1"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=1
                            amr_ref_ratio="2"
                        elif [ $refinement -eq 2 ]; then
                            amr_max_level=2
                            amr_ref_ratio="2 2"
                        elif [ $refinement -eq 4 ]; then
                            amr_max_level=2
                            amr_ref_ratio="2 4"
                        elif [ $refinement -eq 8 ]; then
                            amr_max_level=3
                            amr_ref_ratio="2 4 2"
                        elif [ $refinement -eq 16 ]; then
                            amr_max_level=3
                            amr_ref_ratio="2 4 4"
                        elif [ $refinement -eq 32 ]; then
                            amr_max_level=4
                            amr_ref_ratio="2 4 4 2"
                        elif [ $refinement -eq 64 ]; then
                            amr_max_level=4
                            amr_ref_ratio="2 4 4 4"
                        elif [ $refinement -eq 128 ]; then
                            amr_max_level=5
                            amr_ref_ratio="2 4 4 4 2"
                        elif [ $refinement -eq 256 ]; then
                            amr_max_level=5
                            amr_ref_ratio="2 4 4 4 4"
                        elif [ $refinement -eq 512 ]; then
                            amr_max_level=6
                            amr_ref_ratio="2 4 4 4 4 2"
                        elif [ $refinement -eq 1024 ]; then
                            amr_max_level=6
                            amr_ref_ratio="2 4 4 4 4 4"
                        elif [ $refinement -eq 2048 ]; then
                            amr_max_level=7
                            amr_ref_ratio="2 4 4 4 4 4 2"
                        elif [ $refinement -eq 4096 ]; then
                            amr_max_level=7
                            amr_ref_ratio="2 4 4 4 4 4 4"
                        elif [ $refinement -eq 8192 ]; then
                            amr_max_level=8
                            amr_ref_ratio="2 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 16384 ]; then
                            amr_max_level=8
                            amr_ref_ratio="2 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 32768 ]; then
                            amr_max_level=9
                            amr_ref_ratio="2 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 65536 ]; then
                            amr_max_level=9
                            amr_ref_ratio="2 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 131072 ]; then
                            amr_max_level=10
                            amr_ref_ratio="2 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 262144 ]; then
                            amr_max_level=10
                            amr_ref_ratio="2 4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 524288 ]; then
                            amr_max_level=11
                            amr_ref_ratio="2 4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 1048576 ]; then
                            amr_max_level=11
                            amr_ref_ratio="2 4 4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 2097152 ]; then
                            amr_max_level=12
                            amr_ref_ratio="2 4 4 4 4 4 4 4 4 4 4 2"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                elif [ $stellar_refinement -eq 4 ]; then

                    max_stellar_tagging_level="1"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=1
                            amr_ref_ratio="4"
                        elif [ $refinement -eq 2 ]; then
                            amr_max_level=2
                            amr_ref_ratio="4 2"
                        elif [ $refinement -eq 4 ]; then
                            amr_max_level=2
                            amr_ref_ratio="4 4"
                        elif [ $refinement -eq 8 ]; then
                            amr_max_level=3
                            amr_ref_ratio="4 4 2"
                        elif [ $refinement -eq 16 ]; then
                            amr_max_level=3
                            amr_ref_ratio="4 4 4"
                        elif [ $refinement -eq 32 ]; then
                            amr_max_level=4
                            amr_ref_ratio="4 4 4 2"
                        elif [ $refinement -eq 64 ]; then
                            amr_max_level=4
                            amr_ref_ratio="4 4 4 4"
                        elif [ $refinement -eq 128 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 4 4 2"
                        elif [ $refinement -eq 256 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 4 4 4"
                        elif [ $refinement -eq 512 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 4 2"
                        elif [ $refinement -eq 1024 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 4 4"
                        elif [ $refinement -eq 2048 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 4 2"
                        elif [ $refinement -eq 4096 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 4 4"
                        elif [ $refinement -eq 8192 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 16384 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 32768 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 65536 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 131072 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 262144 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 524288 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 1048576 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 4"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                elif [ $stellar_refinement -eq 8 ]; then

                    max_stellar_tagging_level="2"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=2
                            amr_ref_ratio="4 2"
                        elif [ $refinement -eq 2 ]; then
                            amr_max_level=3
                            amr_ref_ratio="4 2 2"
                        elif [ $refinement -eq 4 ]; then
                            amr_max_level=3
                            amr_ref_ratio="4 2 4"
                        elif [ $refinement -eq 8 ]; then
                            amr_max_level=4
                            amr_ref_ratio="4 2 4 2"
                        elif [ $refinement -eq 16 ]; then
                            amr_max_level=4
                            amr_ref_ratio="4 2 4 4"
                        elif [ $refinement -eq 32 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 2 4 4 2"
                        elif [ $refinement -eq 64 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 2 4 4 4"
                        elif [ $refinement -eq 128 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 2 4 4 4 2"
                        elif [ $refinement -eq 256 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 2 4 4 4 4"
                        elif [ $refinement -eq 512 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 2 4 4 4 4 2"
                        elif [ $refinement -eq 1024 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 2 4 4 4 4 4"
                        elif [ $refinement -eq 2048 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 2 4 4 4 4 4 2"
                        elif [ $refinement -eq 4096 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 2 4 4 4 4 4 4"
                        elif [ $refinement -eq 8192 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 2 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 16384 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 2 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 32768 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 2 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 65536 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 2 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 131072 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 2 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 262144 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 2 4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 524288 ]; then
                            amr_max_level=12
                            amr_ref_ratio="4 2 4 4 4 4 4 4 4 4 4 2"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                elif [ $stellar_refinement -eq 16 ]; then

                    max_stellar_tagging_level="2"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=2
                            amr_ref_ratio="4 4"
                        elif [ $refinement -eq 2 ]; then
                            amr_max_level=3
                            amr_ref_ratio="4 4 2"
                        elif [ $refinement -eq 4 ]; then
                            amr_max_level=3
                            amr_ref_ratio="4 4 4"
                        elif [ $refinement -eq 8 ]; then
                            amr_max_level=4
                            amr_ref_ratio="4 4 4 2"
                        elif [ $refinement -eq 16 ]; then
                            amr_max_level=4
                            amr_ref_ratio="4 4 4 4"
                        elif [ $refinement -eq 32 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 4 4 2"
                        elif [ $refinement -eq 64 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 4 4 4"
                        elif [ $refinement -eq 128 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 4 2"
                        elif [ $refinement -eq 256 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 4 4"
                        elif [ $refinement -eq 512 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 4 2"
                        elif [ $refinement -eq 1024 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 4 4"
                        elif [ $refinement -eq 2048 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 4096 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 8192 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 16384 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 32768 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 65536 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 131072 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 262144 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 4"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                elif [ $stellar_refinement -eq 32 ]; then

                    max_stellar_tagging_level="3"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=3
                            amr_ref_ratio="4 4 2"
                        elif [ $refinement -eq 2 ]; then
                            amr_max_level=4
                            amr_ref_ratio="4 4 2 2"
                        elif [ $refinement -eq 4 ]; then
                            amr_max_level=4
                            amr_ref_ratio="4 4 2 4"
                        elif [ $refinement -eq 8 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 2 4 2"
                        elif [ $refinement -eq 16 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 2 4 4"
                        elif [ $refinement -eq 32 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 2 4 4 2"
                        elif [ $refinement -eq 64 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 2 4 4 4"
                        elif [ $refinement -eq 128 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 2 4 4 4 2"
                        elif [ $refinement -eq 256 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 2 4 4 4 4"
                        elif [ $refinement -eq 512 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 2 4 4 4 4 2"
                        elif [ $refinement -eq 1024 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 2 4 4 4 4 4"
                        elif [ $refinement -eq 2048 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 2 4 4 4 4 4 2"
                        elif [ $refinement -eq 4096 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 2 4 4 4 4 4 4"
                        elif [ $refinement -eq 8192 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 2 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 16384 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 2 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 32768 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 2 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 65536 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 2 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 131072 ]; then
                            amr_max_level=12
                            amr_ref_ratio="4 4 2 4 4 4 4 4 4 4 4 2"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                elif [ $stellar_refinement -eq 64 ]; then

                    max_stellar_tagging_level="3"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=3
                            amr_ref_ratio="4 4 4"
                        elif [ $refinement -eq 2 ]; then
                            amr_max_level=4
                            amr_ref_ratio="4 4 4 2"
                        elif [ $refinement -eq 4 ]; then
                            amr_max_level=4
                            amr_ref_ratio="4 4 4 4"
                        elif [ $refinement -eq 8 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 4 4 2"
                        elif [ $refinement -eq 16 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 4 4 4"
                        elif [ $refinement -eq 32 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 4 2"
                        elif [ $refinement -eq 64 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 4 4"
                        elif [ $refinement -eq 128 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 4 2"
                        elif [ $refinement -eq 256 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 4 4"
                        elif [ $refinement -eq 512 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 1024 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 2048 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 4096 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 8192 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 16384 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 32768 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 65536 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 4"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                elif [ $stellar_refinement -eq 128 ]; then

                    max_stellar_tagging_level="4"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=4
                            amr_ref_ratio="4 4 4 2"
                        elif [ $refinement -eq 2 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 4 2 2"
                        elif [ $refinement -eq 4 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 4 2 4"
                        elif [ $refinement -eq 8 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 2 4 2"
                        elif [ $refinement -eq 16 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 2 4 4"
                        elif [ $refinement -eq 32 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 2 4 4 2"
                        elif [ $refinement -eq 64 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 2 4 4 4"
                        elif [ $refinement -eq 128 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 2 4 4 4 2"
                        elif [ $refinement -eq 256 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 2 4 4 4 4"
                        elif [ $refinement -eq 512 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 2 4 4 4 4 2"
                        elif [ $refinement -eq 1024 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 2 4 4 4 4 4"
                        elif [ $refinement -eq 2048 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 2 4 4 4 4 4 2"
                        elif [ $refinement -eq 4096 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 2 4 4 4 4 4 4"
                        elif [ $refinement -eq 8192 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 2 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 16384 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 2 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 32768 ]; then
                            amr_max_level=12
                            amr_ref_ratio="4 4 4 2 4 4 4 4 4 4 4 2"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                elif [ $stellar_refinement -eq 256 ]; then

                    max_stellar_tagging_level="4"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=4
                            amr_ref_ratio="4 4 4 4"
                        elif [ $refinement -eq 2 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 4 4 2"
                        elif [ $refinement -eq 4 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 4 4 4"
                        elif [ $refinement -eq 8 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 4 2"
                        elif [ $refinement -eq 16 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 4 4"
                        elif [ $refinement -eq 32 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 4 2"
                        elif [ $refinement -eq 64 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 4 4"
                        elif [ $refinement -eq 128 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 256 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 512 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 1024 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 2048 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 4096 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 8192 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 16384 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 4"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                elif [ $stellar_refinement -eq 512 ]; then

                    max_stellar_tagging_level="5"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 4 4 2"
                        elif [ $refinement -eq 2 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 2 2"
                        elif [ $refinement -eq 4 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 2 4"
                        elif [ $refinement -eq 8 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 2 4 2"
                        elif [ $refinement -eq 16 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 2 4 4"
                        elif [ $refinement -eq 32 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 2 4 4 2"
                        elif [ $refinement -eq 64 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 2 4 4 4"
                        elif [ $refinement -eq 128 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 2 4 4 4 2"
                        elif [ $refinement -eq 256 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 2 4 4 4 4"
                        elif [ $refinement -eq 512 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 2 4 4 4 4 2"
                        elif [ $refinement -eq 1024 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 2 4 4 4 4 4"
                        elif [ $refinement -eq 2048 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 2 4 4 4 4 4 2"
                        elif [ $refinement -eq 4096 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 2 4 4 4 4 4 4"
                        elif [ $refinement -eq 8192 ]; then
                            amr_max_level=12
                            amr_ref_ratio="4 4 4 4 2 4 4 4 4 4 4 2"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                elif [ $stellar_refinement -eq 1024 ]; then

                    max_stellar_tagging_level="5"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=5
                            amr_ref_ratio="4 4 4 4 4"
                        elif [ $refinement -eq 2 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 4 2"
                        elif [ $refinement -eq 4 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 4 4"
                        elif [ $refinement -eq 8 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 4 2"
                        elif [ $refinement -eq 16 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 4 4"
                        elif [ $refinement -eq 32 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 64 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 128 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 256 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 512 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 1024 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 2048 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 4096 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 4"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                elif [ $stellar_refinement -eq 2048 ]; then

                    max_stellar_tagging_level="6"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 4 2"
                        elif [ $refinement -eq 2 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 2 2"
                        elif [ $refinement -eq 4 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 2 4"
                        elif [ $refinement -eq 8 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 2 4 2"
                        elif [ $refinement -eq 16 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 2 4 4"
                        elif [ $refinement -eq 32 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 2 4 4 2"
                        elif [ $refinement -eq 64 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 2 4 4 4"
                        elif [ $refinement -eq 128 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 2 4 4 4 2"
                        elif [ $refinement -eq 256 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 2 4 4 4 4"
                        elif [ $refinement -eq 512 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 2 4 4 4 4 2"
                        elif [ $refinement -eq 1024 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 2 4 4 4 4 4"
                        elif [ $refinement -eq 2048 ]; then
                            amr_max_level=12
                            amr_ref_ratio="4 4 4 4 4 2 4 4 4 4 4 2"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                elif [ $stellar_refinement -eq 4096 ]; then

                    max_stellar_tagging_level="6"

                    if [ ! -z $refinement ]; then

                        if [ $refinement -eq 1 ]; then
                            amr_max_level=6
                            amr_ref_ratio="4 4 4 4 4 4"
                        elif [ $refinement -eq 2 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 4 2"
                        elif [ $refinement -eq 4 ]; then
                            amr_max_level=7
                            amr_ref_ratio="4 4 4 4 4 4 4"
                        elif [ $refinement -eq 8 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 16 ]; then
                            amr_max_level=8
                            amr_ref_ratio="4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 32 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 64 ]; then
                            amr_max_level=9
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 128 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 256 ]; then
                            amr_max_level=10
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4"
                        elif [ $refinement -eq 512 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 2"
                        elif [ $refinement -eq 1024 ]; then
                            amr_max_level=11
                            amr_ref_ratio="4 4 4 4 4 4 4 4 4 4 4"
                        else
                            echo "Unknown refinement factor: "$refinement"; exiting."
                            exit
                        fi

                    fi

                fi

                # Allow stellar refinement all the way up to the max level.

                if [ ! -z $full_stellar_refinement ]; then
                    if [ $full_stellar_refinement -eq 1 ]; then
                        max_stellar_tagging_level=$amr_max_level
                    fi
                fi

            fi

        fi

    else

        echoerr "This machine is not set up for this job."

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

                    if [ $DIM -eq 2 ]; then

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

                    fi

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

DIM=1

problem_dir=$CASTRO_HOME/Exec/science/Detonation

use_first_castro_ex="1"

# Use the iso7 network for all the 1D tests.

NETWORK_DIR="iso7"

# Get the right inputs and probin files.

inputs="inputs-collision"
probin="probin-collision"

# Set plotfile details.

amr_plot_per="0.1"
amr_plot_vars="ALL"
amr_derive_plot_vars="ALL"

amr_small_plot_per="0.01"
amr_small_plot_vars="density Temp rho_e rho_c12 rho_o16 rho_si28 rho_ni56 enuc"
amr_derive_small_plot_vars="pressure soundspeed x_velocity y_velocity t_sound_t_enuc"

# Ensure that plotting intervals are hit exactly.
# This helps in resolving the detonation point.

castro_plot_per_is_exact="1"
castro_small_plot_per_is_exact="1"

# Checkpoint save rate.

amr_check_per="0.1"

# Checkpoints should not require plotfiles.

amr_write_plotfile_with_checkpoint="0"

# Set the maximum number of subcycles in a retry.

castro_max_subcycles="128"

# Initial timestep shortening factor.

castro_init_shrink="0.1"

# Set the interval for doing diagnostic global sums.

castro_sum_interval="0"

# Set verbosity.

castro_v="0"
gravity_v="0"
amr_v="1"

# Enable reactions and gravity.

castro_do_react="1"
castro_do_grav="1"

# Do not call the EOS in the burner.

call_eos_in_rhs="F"

# The timesteps can get quite small if you're fully 
# resolving the burning, so allow for this.

castro_dt_cutoff="1.0e-20"

# Allow the timestep to change by up to 5% per advance.

castro_change_max="1.05"

# Some defaults.

castro_small_temp="1.0e7"

spec_tol_default="1.d-6"
temp_tol_default="1.d-6"
enuc_tol_default="1.d-6"

rtol_spec=$spec_tol_default
atol_spec=$spec_tol_default
rtol_temp=$temp_tol_default
atol_temp=$temp_tol_default
rtol_enuc=$enuc_tol_default
atol_enuc=$enuc_tol_default

# The flag we will use to determine whether to run the job.

to_run=1

results_dir="results/1D"

castro_output_at_completion="1"

castro_use_stopping_criterion="1"
castro_T_stopping_criterion="4.0e9"
castro_ts_te_stopping_criterion="1.0e200"

amr_subcycling_mode="None"

stop_time="3.5"
max_step="10000000"

geometry_prob_lo="0.0"
geometry_prob_hi="1.6384e9"

ncell_list="64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144"

burning_mode_list="self-heat suppressed"

dxnuc_list="1.0e-1"
tempgrad_rel_list="0.5"

amr_n_error_buf="2"

gravity_const_grav="-1.1e8"
T_l="1.0d7"
T_r="1.0d7"
cfrac="0.5"
ofrac="0.45"
dens="5.0d6"
vel="2.0d8"

for burning_mode_str in $burning_mode_list
do

    if [ $burning_mode_str == "self-heat" ]; then
        burning_mode="1"
    elif [ $burning_mode_str == "suppressed" ]; then
        burning_mode="3"
    fi

    # Sanity check: X(C) + X(O) <= 1.0.

    if [ $(echo "if ($cfrac + $ofrac > 1.0) 1 else 0" | bc -l) -eq 1 ]; then
        continue
    fi

    for ncell in $ncell_list
    do

        # Low resolution runs tend to skip right through the burning
        # and have a very rough time with it as a result. Severely
        # limit the CFL number for these runs to compensate.

        if [ $ncell -le 128 ]; then
            castro_cfl="0.01"
        else
            castro_cfl="0.5"
        fi

        tempgrad_r_list="1 4 16 64"
        dxnuc_r_list="1 4 16 64"

        # Only do refinement for certain ncell values.

        if [ $ncell -ne 4096 ]; then
            tempgrad_r_list="1"
            dxnuc_r_list="1"
        fi

        # Only do refinement for the self-heating burns.

        if [ $burning_mode -ne 1 ]; then
            tempgrad_r_list="1"
            dxnuc_r_list="1"
        fi

        for tempgrad_rel in $tempgrad_rel_list
        do

            for dxnuc in $dxnuc_list
            do

                castro_dxnuc=$dxnuc

                for tempgrad_r in $tempgrad_r_list
                do

                    for dxnuc_r in $dxnuc_r_list
                    do

                        # Only do the suppressed burn for certain parameters.

                        if [ $burning_mode -eq 3 ]; then

                            if [ $ncell -gt 8192 ]; then
                                continue
                            fi

                        fi

                        # We want either dxnuc_r > 1 or tempgrad_r > 1 but not both, since the idea
                        # is to test them both as refinement criteria.

                        if [ $dxnuc_r -gt 1 ] && [ $tempgrad_r -gt 1 ]; then
                            continue
                        fi

                        refinement=1

                        if [ $dxnuc_r -gt $refinement ]; then
                            refinement=$dxnuc_r
                        fi

                        if [ $tempgrad_r -gt $refinement ]; then
                            refinement=$tempgrad_r
                        fi

                        # For the suppressed burns, the goal is to explicitly compare
                        # the simulation at a given time with the normal self-heating
                        # burn. So instead of using the same stopping criterion, we
                        # set the stopping time for the suppressed  burn to be the
                        # same as the stopping time for the self-heating burn, which
                        # requires the latter to have already finished.

                        dir_end=v$vel/c$cfrac/o$ofrac/n$ncell/tempgrad$tempgrad_rel/dxnuc$dxnuc/tempgrad_r$tempgrad_r/dxnuc_r$dxnuc_r

                        if [ $burning_mode -eq 3 ]; then

                            dir=$results_dir/self-heat/$dir_end

                            if [ $(is_dir_done) -eq 1 ]; then
                                checkpoint=$(get_last_checkpoint $dir)
                                stop_time=$(awk 'NR==3' $dir/$checkpoint/Header)
                                castro_T_stopping_criterion="1.0e200"

                                # For one particular run, we want to demonstrate
                                # what happens in the suppressed burn if you let
                                # it continue to develop. So what we'll do is wait
                                # until the job has reached the same stop time
                                # as the self-heating burn, so we can make a fair
                                # comparison with a plotfile at the same simulation
                                # time. Then we'll continue the job.

                                if [ $ncell -eq 64 ]; then

                                    dir=$results_dir/suppressed/$dir_end

                                    if [ -d $dir ]; then
                                        stop_time=$(get_inputs_var "stop_time")
                                        if [ $(is_dir_done) -eq 1 ] && [ "$stop_time" != "3.5" ]; then
                                            stop_time=3.5
                                            replace_inputs_var "stop_time"
                                            rm -f $dir/jobIsDone
                                        fi
                                    fi

                                fi
                            else
                                stop_time=3.5
                                castro_T_stopping_criterion="4.0e9"
                                continue
                            fi

                        fi

                        dir=$results_dir/$burning_mode_str/$dir_end
                        set_run_opts

                        if [ $to_run -eq 1 ]; then
                            run
                        fi

                        to_run=1

                    done # dxnuc_r

                done # tempgrad_r

            done # dxnuc

        done # tempgrad

    done # ncell

done # burning mode
