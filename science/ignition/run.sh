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
                        nprocs="2048"
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

amr_plot_per="-0.01"
amr_plot_vars="ALL"
amr_derive_plot_vars="ALL"

amr_small_plot_per="-0.01"
amr_small_plot_vars="density Temp rho_e rho_c12 rho_o16 rho_si28 rho_ni56 enuc"
amr_derive_small_plot_vars="pressure soundspeed x_velocity y_velocity t_sound_t_enuc"

# Ensure that plotting intervals are hit exactly.
# This helps in resolving the detonation point.

castro_plot_per_is_exact="1"
castro_small_plot_per_is_exact="1"

# Checkpoint save rate.

amr_check_per="-1.0"

# Checkpoints should not require plotfiles.

amr_write_plotfile_with_checkpoint="0"

# Initial timestep shortening factor.

castro_init_shrink="0.1"

# Set the interval for doing diagnostic global sums.

castro_sum_interval="1"

# Enable reactions.

castro_do_react="1"

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

stop_time="20.0"
max_step="10000000"

geometry_prob_lo="0.0"
geometry_prob_hi="8.192e8"

ncell_list="64"

burning_mode_list="self-heat suppressed"

dxnuc_list="1.0e-1"
tempgrad_rel_list="0.5"

amr_n_error_buf="2"

gravity_const_grav="-1.1e8"
T_l="1.0d7"
T_r="1.0d7"
cfrac="0.5"
ofrac="0.5"
dens="5.0d6"
gravity_const_grav="-1.1e8"
vel="2.0d8"

for burning_mode_str in $burning_mode_list
do

    if [ $burning_mode_str == "self-heat" ]; then
        burning_mode="1"
    elif [ $burning_mode_str == "suppressed" ]; then
        burning_mode="3"
    fi

    if [ $burning_mode -eq 3 ]; then
        continue
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

        tempgrad_r_list="1"
        dxnuc_r_list="1"

        # Only allow dxnuc_r > tempgrad_r for certain ncell values.

        allow_extra_dxnuc=0

        if [ $ncell -eq 64 ] || [ $ncell -eq 256 ] || [ $ncell -eq 8192 ]; then
            allow_extra_dxnuc=1
        fi

        # Only allow tempgrad_r > 1 for certain ncell values.

        if [ $ncell -gt 8192 ]; then
            tempgrad_r_list="1"
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

                            if [ "$tempgrad_rel" != "0.5" ]; then
                                continue
                            fi

                        fi


                        # Skip runs where dxnuc_r < tempgrad_r; these add no real extra information.

                        if [ $dxnuc_r -lt $tempgrad_r ]; then
                            continue
                        fi

                        # Set the refinement variable so we know how many levels to add.
                        # We'll make it equal to dxnuc_r since this will always be at
                        # least as large as tempgrad_r.

                        refinement=$dxnuc_r

                        # Skip runs where dxnuc_r > tempgrad_r except in the cases we've explicitly permitted,
                        # since we are only intending to demonstrate the benefit from the extra dxnuc
                        # in a few cases.

                        if [ $dxnuc_r -gt $tempgrad_r ] && [ $allow_extra_dxnuc -eq 0 ]; then
                            continue
                        fi

                        # For the suppressed burns, we cannot use a temperature
                        # stopping criterion in the same way as for the self-heating
                        # burns, since the suppressed burn gives a qualitatively
                        # different result. To make the comparison between the two
                        # fair, we will set the stopping time for the suppressed
                        # burn to be the same as the stopping time for the self-heating
                        # burn, which requires the latter to have already finished.

                        dir_end=g$g/v$vel/c$cfrac/o$ofrac/n$ncell/tempgrad$tempgrad_rel/dxnuc$dxnuc/tempgrad_r$tempgrad_r/dxnuc_r$dxnuc_r

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
                                        if [ $(is_dir_done) -eq 1 ] && [ "$stop_time" != "3.0" ]; then
                                            stop_time=3.0
                                            replace_inputs_var "stop_time"
                                            rm -f $dir/jobIsDone
                                        fi
                                    fi

                                fi
                            else
                                stop_time=3.0
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





# Now do the 2D collision problem.

# Specify the problem directory.

problem_dir=$CASTRO_HOME/Exec/science/wdmerger

# Needed for the makefile: we want to compile in 2D for these tests.

DIM="2"

# Use the iso7 network for all the 2D tests.

NETWORK_DIR="iso7"

# Get the right inputs and probin files.

inputs="inputs_2d"
probin="probin"

# Variables we need to set up the collision.

problem="0"
collision_separation="2.0"
collision_impact_parameter="0.0"

# Disable refinement based on non-burning criteria.
# These can be re-enabled later for specific tests.

max_stellar_tagging_level="0"
max_temperature_tagging_level="0"
max_center_tagging_level="0"

# Make a full plotfile every second.

amr_plot_per="1.0"
amr_plot_vars="ALL"
amr_derive_plot_vars="ALL"

# Make small plotfiles rapidly.

amr_small_plot_per="0.01"
amr_small_plot_vars="density Temp rho_e rho_c12 rho_o16 rho_si28 rho_ni56 enuc"
amr_derive_small_plot_vars="pressure soundspeed x_velocity y_velocity t_sound_t_enuc"

castro_plot_per_is_exact="0"
castro_small_plot_per_is_exact="0"

castro_output_at_completion="1"

# Save checkpoints every second.

amr_check_per="1.0"

# Initial timestep shortening factor.

castro_init_shrink="0.1"

# Set the interval for doing diagnostic global sums.

castro_sum_interval="1"

# Enable reactions.

castro_do_react="1"

# Disable rotation.

castro_do_rotation="0"

# Disable source term update after refluxes.
# These are relatively expensive, but not needed
# for this problem, since we do not need to be
# very accurate with the gravity source.

castro_update_sources_after_reflux="0"

# Do not solve gravity on the fine levels.

gravity_max_solve_level="0"

# Do not regrid in the middle of a step.
# Regrid every timestep to partially compensate.

castro_use_post_step_regrid="0"
amr_regrid_int=1

# To save on cycles, only worry about negative
# energies, not small energies, since we have
# to do an EOS call to guard against them.

castro_allow_small_energy="1"
castro_allow_negative_energy="0"

# The collision papers in the literature all use an equal C/O ratio 
# by mass in the initial white dwarfs. We will do this too for 
# comparison purposes.

co_wd_c_frac="0.5d0"
co_wd_o_frac="0.5d0"

# Use the sponge to damp out noise in the ambient medium,
# which can be rather violent due to how much churn
# the colliding white dwarfs create.

castro_do_sponge="1"
sponge_lower_density="1.0d-1"
sponge_upper_density="1.0d0"

# The timesteps can get quite small if you're fully 
# resolving the burning, so allow for this.

castro_dt_cutoff="1.0e-12"

# Allow the timestep to change by up to 5% per advance.
# We don't need to be super strict on the timestep control
# for this problem, we have empirically found.

castro_change_max="1.05"

# Some defaults.

ncell_default="256"
dtnuc_e_default="1.e200"
dtnuc_X_default="1.e200"
dxnuc_default="1.0e200"
dxnuc_max_default="1.0e200"
mass_P_default="0.64"
mass_S_default="0.64"

ncell=$ncell_default
mass_P=$mass_P_default
mass_S=$mass_S_default
castro_dtnuc_e=$dtnuc_e_default
castro_dtnuc_X=$dtnuc_X_default
castro_dxnuc=$dxnuc_default
castro_dxnuc_max=$dxnuc_max_default
castro_react_T_min="1.0e8"
castro_react_rho_min="1.0e6"
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

burning_mode=1


# The flag we will use to determine whether to run the job.

to_run=1

results_dir="results/2D"

mass_P=0.64
mass_S=0.64

ncell_list=""

ncell_list="64"
stop_time="9.0"
prob_lo="0.0"
prob_hi="1.6384e9"
stellar_density_threshold_list="5.0d6"
burning_mode_list="self-heat suppressed"
castro_T_stopping_criterion="4.0e9"
castro_dxnuc="1.0e-2"

for ncell in $ncell_list
do

    # The tolerance on the gravity solve needs to
    # be loosened with resolution, we have found
    # experimentally. This is mostly a problem with
    # the cylindrical grid we are using.

    if   [ $ncell -ge 8192 ]; then
        gravity_abs_tol="1.e-8"
    elif [ $ncell -eq 4096 ]; then
        gravity_abs_tol="1.e-9"
    elif [ $ncell -eq 2048 ]; then
        gravity_abs_tol="5.e-9"
    else
        gravity_abs_tol="1.e-10"
    fi

    # Do a run with no refinement.

    refinement=1

    base_dir=$results_dir/n$ncell
    start_dir=$base_dir/uniform/start
    end_dir=$base_dir/uniform/finish

    start_done=0

    if [ -d $start_dir ]; then
        dir=$start_dir
        start_done=$(is_dir_done)
    fi

    if [ $start_done -eq 0 ]; then

        dir=$start_dir
        set_run_opts

        castro_ts_te_stopping_criterion="1.0e-2"

        if [ $to_run -eq 1 ]; then
            run
        fi

        unset castro_ts_te_stopping_criterion

    else

        dir=$end_dir
        set_run_opts

        copy_checkpoint

        if [ $to_run -eq 1 ]; then
            run
        fi

    fi



    # Now refine on the density.

    for stellar_density_threshold in $stellar_density_threshold_list
    do

        stellar_refinement_list=""

        if   [ $ncell -eq 64 ]; then
            stellar_refinement_list="1"
        elif [ $ncell -eq 128 ]; then
            stellar_refinement_list="1"
        elif [ $ncell -eq 256 ]; then
            stellar_refinement_list="1"
        elif [ $ncell -eq 512 ]; then
            stellar_refinement_list="1"
        elif [ $ncell -eq 1024 ]; then
            stellar_refinement_list="1"
        elif [ $ncell -eq 2048 ]; then
            stellar_refinement_list="1"
        elif [ $ncell -eq 4096 ]; then
            stellar_refinement_list="1"
        else
            stellar_refinement_list="1"
        fi

        for stellar_refinement in $stellar_refinement_list
        do

            cur_dir=$base_dir/stellar$stellar_density_threshold/r$stellar_refinement
            start_dir=$cur_dir/start
            end_dir=$cur_dir/finish

            # For the stellar_refinement == 1 case, create a symlink to the uniform run
            # and then skip the actual calculation. This is helpful later on for analysis.

            if [ $stellar_refinement -eq 1 ]; then

                mkdir -p $cur_dir

                if [ ! -L $start_dir ]; then
                    target_dir=$(pwd)/$base_dir/uniform/start
                    if [ -d $target_dir ]; then
                        ln -s $target_dir $start_dir
                    fi
                fi

                if [ ! -L $end_dir ]; then
                    target_dir=$(pwd)/$base_dir/uniform/finish
                    if [ -d $target_dir ]; then
                        ln -s $target_dir $end_dir
                    fi
                fi

            else

                refinement=1

                start_done=0

                if [ -d $start_dir ]; then
                    dir=$start_dir
                    start_done=$(is_dir_done $start_dir)
                fi

                if [ $start_done -eq 0 ]; then

                    dir=$start_dir
                    set_run_opts

                    castro_ts_te_stopping_criterion="1.0e-2"

                    if [ $to_run -eq 1 ]; then
                        run
                    fi

                    unset castro_ts_te_stopping_criterion

                else

                    dir=$end_dir
                    set_run_opts

                    copy_checkpoint

                    if [ $to_run -eq 1 ]; then
                        run
                    fi

                fi

            fi



            # Do runs with burning-based AMR up to a selected level.

            refinement_list=1

            if [ $ncell -eq 64 ]; then
                if [ $stellar_refinement -eq 64 ]; then
                    refinement_list="1"
                fi
            fi

            # Do not subcycle this test.

            amr_subcycling_mode="None"

            for refinement in $refinement_list
            do

                cur_dir=$base_dir/stellar$stellar_density_threshold/r$stellar_refinement/dxnuc/r$refinement
                start_dir=$cur_dir/start
                end_dir=$cur_dir/finish

                mkdir -p $cur_dir

                if [ ! -L $start_dir ]; then
                    target_dir=$(pwd)/$base_dir/stellar$stellar_density_threshold/r$stellar_refinement/start
                    if [ -d $target_dir ]; then
                        ln -s $target_dir $start_dir
                    fi
                fi

                if [ $refinement -eq 1 ]; then

                    if [ ! -L $end_dir ]; then
                        target_dir=$(pwd)/$base_dir/stellar$stellar_density_threshold/r$stellar_refinement/finish
                        if [ -d $target_dir ]; then
                            ln -s $target_dir $end_dir
                        fi
                    fi

                else

                    start_done=0

                    if [ -d $start_dir ]; then
                        dir=$start_dir
                        start_done=$(is_dir_done)
                    fi

                    if [ $start_done -ne 0 ]; then

                        dir=$end_dir
                        set_run_opts

                        copy_checkpoint

                        if [ $to_run -eq 1 ]; then
                            run
                        fi

                    fi

                fi

            done

            unset amr_subcycling_mode
            unset refinement

        done # stellar_refinement

        unset stellar_refinement

    done # stellar_density_threshold

done # ncell
