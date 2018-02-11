#!/bin/bash

source $WDMERGER_HOME/job_scripts/run_utils.sh

function set_run_opts {

    if [ -z $ncell ]; then
	echoerr "ncell not set; exiting."
	exit
    fi

    to_run=1

    # The following assumes we are using Titan.

    if [ $MACHINE == "TITAN" ]; then

        queue="batch"

        threads_per_task=4

        if [ $DIM -eq 2 ]; then

          if   [ $ncell -eq 16384 ]; then
              nprocs="2048"
              walltime="6:00:00"
          elif [ $ncell -eq 8192 ]; then
              nprocs="2048"
              walltime="6:00:00"
          elif [ $ncell -eq 4096 ]; then
              nprocs="1024"
              walltime="2:00:00"

              if [ $stellar_refinement -eq 8 ]; then
                  nprocs="2048"
                  walltime="6:00:00"
              elif [ $stellar_refinement -eq 16 ]; then
                  nprocs="2048"
                  walltime="6:00:00"
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
                  elif [ $stellar_refinement -eq 32 ]; then
                      nprocs="256"
                      walltime="2:00:00"
                  fi
              fi
          else
              echoerr "Unknown number of cells per dimension."
          fi

        elif [ $DIM -eq 1 ]; then

            nprocs="16"
            walltime="2:00:00"

            if [ $ncell -eq 256 ]; then
                if [ ! -z $stellar_refinement ]; then
                    if [ $stellar_refinement -eq 2048 ]; then
                        nprocs="32"
                    elif [ $stellar_refinement -eq 4096 ]; then
                        nprocs="64"
                    fi
                fi
            fi

        fi

    else

        echoerr "This machine is not set up for this job."

    fi

    # Set up the geometry appropriately.

    if [ $DIM -eq 2 ]; then

        amr_n_cell="$(echo "$ncell / 2" | bc) $ncell"
        geometry_prob_lo="0.0 $prob_lo"
        geometry_prob_hi="$prob_hi $prob_hi"

    elif [ $DIM -eq 1 ]; then

        amr_n_cell="$ncell"
        geometry_prob_lo="$prob_lo"
        geometry_prob_hi="$prob_hi"

    fi

    if [ $ncell -eq 256 ]; then

        amr_blocking_factor="32"
	amr_max_grid_size="64 64 64 64 128 128 256 256 512 512"

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

# Needed for the makefile: we want to compile in 2D for these tests.

DIM="2"

# Use the aprox13 network for all following tests.

NETWORK_DIR="aprox13"

# Get the right inputs and probin files.

inputs="inputs_2d"
probin="probin"

# Variables we need to set up the collision.

problem="0"
collision_separation="4.0"
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

# Ensure that plotting intervals are hit exactly.
# This helps in resolving the detonation point.

castro_plot_per_is_exact="1"
castro_small_plot_per_is_exact="1"

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



# The flag we will use to determine whether to run the job.

to_run=1



mass_list="0.20 0.25 0.30 0.64"

for mass in $mass_list
do

    mass_P=$mass
    mass_S=$mass

    ncell_list=""

    if   [ $mass == "0.64" ]; then

        ncell_list="256 512 1024 2048 4096 8192"
        stop_time="9.0"
        castro_react_rho_min="1.0e6"
        prob_lo="-4e9"
        prob_hi="4e9"

    elif [ $mass == "0.30" ]; then

        ncell_list="256 512 1024"
        stop_time="34.0"
        castro_react_rho_min="1.0e5"
        prob_lo="-8e9"
        prob_hi="8e9"

    elif [ $mass == "0.25" ]; then

        ncell_list="256 512"
        stop_time="54.0"
        castro_react_rho_min="1.0e0"
        prob_lo="-8e9"
        prob_hi="8e9"

    elif [ $mass == "0.20" ]; then

        ncell_list="256 512"
        stop_time="73.0"
        castro_react_rho_min="1.0e0"
        prob_lo="-8e9"
        prob_hi="8e9"

    fi

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

       stellar_refinement_list=""

       if   [ $mass == "0.64" ]; then

           if   [ $ncell -eq 256 ]; then
               stellar_refinement_list="1"
           elif [ $ncell -eq 512 ]; then
               stellar_refinement_list="1"
           elif [ $ncell -eq 1024 ]; then
               stellar_refinement_list="1 4 16"
           elif [ $ncell -eq 2048 ]; then
               stellar_refinement_list="1 16"
           elif [ $ncell -eq 4096 ]; then
               stellar_refinement_list="1"
           else
               stellar_refinement_list="1"
           fi

       elif [ $mass == "0.30" ]; then

           if   [ $ncell -eq 256 ]; then
               stellar_refinement_list="1"
           elif [ $ncell -eq 512 ]; then
               stellar_refinement_list="1"
           elif [ $ncell -eq 1024 ]; then
               stellar_refinement_list="1 2 4"
           fi

       elif [ $mass == "0.25" ]; then

           if   [ $ncell -eq 256 ]; then
               stellar_refinement_list="1"
           elif [ $ncell -eq 512 ]; then
               stellar_refinement_list="1 2 4"
           fi

       elif [ $mass == "0.20" ]; then

           if   [ $ncell -eq 256 ]; then
               stellar_refinement_list="1"
           elif [ $ncell -eq 512 ]; then
               stellar_refinement_list="1"
           fi

       fi

       for stellar_refinement in $stellar_refinement_list
       do

           base_dir=$results_dir/collision_2D/mass_P_$mass_P/mass_S_$mass_S/n$ncell/r$stellar_refinement
           start_dir=$base_dir/start

           start_done="0"

           if [ -d $start_dir ]; then
               dir=$start_dir
               start_done=$(is_dir_done)
           fi

           if [ $start_done -ne 1 ]; then

               refinement=1

               dir=$start_dir
               set_run_opts

               # First, we need to do the initial run, up to the point
               # when burning starts.

               castro_ts_te_stopping_criterion="9.0e-3"

               if   [ $mass == "0.64" ]; then
                   castro_ts_te_stopping_criterion="1.0e-6"
               fi

               if [ $to_run -eq 1 ]; then
                   run
               fi

               unset castro_ts_te_stopping_criterion
               unset refinement

           else

               # At this point we should have our starting checkpoint in $start_dir.
               # Now we can do the runs by copying the checkpoint.

               burning_mode_list="self-heat suppressed"

               for burning_mode_str in $burning_mode_list
               do

                   if [ $burning_mode_str == "self-heat" ]; then

                       burning_mode="1"

                       rtol_spec="1.d-6"
                       atol_spec="1.d-6"
                       rtol_temp="1.d-6"
                       atol_temp="1.d-6"
                       rtol_enuc="1.d-6"
                       atol_enuc="1.d-6"

                   elif [ $burning_mode_str == "suppressed" ]; then

                       burning_mode="3"

                       # Use tighter tolerances for the suppressed burn;
                       # this seems to prevent integration failures.

                       rtol_spec="1.d-8"
                       atol_spec="1.d-8"
                       rtol_temp="1.d-8"
                       atol_temp="1.d-8"
                       rtol_enuc="1.d-8"
                       atol_enuc="1.d-8"

                   fi

                   # Only do the suppressed burn at certain resolutions.

                   if [ $burning_mode_str == "suppressed" ]; then
                       if [ $mass != "0.64" ]; then
                           continue
                       fi
                       if [ $ncell -gt 4096 ]; then
                           continue
                       fi
                       if [ $stellar_refinement -gt 1 ]; then
                           continue
                       fi
                   fi

                   # Complete the run with no special options.

                   refinement=1

                   dir=$base_dir/$burning_mode_str/finish
                   set_run_opts

                   copy_checkpoint

                   if [ $to_run -eq 1 ]; then
                       run
                   fi



                   # The above is the only run we want to do with the suppressed burn.

                   if [ $burning_mode_str == "suppressed" ]; then
                       continue
                   fi



                   # Test the effect of the burning timestep limiter parameter values.

                   dtnuc_e_list=""
                   dtnuc_X_list=""
                   dtnuc_eX_list=""

                   if [ $mass == "0.64" ]; then

                       if [ $ncell -eq 256 ] && [ $stellar_refinement -eq 1 ]; then

                          dtnuc_e_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 5.0e0 2.0e0 1.0e0 5.0e-1 2.0e-1 1.0e-1 5.0e-2 2.0e-2 1.0e-2"
                          dtnuc_X_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 5.0e0 2.0e0 1.0e0 5.0e-1 2.0e-1 1.0e-1"
                          dtnuc_eX_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 5.0e0 2.0e0 1.0e0 5.0e-1 2.0e-1 1.0e-1"

                      elif [ $ncell -eq 512 ] && [ $stellar_refinement -eq 1 ]; then

                          dtnuc_e_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0 5.0e-1 2.0e-1 1.0e-1"
                          dtnuc_X_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0 5.0e-1 2.0e-1 1.0e-1"
                          dtnuc_eX_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0 5.0e-1 2.0e-1 1.0e-1"

                      elif [ $ncell -eq 1024 ] && [ $stellar_refinement -eq 1 ]; then

                          dtnuc_e_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0"
                          dtnuc_X_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0"
                          dtnuc_eX_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0"

                      elif [ $ncell -eq 2048 ] && [ $stellar_refinement -eq 1 ]; then

                          dtnuc_e_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0"
                          dtnuc_X_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0"
                          dtnuc_eX_list="1.0e6 1.0e5 1.0e4 1.0e3 1.0e2 1.0e1 1.0e0"

                      fi

                   fi

                   for dtnuc in $dtnuc_e_list
                   do

                       castro_dtnuc_e=$dtnuc

                       # For these tests only, disable timestep limiting for the plotting,
                       # since the whole point is to test the effect of the timestep.
                       # We'll still get enough plotfiles to examine the evolution if
                       # we choose to do so.

                       castro_plot_per_is_exact="0"
                       castro_small_plot_per_is_exact="0"

                       dir=$base_dir/$burning_mode_str/burning_limiter_e/dt$dtnuc
                       set_run_opts

                       copy_checkpoint

                       if [ $to_run -eq 1 ]; then
                           run
                       fi

                       castro_plot_per_is_exact="1"
                       castro_small_plot_per_is_exact="1"

                   done

                   castro_dtnuc_e=$dtnuc_e_default



                   for dtnuc in $dtnuc_X_list
                   do

                       castro_dtnuc_X=$dtnuc

                       castro_plot_per_is_exact="0"
                       castro_small_plot_per_is_exact="0"

                       dir=$base_dir/$burning_mode_str/burning_limiter_X/dt$dtnuc
                       set_run_opts

                       copy_checkpoint

                       if [ $to_run -eq 1 ]; then
                           run
                       fi

                       castro_plot_per_is_exact="1"
                       castro_small_plot_per_is_exact="1"

                   done

                   castro_dtnuc_X=$dtnuc_X_default



                   for dtnuc in $dtnuc_eX_list
                   do

                       castro_dtnuc_e=$dtnuc
                       castro_dtnuc_X=$dtnuc

                       castro_plot_per_is_exact="0"
                       castro_small_plot_per_is_exact="0"

                       dir=$base_dir/$burning_mode_str/burning_limiter_eX/dt$dtnuc
                       set_run_opts

                       copy_checkpoint

                       if [ $to_run -eq 1 ]; then
                           run
                       fi

                       castro_plot_per_is_exact="1"
                       castro_small_plot_per_is_exact="1"

                   done

                   castro_dtnuc_e=$dtnuc_e_default
                   castro_dtnuc_X=$dtnuc_X_default



                   # Do runs with stellar refinement up to a selected level.

                   refinement_list=""

                   if [ $mass == "0.64" ]; then
                       if [ $ncell -eq 256 ]; then
                           if [ $stellar_refinement -eq 1 ]; then
                               refinement_list="2 4 8 16"
                           fi
                       fi
                   fi

                   for refinement in $refinement_list
                   do

                       full_stellar_refinement=1

                       dir=$base_dir/$burning_mode_str/stellar/r$refinement
                       set_run_opts

                       copy_checkpoint

                       if [ $to_run -eq 1 ]; then
                           run
                       fi

                       unset full_stellar_refinement

                   done



                   # Do the runs with temperature-based refinement.

                   temp_list="1.0d9"

                   for temperature_tagging_threshold in $temp_list
                   do

                       refinement_list=""

                       if [ $mass == "0.64" ]; then
                           if [ $ncell -eq 1024 ]; then
                               if [ $stellar_refinement -eq 16 ]; then
                                   refinement_list="4"
                               fi
                           fi
                       fi

                       for refinement in $refinement_list
                       do

                           dir=$base_dir/$burning_mode_str/temperature/r$refinement
                           set_run_opts

                           max_temperature_tagging_level=$amr_max_level

                           copy_checkpoint

                           if [ $to_run -eq 1 ]; then
                               run
                           fi

                           unset max_temperature_tagging_level

                       done

                   done

                   max_temperature_tagging_level="0"



                   # Do runs with burning-based AMR up to a selected level.

                   dxnuc_list="1.0e-2"

                   for castro_dxnuc in $dxnuc_list
                   do

                       refinement_list=""

                       if   [ $mass == "0.64" ]; then

                          if   [ $ncell -eq 256 ]; then
                              if [ $stellar_refinement -eq 1 ]; then
                                  refinement_list="2 4 8 16 32"
                              fi
                          elif [ $ncell -eq 512 ]; then
                              if [ $stellar_refinement -eq 1 ]; then
                                  refinement_list="2 4 8 16 32"
                              fi
                          elif [ $ncell -eq 1024 ]; then
                              if [ $stellar_refinement -eq 1 ]; then
                                  refinement_list="2 4 8 16 32 64 128"
                              elif [ $stellar_refinement -eq 4 ]; then
                                  refinement_list="2 4 8 16 32 64"
                              elif [ $stellar_refinement -eq 16 ]; then
                                  refinement_list="2 4 8 16 32"
                              fi
                          elif [ $ncell -eq 2048 ]; then
                              if [ $stellar_refinement -eq 1 ]; then
                                  refinement_list="2 4 8 16 32 64"
                              fi
                          fi

                       elif [ $mass == "0.30" ]; then

                           if   [ $ncell -eq 256 ]; then
                               if [ $stellar_refinement -eq 1 ]; then
                                   refinement_list="2 4 8 16 32 64"
                               fi
                           elif [ $ncell -eq 512 ]; then
                               if [ $stellar_refinement -eq 1 ]; then
                                   refinement_list="2 4 8 16"
                               fi
                           elif [ $ncell -eq 1024 ]; then
                               if   [ $stellar_refinement -eq 1 ]; then
                                   refinement_list="2 4 8 16"
                               elif [ $stellar_refinement -eq 2 ]; then
                                   refinement_list="2 4 8 16 32 64"
                               fi
                           fi

                       elif [ $mass == "0.25" ]; then

                           if   [ $ncell -eq 256 ]; then
                               if [ $stellar_refinement -eq 1 ]; then
                                   refinement_list="2 4 8 16"
                               fi
                           elif [ $ncell -eq 512 ]; then
                               if   [ $stellar_refinement -eq 1 ]; then
                                   refinement_list="2 4 8 16"
                               elif [ $stellar_refinement -eq 2 ]; then
                                   refinement_list="2 4 8 16"
                               fi
                           fi

                       fi

                       for refinement in $refinement_list
                       do

                           dir=$base_dir/$burning_mode_str/dxnuc/f$castro_dxnuc/r$refinement/
                           set_run_opts

                           copy_checkpoint

                           if [ $to_run -eq 1 ]; then
                               run
                           fi

                       done

                   done

                   castro_dxnuc=$dxnuc_default



                   # Do runs with static refinement in the center up to a selected level.

                   center_tagging_radius="2.0d8"

                   refinement_list=""

                   if [ $mass == "0.64" ]; then
                       if [ $ncell -eq 1024 ]; then
                           if [ $stellar_refinement -eq 16 ]; then
                               refinement_list="2"
                           fi
                       fi
                   fi

                   for refinement in $refinement_list
                   do

                       max_center_tagging_level=$amr_max_level

                       dir=$base_dir/$burning_mode_str/center/d$center_tagging_radius/r$refinement/
                       set_run_opts

                       copy_checkpoint

                       if [ $to_run -eq 1 ]; then
                           run
                       fi

                       unset max_center_tagging_level

                   done

                   unset center_tagging_radius

               done

               burning_mode=$burning_mode_default

               rtol_spec=$spec_tol_default
               atol_spec=$spec_tol_default
               rtol_temp=$temp_tol_default
               atol_temp=$temp_tol_default
               rtol_enuc=$enuc_tol_default
               atol_enuc=$enuc_tol_default

           fi

       done

   done

done








# Now set up the 1D problem

DIM="1"

inputs="inputs_3d"
probin="probin"

collision_separation="2.0"
collision_impact_parameter="0.0"

max_stellar_tagging_level="0"
max_temperature_tagging_level="0"
max_center_tagging_level="0"

castro_react_T_min="1.0e8"
castro_react_rho_min="1.0e0"

mass_list="0.50 0.60 0.70 0.80 0.90 1.00"

prob_lo="-2.56e9"
prob_hi="2.56e9"

# Disable gravity for these 1D tests, because
# it doesn't make much sense anyway in 1D Cartesian.

castro_do_grav="0"

# Disable subcycling, so that we can watch the evolution
# of a detonation, and so that timesteps don't get too long.

amr_subcycling_mode="None"

castro_output_at_completion="1"

for mass in $mass_list
do

    mass_P=$mass
    mass_S=$mass

    ncell_list="256"

    for ncell in $ncell_list
    do

       stellar_refinement_list=""

       if   [ $mass == "0.70" ]; then

           if   [ $ncell -eq 256 ]; then
               stellar_refinement_list="1 2 4 8 16 32 64 128 256 512"
           fi

       elif [ $mass == "1.00" ]; then

           if   [ $ncell -eq 256 ]; then
               stellar_refinement_list="1 2 4 8 16 32 64 128 256 512 1024 2048 4096"
           fi

       fi

       for stellar_refinement in $stellar_refinement_list
       do

           base_dir=$results_dir/collision_1D/mass_P_$mass_P/mass_S_$mass_S/n$ncell/r$stellar_refinement
           start_dir=$base_dir/start

           start_done="0"

           if [ $mass == "1.00" ]; then

               if [ $stellar_refinement -le 128 ]; then
                   stop_time="0.25"
               elif [ $stellar_refinement -le 256 ]; then
                   stop_time="0.245"
               elif [ $stellar_refinement -eq 512 ]; then
                   stop_time="0.26"
               elif [ $stellar_refinement -eq 1024 ]; then
                   stop_time="0.285"
               elif [ $stellar_refinement -eq 2048 ]; then
                   stop_time="0.295"
               elif [ $stellar_refinement -eq 4096 ]; then
                   stop_time="0.3035"
               fi

           elif [ $mass == "0.70" ]; then

               if   [ $stellar_refinement -le 4 ]; then
                   stop_time="0.92"
               elif [ $stellar_refinement -eq 8 ]; then
                   stop_time="0.89"
               elif [ $stellar_refinement -eq 16 ]; then
                   stop_time="0.85"
               elif [ $stellar_refinement -eq 32 ]; then
                   stop_time="0.82"
               elif [ $stellar_refinement -eq 64 ]; then
                   stop_time="0.795"
               elif [ $stellar_refinement -eq 128 ]; then
                   stop_time="0.79"
               elif [ $stellar_refinement -eq 256 ]; then
                   stop_time="0.84"
               elif [ $stellar_refinement -eq 512 ]; then
                   stop_time="0.80"
               fi

           elif [ $mass == "0.50" ]; then

               if [ $stellar_refinement -le 64 ]; then
                   stop_time="1.15"
               elif [ $stellar_refinement -le 256 ]; then
                   stop_time="1.35"
               else
                   stop_time="1.45"
               fi

           fi

           if [ -d $start_dir ]; then
               dir=$start_dir
               start_done=$(is_dir_done)
           fi

           if [ $start_done -ne 1 ]; then

               refinement="1"

               dir=$start_dir
               set_run_opts

               if [ $to_run -eq 1 ]; then
                   run
               fi

               unset refinement

           else

               stop_time=$(echo "$stop_time + 1.0" | bc -l)

               burning_mode_list="self-heat suppressed"

               for burning_mode_str in $burning_mode_list
               do

                   if [ $burning_mode_str == "self-heat" ]; then

                       burning_mode="1"

                       rtol_spec="1.d-6"
                       atol_spec="1.d-6"
                       rtol_temp="1.d-6"
                       atol_temp="1.d-6"
                       rtol_enuc="1.d-6"
                       atol_enuc="1.d-6"

                   elif [ $burning_mode_str == "suppressed" ]; then

                       burning_mode="3"

                       # Use tighter tolerances for the suppressed burn;
                       # this seems to prevent integration failures.

                       rtol_spec="1.d-8"
                       atol_spec="1.d-8"
                       rtol_temp="1.d-8"
                       atol_temp="1.d-8"
                       rtol_enuc="1.d-8"
                       atol_enuc="1.d-8"

                   fi

                   if [ $burning_mode_str == "suppressed" ]; then
                       continue
                   fi



                   # Complete the run with no special options.

                   refinement=1

                   dir=$base_dir/$burning_mode_str/finish
                   set_run_opts

                   copy_checkpoint

                   if [ $to_run -eq 1 ]; then
                       run
                   fi



                   # Do runs with burning-based AMR up to a selected level.

                   dxnuc_list="1.0e-1 5.0e-2 1.0e-2"

                   refinement_list=""

                   for castro_dxnuc in $dxnuc_list
                   do

                       if   [ $ncell -eq 256 ]; then

                           if   [ $stellar_refinement -eq 16 ]; then
                               refinement_list="131072"
                           elif [ $stellar_refinement -eq 32 ]; then
                               refinement_list="2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536"
                           elif [ $stellar_refinement -eq 64 ]; then
                               refinement_list="2 4 8 16 32768"
                           elif [ $stellar_refinement -eq 128 ]; then
                               refinement_list="16384"
                           elif [ $stellar_refinement -eq 256 ]; then
                               refinement_list="8192"
                           elif [ $stellar_refinement -eq 512 ]; then
                               refinement_list="4 16 64 256 1024 4096"
                           elif [ $stellar_refinement -eq 1024 ]; then
                               refinement_list="4 16 64 256 1024 2048"
                           elif [ $stellar_refinement -eq 2048 ]; then
                               refinement_list="4 16 64 256 1024"
                           elif [ $stellar_refinement -eq 4096 ]; then
                               refinement_list="512"
                           fi

                       fi

                       for refinement in $refinement_list
                       do

                           dir=$base_dir/$burning_mode_str/dxnuc/f$castro_dxnuc/r$refinement/
                           set_run_opts

                           copy_checkpoint

                           if [ $to_run -eq 1 ]; then
                               run
                           fi

                       done

                       castro_dxnuc=$dxnuc_default

                       to_run="1"

                   done

                   burning_mode=$burning_mode_default

                   rtol_spec=$spec_tol_default
                   atol_spec=$spec_tol_default
                   rtol_temp=$temp_tol_default
                   atol_temp=$temp_tol_default
                   rtol_enuc=$enuc_tol_default
                   atol_enuc=$enuc_tol_default

               done

           fi

       done

    done

done
