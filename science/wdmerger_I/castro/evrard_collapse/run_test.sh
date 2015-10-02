source $WDMERGER_HOME/job_scripts/run_utils.sh

TEST_DIR=$CASTRO_DIR/Exec/evrard_collapse

cp $TEST_DIR/inputs $compile_dir/
cp $TEST_DIR/probin $compile_dir/

# Loop over the gravity source options

for gs in 1 2 3 4
do
  dir=$results_dir/gs$gs

  castro_grav_source_type=$gs

  if [ $MACHINE == "BLUE_WATERS" ]; then

    if [ $ncell -eq 64 ]; then
      nprocs=32
      walltime=1:00:00
    elif [ $ncell -eq 128 ]; then
      nprocs=128
      walltime=2:00:00
    fi

  fi

  run $dir $nprocs $walltime
done
