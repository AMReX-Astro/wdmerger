source $WDMERGER_HOME/job_scripts/run_utils.sh

# Set up problem-specific inputs

mass_P=0.90
mass_S=0.75

castro_do_hydro=0
castro_do_rotation=0
castro_do_react=0

max_step=0

fab_format="NATIVE"

# Now set up the runs

minell=0
maxell=20

# Loop over the multipole orders we want to examine

if [ $MACHINE == "BLUE_WATERS" ]; then
    nprocs="32"
    walltime="1:00:00"
elif [ $MACHINE == "LIRED" ]; then
    nprocs="24"
    walltime="1:00:00"
fi

for l in $(seq $minell $maxell)
do
    dir=$results_dir/$l
    gravity_direct_sum_bcs=0
    gravity_max_multipole_order=$l
    run
done

if [ $MACHINE == "BLUE_WATERS" ]; then
    nprocs="1024"
    walltime="1:00:00"
elif [ $MACHINE == "LIRED" ]; then
    nprocs="576"
    walltime="2:00:00"
fi

# Now do the 'exact' direct summation, for comparison purposes

dir=$results_dir/true
gravity_direct_sum_bcs=1

run
