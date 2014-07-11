mkdir $1/
mkdir $1/inputs
mkdir $1/model_files
mkdir $1/job_scripts

cp GNUmakefile $1/
cp Make.package $1/
cp *.f90 $1/
cp *.H $1/
cp *.cpp $1/
cp inputs/* $1/inputs/
cp model_files/* $1/model_files/
cp job_scripts/* $1/job_scripts/
cp copy_for_run.sh $1/

