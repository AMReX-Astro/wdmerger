mkdir $1/inputs
mkdir $1/model_files
mkdir $1/job_scripts
cp *.ex $1/
cp inputs/* $1/inputs/
cp model_files/* $1/model_files/
cp job_scripts/* $1/job_scripts/
cp $CASTRO_DIR/EOS/helmeos/helm_table.dat $1/

