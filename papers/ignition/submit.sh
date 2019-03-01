submit_dir="arxiv"

# Build the PDF so that we have the bibilography

make clean
make

# Move all the source files to the submit directory

mkdir -p $submit_dir

cp -f paper.tex $submit_dir
cp -f paper.bbl $submit_dir
cp -f plots/*.eps $submit_dir
cp -f ../aasjournal.bst $submit_dir
cp -f ../aastex62.cls $submit_dir

# Get rid of all references to the plots subdirectory

sed -i "s/plots\///" $submit_dir/paper.tex

# Get rid of all references to the parent directory

sed -i "s/\.\.\///" $submit_dir/paper.tex

# Replace the bibliography command with a direct input
# of the bibliography file

sed -i "s/bibliography{refs}/input{paper.bbl}/" $submit_dir/paper.tex

tar -czvf $submit_dir.tar $submit_dir
