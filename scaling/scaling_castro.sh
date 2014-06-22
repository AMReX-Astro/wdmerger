mpi=`grep "MPI initialized with" $1 | awk '{print $4}'`
omp=`grep "OMP initialized with" $1 | awk '{print $4}'`


total=`grep "Coarse TimeStep time:" $1 | awk '{sum += $4; count += 1} END {print sum/count}'`

printf "# %10s  %10s  %10s \n" \
  "MPI" "threads" "total"

printf " %10d  %10d  %10f \n" \
   $mpi   $omp  $total

