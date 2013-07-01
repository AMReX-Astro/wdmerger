export OMP_NUM_THREADS=2

Castro='Castro3d.Linux.g++.gfortran.MPI.OMP.ex'
inputs='inputs_3d'

if [ ! -d results/ ]; then
  mkdir results
fi

if [ ! -d results/true/ ]; then
  mkdir results/true
  echo "Now computing exact solution"
  if [ ! grep -Fq "gravity.direct_sum_bcs" inputs_3d ];
  then
    echo "" >> inputs_3d
    echo "gravity.direct_sum_bcs = 1" >> inputs_3d
  fi
  sed -i "/gravity.direct_sum_bcs/c gravity.direct_sum_bcs = 1" inputs_3d
  mpiexec -n 8 $Castro $inputs >> info.out
  sh results_copy true
fi

sed -i "/gravity.direct_sum_bcs/c gravity.direct_sum_bcs = 0" inputs_3d

for l in {0..30}
do
  if [ ! -d results/$l ]; then
    echo "Now doing l =" $l
    rm -rf 
    sed -i "/gravity.max_multipole_order/c gravity.max_multipole_order = $l" inputs_3d
    mpiexec -n 8 $Castro $inputs >> info.out
    sh results_copy $l
  fi
done

#sh compare_results
