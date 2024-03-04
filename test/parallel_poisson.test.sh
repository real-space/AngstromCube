#!/usr/bin/env bash

. ../source_this_on_JUSUF.sh
exe=../build/src/green

out=parallel_poisson.out
rm -f $out
echo -n "# " > $out
date >> $out

for s in 2; do ### s=2:double, s=4: float
  for nprocs in {1..65}; do
    echo >> $out
    echo "# test parallel_poisson with $nprocs MPI processes" >> $out
    srun -n $nprocs \
    $exe  -t parallel_poisson \
            +parallel_poisson.select.test=$s \
            +parallel_poisson.plot.compressed=1e-6 \
            +verbosity=8 \
            >> $out
  done
done
