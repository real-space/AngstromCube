#!/usr/bin/env bash

exe=../src/a43

out=sho_hamiltonian.out
rm -f $out; echo -n "# " > $out; date >> $out

for asy in 1; do

(cd ../src/ && make -j) && \
  $exe +verbosity=10 \
    -test sho_hamiltonian. \
    +sho_hamiltonian.test.numax=2 \
    +sho_hamiltonian.test.sigma=1.5 \
    +sho_hamiltonian.test.sigma.asymmetry=$asy \
    +sho_hamiltonian.scale.kinetic=1 \
    +sho_hamiltonian.scale.potential=1 \
    +sho_hamiltonian.reduce.centers=q \
    +dense_solver.test.overlap.eigvals=1 \
    +dense_solver.test.green.function=20 \
   >> $out

done
grep 'deviation' $out
