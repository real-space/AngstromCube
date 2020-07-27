#!/usr/bin/env bash

exe=../src/a43

out=sho_hamiltonian.out
rm -f $out; echo -n "# " > $out; date >> $out

for asy in 1 1.1; do

(cd ../src/ && make -j) && \
  $exe +verbosity=8 \
    -test sho_hamiltonian. \
    +sho_hamiltonian.test.numax=1 \
    +sho_hamiltonian.test.sigma=1.5 \
    +sho_hamiltonian.test.sigma.asymmetry=$asy \
    +sho_hamiltonian.scale.kinetic=1 \
    +sho_hamiltonian.scale.potential=1 \
    +sho_hamiltonian.test.overlap.eigvals=1 \
   >> $out

done
grep 'deviation' $out
