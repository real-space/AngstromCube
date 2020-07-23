#!/usr/bin/env bash

exe=../src/a43

out=sho_hamiltonian.out
rm -f $out; echo -n "# " > $out; date >> $out

for asy in 1; do

(cd ../src/ && make -j) && \
  $exe +verbosity=9 \
    -test sho_hamiltonian. \
    +sho_hamiltonian.test.numax=3 \
    +sho_hamiltonian.test.sigma=2.0 \
    +sho_hamiltonian.test.method=6 \
    +sho_hamiltonian.test.sigma.asymmetry=$asy \
    +sho_hamiltonian.test.show.potential.only=0 \
   >> $out

done
grep 'deviation' $out
