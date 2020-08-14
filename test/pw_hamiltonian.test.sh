#!/usr/bin/env bash

exe=../src/a43

out=pw_hamiltonian.out
rm -f $out; echo -n "# " > $out; date >> $out

for asy in 1; do

(cd ../src/ && make -j) && \
  $exe +verbosity=7 \
    -test pw_hamiltonian. \
    +pw_hamiltonian.cutoff.energy=3.5 \
    +pw_hamiltonian.scale.kinetic=1 \
    +pw_hamiltonian.scale.potential=1 \
    +pw_hamiltonian.scale.nonlocal.h=0 \
    +pw_hamiltonian.scale.nonlocal.s=0 \
    +pw_hamiltonian.test.overlap.eigvals=1 \
    +pw_hamiltonian.test.kpoints=65 \
    >> $out

done
grep ' spectrum  ' $out | sed -e 's/#//g' -e 's/spectrum//g' -e 's/ Ha//g' -e 's/ eV//g' -e 's/\.\.\./  /g'
