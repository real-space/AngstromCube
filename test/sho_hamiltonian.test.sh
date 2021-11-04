#!/usr/bin/env bash

exe=../src/a43

out=sho_hamiltonian.out
rm -f $out; echo -n "# " > $out; date >> $out

for asy in 1; do

  echo "# " >> $out
  echo "# " >> $out
  echo "# Assume that atoms.xyz and vtot.dat belong to each other!" >> $out
  echo "# " >> $out
  

(cd ../src/ && make -j) && \
  $exe +verbosity=10 \
    -test sho_hamiltonian \
    +sho_hamiltonian.test.numax=5 \
    +sho_hamiltonian.test.sigma=1.0 \
    +sho_hamiltonian.test.sigma.asymmetry=$asy \
    +hamiltonian.test.kpoints=17 \
    +hamiltonian.scale.kinetic=1 \
    +hamiltonian.scale.potential=1 \
    +hamiltonian.scale.nonlocal.s=1 \
    +hamiltonian.scale.nonlocal.h=1 \
    +sho_hamiltonian.reduce.centers=q \
    +dense_solver.test.overlap.eigvals=1 \
    +dense_solver.test.green.function=0 \
    +output.energy.unit=eV \
   >> $out

   ./spectrum.sh $out

done
grep 'deviation' $out
