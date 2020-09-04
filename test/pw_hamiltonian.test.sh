#!/usr/bin/env bash

exe=../src/a43

out=pw_hamiltonian.out
rm -f $out; echo -n "# " > $out; date >> $out

for asy in 1; do

(cd ../src/ && make -j) && \
  $exe +verbosity=7 \
    -test pw_hamiltonian. \
    +pw_hamiltonian.cutoff.energy=2.5 \
    +pw_hamiltonian.scale.kinetic=1 \
    +pw_hamiltonian.scale.potential=1 \
    +pw_hamiltonian.scale.nonlocal.h=0 \
    +pw_hamiltonian.scale.nonlocal.s=0 \
    +dense_solver.test.overlap.eigvals=1 \
    +pw_hamiltonian.test.kpoints=65 \
    >> $out

done
grep ' spectrum  ' $out | sed -e 's/#//g' -e 's/spectrum//g' -e 's/ Ha//g' -e 's/ eV//g' -e 's/\.\.\./  /g'
grep 'average number of plane waves' $out

# ecut avg_N_plane_waves for #cell 4.233418 4.233418 8.466836, 65 k-points in x-direction
# 0.5   13.723 +/- 1.705 min  11 max  16   diff  5
# 1.5   87.692 +/- 4.462 min  77 max  92   diff 15
# 2.5  192.646 +/- 3.900 min 183 max 198   diff 15
# 3.5  325.692 +/- 8.341 min 314 max 345   diff 31
# 4.5  467.815 +/- 9.567 min 456 max 483   diff 27
# 5.5  632.092 +/- 4.676 min 624 max 643   diff 19
# 6.5  812.431 +/- 9.721 min 790 max 839   diff 49

