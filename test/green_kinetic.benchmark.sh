#!/usr/bin/env bash

exe=./green

## How to generate a test potential Hmt.xml for the Green function module?
##        +hamiltonian.export.format=xml +hamiltonian.export=-1
## -1 --> die after writing the file Hmt.xml

## e.g. ./a43 -t self_consistency +hamiltonian.export.format=xml +hamiltonian.export=-1 +bands.per.atom=0 +geometry.file=Cu40Zr22.xyz

for action in 412 812 422 822; do
  echo "#  action= $action (400:float, 800:double; 10:Noco=1, 20:Noco=2; 1:real, 2:complex)"
# for rtrunc in `seq 0 0.1 19`; do
#### To test a truncation radius > 64, increase the potential extend beyond 128 or change to Vacuum_Boundary (not well tested)
  out=green_kinetic.benchmark.$action.out
  rm -f $out
  touch $out

for rtrunc in {15..15}; do
  echo -n "rtrunc=$rtrunc "
  #   (cd ../src/ && make -j) && \
  #   (cd ../green/ && make -j) && \
  $exe -test green_function \
        +verbosity=6 \
        +green_function.benchmark.noco=2 \
        +green_function.benchmark.action=$action \
        +green_function.benchmark.iterations=-9 \
        +green_function.source.cube=1 \
        +green_function.truncation.radius="$rtrunc Bohr" \
        +green_function.potential.exchange=0 \
        +green_function.hamiltonian.file=Hmt.empty-iso.xml \
        +green_input.empty.potential=1 \
        "$@" \
       >> $out

done # rtrunc
done # action
