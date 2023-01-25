#!/usr/bin/env bash

exe=./a43

## How to generate a test potential Hmt.xml for the Green function module?
##        +hamiltonian.export.format=xml +hamiltonian.export=-1
## -1 --> die after writing the file Hmt.xml

## e.g. ./a43 -t self_consistency +hamiltonian.export.format=xml +hamiltonian.export=-1 +bands.per.atom=0 +geometry.file=Cu40Zr22.xyz

# for rtrunc in `seq 0 0.1 19`; do
for rtrunc in {25..25}; do
  echo -n "rtrunc=$rtrunc "
#   (cd ../src/ && make -j) && \
  $exe -test green_function \
        +verbosity=6 \
        +green_input.empty.potential=1 \
        +green_function.benchmark.iterations=-9 \
        +green_function.source.cube=1 \
        +green_function.truncation.radius="$rtrunc Bohr" \
        +green_function.potential.exchange=0 \
        +green_function.hamiltonian.file=Hmt.empty-iso.xml \
        "$@" \
       > green_kinetic.benchmark.out

done
