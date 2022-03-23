#!/usr/bin/env bash

exe=./a43
project_base=green_function

## How to generate a test potential Hmt.xml for the Green function module?
##        +hamiltonian.export.format=xml +hamiltonian.export=-1
## -1 --> die after writing the file Hmt.xml

## e.g. ./a43 -t self_consistency +hamiltonian.export.format=xml +hamiltonian.export=-1 +bands.per.atom=0 +geometry.file=Cu40Zr22.xyz

# for rtrunc in `seq 0 0.1 19`; do
for rtrunc in 11; do
  echo -n "rtrunc=$rtrunc "
#   project=$project_base.r$rtrunc
#   (cd ../src/ && make -j) && \
  $exe -test green_function \
        +verbosity=17 \
        +green_input.empty.potential=1 \
        +green_function.source.cube=2 \
        +green_function.truncation.radius=$rtrunc \
        +green_function.scale.grid.spacing.x=0 \
        +green_function.scale.grid.spacing.y=0 \
        +green_function.benchmark.iterations=1 \
       > green_function.out

#         | grep 'target blocks per source block'
#         $1 > $project.out
done
#      +green_function.scale.grid.spacing.z=0 \
