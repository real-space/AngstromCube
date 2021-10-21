#!/usr/bin/env bash

exe=../src/a43
project_base=green_function

for rcut in `seq 0 0.1 9`; do
  echo -n "rcut=$rcut "
#   project=$project_base.r$rcut
#   (cd ../src/ && make -j) && \
  $exe -test green_function \
        +output.energy.unit=eV \
        +output.length.unit=Bohr \
        +verbosity=7 \
        +green.function.benchmark.iterations=-1 \
        +green.function.truncation.radius=$rcut \
        | grep 'target blocks per source block'
#         $1 > $project.out
done
