#!/usr/bin/env bash

exe=../src/a43
project_base=green_function

for rtrunc in `seq 0 0.1 19`; do
  echo -n "rtrunc=$rtrunc "
#   project=$project_base.r$rtrunc
#   (cd ../src/ && make -j) && \
  $exe -test green_function \
        +verbosity=7 \
        +green_function.benchmark.iterations=-1 \
        +green_function.truncation.radius=$rtrunc \
        | grep 'target blocks per source block'
#         $1 > $project.out
done
