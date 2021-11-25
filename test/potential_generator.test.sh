#!/usr/bin/env bash

exe=../src/a43
make -j -C ../src/

for sigma in {1..4}; do
  for ellmax in {0..9}; do
    project=$project_base.grid$spacing
    echo "# start calculation sigma= $sigma Bohr, ellmax= $ellmax"
    $exe -test potential_generator \
          +potential_generator.test.generalized.gaussian.sigma=$sigma \
          +potential_generator.test.generalized.gaussian.ellmax=$ellmax
  done
done
