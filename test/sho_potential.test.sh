#!/usr/bin/env bash

exe=../src/a43

(cd ../src/ && make -j) && \
  $exe +verbosity=7 \
    -test sho_potential. \
    +sho_potential.test.numax=2 \
    +sho_potential.test.sigma=2.0 \
    +sho_potential.test.method=7 \
    -v

