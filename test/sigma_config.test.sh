#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

$exe +verbosity=17 -test sigma_config \
     +element_Lu="5s* 2 6s 2 5p 6 6p 2e-99 5d* 1 .1 | 2.4 sigma .6 Z 71.1 3pCore 3 2.5 " \
     "$@" > sigma_config.out

     ### should launch this warning:
     ### sigma_config.cxx:549 warn("PAW setup for Lu (Z=71.1) is charged with -0.5 electrons")
