#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

## make sure A43's single_atom.cxx has been compiled with -D HAS_AUTO_CONF

## experiment: check the semicore feature with sodium Z=11 Na
$exe  +verbosity=7 \
    -t single_atom \
      +single_atom.test.Z=$1 \
      +single_atom.config=element_config \
      +single_atom.core.state.localization=-1 \
      +element_config.core.semicore=-4.5 \
      +element_config.semicore.valence=-1.5 \
      +element_Na.rcut=2.27 \
      +element_Na.sigma=0.69 \
      +element_Na.numax=1 \
      +control.show=1 \
    >  element_config.semicore.out

## We can use the energy criterion to separate core and valence states
## with an energy window of 0 for semicore states by using this:
#     +element_config.core.valence=-2.2 \

